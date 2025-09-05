theme_Publication <- function(base_size=16, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family = "")
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "top",
            legend.direction = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

#' Read FlowJo CSV export into a SummarizedExperiment
#' @param path character, path to FlowJo CSV
#' @param assay_name character, name for assay (default "exprs")
#' @param cells_as_columns logical, TRUE => channels x cells (default)
#' @return SummarizedExperiment
read_flowjo_csv_to_SE <- function(path, assay_name = "exprs", cells_as_columns = TRUE) {
  stopifnot(file.exists(path))
  lines <- readLines(path, warn = FALSE)
  
  # locate the first quoted, comma-separated column header line
  hdr_idx <- suppressWarnings(which(grepl('^"', lines) & grepl('","', lines))[1])
  if (is.na(hdr_idx)) {
    # fallback: first line with many commas that's not a simple KEY,VALUE pair
    comma_count <- lengths(regmatches(lines, gregexpr(",", lines)))
    cand <- which(comma_count > 5 & !grepl("^[A-Z$ _-]+,\\s*[^,]+$", lines))
    hdr_idx <- cand[1]
  }
  if (is.na(hdr_idx)) stop("Could not locate data header in: ", path)
  
  # parse key/value metadata above the data header
  kv_raw <- strsplit(lines[seq_len(hdr_idx - 1)], ",", fixed = TRUE)
  kv_raw <- kv_raw[sapply(kv_raw, length) >= 2]
  keys <- make.unique(vapply(kv_raw, `[`, "", 1), sep = "_")
  vals <- vapply(kv_raw, function(x) paste(x[-1], collapse = ","), "")
  suppressWarnings(vals_num <- as.numeric(vals))
  vals <- ifelse(!is.na(vals_num), vals_num, vals)
  file_meta <- as.list(structure(vals, names = keys))
  
  # read the data block
  df <- utils::read.csv(path, skip = hdr_idx - 1, check.names = FALSE)
  df[] <- lapply(df, function(x) {
    if (is.factor(x)) as.numeric(as.character(x)) else as.numeric(x)
  })
  
  mat <- as.matrix(df)
  if (cells_as_columns) {
    assay <- t(mat)                        # features x cells
    feature_names <- colnames(df)          # channels
    cell_names <- sprintf("cell_%d", seq_len(nrow(df)))
  } else {
    assay <- mat                           # cells x channels (not typical)
    feature_names <- rownames(df)
    cell_names <- colnames(df)
  }
  
  rownames(assay) <- feature_names
  colnames(assay) <- cell_names
  
  rowData <- S4Vectors::DataFrame(marker = feature_names, row.names = feature_names)
  colData <- S4Vectors::DataFrame(event_id = seq_len(ncol(assay)), row.names = cell_names)
  
  SummarizedExperiment::SummarizedExperiment(
    assays   = S4Vectors::SimpleList(setNames(list(assay), assay_name)),
    rowData  = rowData,
    colData  = colData,
    metadata = list(file = basename(path), header_index = hdr_idx, flowjo_meta = file_meta)
  )
}

annotate_antibodies <- function(se, map, field = "antibody") {
  stopifnot(inherits(se, "SummarizedExperiment"))
  rn <- rownames(se)
  ab <- setNames(rep(NA_character_, length(rn)), rn)
  hit <- intersect(names(map), rn)
  ab[hit] <- unname(map[hit])
  if (!field %in% colnames(rowData(se))) rowData(se)[[field]] <- NA_character_
  rowData(se)[[field]] <- ab
  metadata(se)$channel_to_antibody <- map
  se
}

# QC for Flow/Imaging cytometry exported to a SummarizedExperiment
# cells = columns, channels = rows

# --------------------------
# Example usage
# --------------------------
# qc <- qc_flow_se_suite(se_combined, assay_name = "exprs", condition_col = "condition",
#                        transform = "arcsinh", cofactor = 150)
# qc$plots$cells_per_condition
# qc$plots$channel_distributions
# qc$plots$total_signal
# qc$plots$pca
# qc$plots$correlations[["condA"]]   # heatmap for a given condition
# qc$summary                         # per-channel per-condition stats

# QC for cytometry-style SummarizedExperiment (channels x cells)
qc_flow_se_suite <- function(
    se,
    assay_name = NULL,
    condition_col = "condition",
    transform = c("arcsinh","none","log10"),
    cofactor = 150,
    cell_sample_per_condition = 10000,
    corr_max_cells = 50000,
    density_facet_ncol = 5,
    label_field = "antibody"
){
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(condition_col %in% colnames(SummarizedExperiment::colData(se)))
  
  # choose assay
  if (is.null(assay_name)) assay_name <- SummarizedExperiment::assayNames(se)[1]
  stopifnot(assay_name %in% SummarizedExperiment::assayNames(se))
  X <- SummarizedExperiment::assay(se, assay_name)   # channels x cells
  cond <- as.character(SummarizedExperiment::colData(se)[[condition_col]])
  stopifnot(ncol(X) == length(cond))
  
  # labels from rowData; disambiguate duplicates by appending channel id
  get_feature_labels <- function(se, field){
    ids  <- rownames(se)
    rdat <- SummarizedExperiment::rowData(se)
    labs <- if (field %in% colnames(rdat)) as.character(rdat[, field]) else ids
    labs[is.na(labs) | labs == ""] <- ids[is.na(labs) | labs == ""]
    dup <- duplicated(labs) | duplicated(labs, fromLast = TRUE)
    labs[dup] <- paste0(labs[dup], " [", ids[dup], "]")
    stats::setNames(labs, ids)
  }
  feat_lab <- get_feature_labels(se, label_field)
  
  # 1) cells per condition
  counts_df <- as.data.frame(table(condition = cond))
  p_counts <- ggplot2::ggplot(counts_df, ggplot2::aes(condition, Freq)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = Freq), vjust = -0.3, size = 3) +
    ggplot2::labs(x = "Condition", y = "Cells", title = "Cells per condition") +
    ggplot2::theme_bw(base_size = 12)
  
  # transform
  transform <- match.arg(transform)
  Xtr <- switch(transform,
                "arcsinh" = asinh(X / cofactor),
                "log10"   = { tmp <- X; tmp[tmp <= 0] <- NA_real_; log10(tmp) },
                "none"    = X
  )
  
  # sampling index (per condition) for density/PCA
  set.seed(1)
  idx_sample <- unlist(lapply(split(seq_along(cond), cond), function(ix){
    if (length(ix) > cell_sample_per_condition) sample(ix, cell_sample_per_condition) else ix
  }), use.names = FALSE)
  
  # 2) distributions by antibody label
  dens_mat <- t(Xtr[, idx_sample, drop = FALSE])       # cells x channels
  colnames(dens_mat) <- rownames(Xtr)                  # force channel ids as names
  dens_df <- data.frame(dens_mat, check.names = FALSE)
  dens_df[[condition_col]] <- cond[idx_sample]
  dens_long <- dens_df |>
    dplyr::mutate(.cell = dplyr::row_number()) |>
    tidyr::pivot_longer(
      cols = -c(dplyr::all_of(condition_col), .cell),
      names_to = "channel", values_to = "value"
    ) |>
    dplyr::mutate(feature_label = feat_lab[channel])
  
  p_dens <- ggplot2::ggplot(
    dens_long,
    ggplot2::aes(x = value, fill = .data[[condition_col]])
  ) +
    ggplot2::geom_density(alpha = 0.35, na.rm = TRUE) +
    ggplot2::facet_wrap(~ feature_label, scales = "free", ncol = density_facet_ncol) +
    ggplot2::labs(
      x = paste0(
        "Value (", transform,
        if (identical(transform, "arcsinh")) paste0(", cofactor=", cofactor) else "",
        ")"
      ),
      y = "Density", fill = "Condition",
      title = "Channel distributions by condition"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
  
  # 3a) summary stats per channel & condition (includes pct ≥ global 99th)
  q99 <- apply(Xtr, 1, function(v) stats::quantile(v, 0.99, na.rm = TRUE))
  q99_df <- tibble::tibble(channel = names(q99), q99 = unname(q99))
  
  long_values <- data.frame(t(Xtr), check.names = FALSE) |>
    dplyr::mutate(!!condition_col := cond) |>
    tidyr::pivot_longer(
      cols = -dplyr::all_of(condition_col),
      names_to = "channel", values_to = "value"
    ) |>
    dplyr::left_join(q99_df, by = "channel") |>
    dplyr::mutate(antibody = feat_lab[channel])
  
  summary_tbl <- long_values |>
    dplyr::group_by(.data[[condition_col]], antibody, channel) |>
    dplyr::summarise(
      n_cells   = sum(!is.na(value)),
      median    = stats::median(value, na.rm = TRUE),
      mad       = stats::mad(value, na.rm = TRUE, constant = 1),
      mean      = base::mean(value, na.rm = TRUE),
      sd        = stats::sd(value, na.rm = TRUE),
      p01       = stats::quantile(value, 0.01, na.rm = TRUE),
      p99       = stats::quantile(value, 0.99, na.rm = TRUE),
      pct_neg   = base::mean(value < 0, na.rm = TRUE),
      pct_zero  = base::mean(value == 0, na.rm = TRUE),
      pct_ge_q99 = base::mean(value >= q99, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 3b) per-cell total signal
  ts_df <- tibble::tibble(
    condition = cond,
    total_signal = colSums(Xtr, na.rm = TRUE)
  )
  p_total <- ggplot2::ggplot(ts_df, ggplot2::aes(x = condition, y = total_signal)) +
    ggplot2::geom_boxplot(outlier.size = 0.5) +
    ggplot2::labs(x = "Condition",
                  y = paste0("Sum across channels (", transform, ")"),
                  title = "Per-cell total signal") +
    ggplot2::theme_bw(base_size = 12)
  
  # 3c) channel–channel correlation (Spearman) per condition; axes use antibody labels
  cor_plots <- list()
  for (cc in unique(cond)) {
    ix <- which(cond == cc)
    if (length(ix) > corr_max_cells) ix <- sample(ix, corr_max_cells)
    C <- suppressWarnings(stats::cor(
      t(Xtr[, ix, drop = FALSE]),
      method = "spearman",
      use = "pairwise.complete.obs"
    ))
    C[is.na(C)] <- 0
    labs_axes <- unname(feat_lab[rownames(Xtr)])
    rownames(C) <- colnames(C) <- labs_axes
    C_df <- as.data.frame(as.table(C))
    names(C_df) <- c("Var1","Var2","value")
    cor_plots[[cc]] <- ggplot2::ggplot(C_df, ggplot2::aes(Var1, Var2, fill = value)) +
      ggplot2::geom_raster() +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = paste0("Channel Spearman correlation (", cc, ")"),
                    x = NULL, y = NULL, fill = "\u03C1") +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
  
  # 3d) PCA on sampled subset
  Xp <- t(scale(t(Xtr[, idx_sample, drop = FALSE]), center = TRUE, scale = TRUE))
  Xp[is.na(Xp)] <- 0
  pc <- stats::prcomp(t(Xp), center = FALSE, scale. = FALSE)
  var_expl <- pc$sdev^2 / sum(pc$sdev^2)
  p_pca <- tibble::tibble(
    PC1 = pc$x[,1], PC2 = pc$x[,2],
    condition = cond[idx_sample]
  ) |>
    ggplot2::ggplot(ggplot2::aes(PC1, PC2, color = condition)) +
    ggplot2::geom_point(alpha = 0.5, size = 0.7) +
    ggplot2::labs(
      title = sprintf("PCA (subset) — PC1 %.1f%%, PC2 %.1f%%",
                      100*var_expl[1], 100*var_expl[2])
    ) +
    ggplot2::theme_bw(base_size = 12)
  
  list(
    plots = list(
      cells_per_condition    = p_counts,
      channel_distributions  = p_dens,
      total_signal           = p_total,
      pca                    = p_pca,
      correlations           = cor_plots
    ),
    summary = summary_tbl,
    params = list(assay = assay_name, transform = transform, cofactor = cofactor,
                  label_field = label_field)
  )
}

get_antibody_distribution <- function(
    se,
    antibody = "H3",
    assay_name = NULL,
    condition_col = "condition",
    label_field = "antibody",
    transform = c("arcsinh","none","log10"),
    cofactor = 150,
    max_cells_per_condition = 20000,
    match = c("exact","ignore_case","substring"),
    seed = 1
){
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(condition_col %in% colnames(SummarizedExperiment::colData(se)))
  transform <- match.arg(transform)
  match <- match.arg(match)
  
  # assay + metadata
  if (is.null(assay_name)) assay_name <- SummarizedExperiment::assayNames(se)[1]
  X <- SummarizedExperiment::assay(se, assay_name)   # features (channels) x cells
  cond <- as.character(SummarizedExperiment::colData(se)[[condition_col]])
  ids  <- rownames(se)
  labs <- if (label_field %in% colnames(SummarizedExperiment::rowData(se)))
    as.character(SummarizedExperiment::rowData(se)[, label_field]) else ids
  labs[is.na(labs) | labs==""] <- ids[is.na(labs) | labs==""]
  
  # feature selection
  hit <- switch(match,
                exact        = which(labs == antibody),
                ignore_case  = which(tolower(labs) == tolower(antibody)),
                substring    = grep(antibody, labs, ignore.case = TRUE, fixed = FALSE)
  )
  if (length(hit) == 0) stop("No features match antibody label: ", antibody)
  sel_ids <- ids[hit]
  
  # transform
  Xtr <- switch(transform,
                arcsinh = asinh(X / cofactor),
                log10   = { tmp <- X; tmp[tmp <= 0] <- NA_real_; log10(tmp) },
                none    = X
  )
  
  # sample cells per condition (for plotting-scale data sizes)
  set.seed(seed)
  idx <- unlist(lapply(split(seq_along(cond), cond), function(ix){
    if (length(ix) > max_cells_per_condition) sample(ix, max_cells_per_condition) else ix
  }), use.names = FALSE)
  
  # tidy output
  df <- data.frame(t(Xtr[sel_ids, idx, drop = FALSE]), check.names = FALSE)
  df[[condition_col]] <- cond[idx]
  out <- df |>
    tidyr::pivot_longer(cols = tidyselect::all_of(sel_ids),
                        names_to = "channel", values_to = "value") |>
    dplyr::mutate(
      antibody = labs[match(channel, ids)],
      assay = assay_name,
      transform = transform,
      cofactor = ifelse(transform == "arcsinh", cofactor, NA_real_)
    ) |>
    dplyr::relocate(!!condition_col, antibody, channel, value, assay, transform, cofactor)
  
  tibble::as_tibble(out)
}

filter_cells_by_mad_window <- function(
    se,
    channel,
    assay_name = NULL,
    condition_col = "condition",
    label_field = "antibody",
    transform = c("arcsinh","none","log10"),
    cofactor = 150,
    k = 2,
    combine = c("first","any","all"),
    max_cells_per_condition = 20000,
    seed = 1,
    return_plot = TRUE
){
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(condition_col %in% colnames(SummarizedExperiment::colData(se)))
  transform <- match.arg(transform)
  combine <- match.arg(combine)
  
  if (is.null(assay_name)) assay_name <- SummarizedExperiment::assayNames(se)[1]
  X <- SummarizedExperiment::assay(se, assay_name)
  ids <- rownames(se)
  cond <- as.character(SummarizedExperiment::colData(se)[[condition_col]])
  
  rdat <- SummarizedExperiment::rowData(se)
  labs <- if (label_field %in% colnames(rdat)) as.character(rdat[, label_field]) else ids
  labs[is.na(labs) | labs==""] <- ids[is.na(labs) | labs==""]
  
  hit <- which(labs == channel)
  if (!length(hit)) hit <- which(tolower(labs) == tolower(channel))
  if (!length(hit)) hit <- which(ids == channel)
  if (!length(hit)) hit <- grep(channel, labs, ignore.case = TRUE)
  if (!length(hit)) stop("Channel/antibody not found: ", channel)
  sel_ids <- ids[hit]; sel_labs <- labs[hit]
  
  Xtr <- switch(transform,
                "arcsinh" = asinh(X / cofactor),
                "log10"   = { tmp <- X; tmp[tmp <= 0] <- NA_real_; log10(tmp) },
                "none"    = X
  )
  
  feature_windows <- lapply(seq_along(sel_ids), function(i){
    v <- as.numeric(Xtr[sel_ids[i], ])
    med <- stats::median(v, na.rm = TRUE)
    md  <- stats::mad(v, center = med, constant = 1, na.rm = TRUE)
    if (!is.finite(md) || md == 0) {
      iqr <- diff(stats::quantile(v, c(0.25, 0.75), na.rm = TRUE))
      alt <- if (is.finite(iqr) && iqr > 0) iqr/1.349 else stats::sd(v, na.rm = TRUE)
      if (!is.finite(alt) || alt == 0) alt <- 0
      md <- alt
    }
    list(id = sel_ids[i], label = sel_labs[i],
         median = med, mad = md, low = med - k*md, high = med + k*md)
  })
  
  masks <- lapply(seq_along(feature_windows), function(i){
    v <- as.numeric(Xtr[feature_windows[[i]]$id, ])
    (v >= feature_windows[[i]]$low) & (v <= feature_windows[[i]]$high)
  })
  keep_idx <- switch(combine,
                     "first" = masks[[1]],
                     "any"   = Reduce(`|`, masks),
                     "all"   = Reduce(`&`, masks)
  )
  keep_idx[is.na(keep_idx)] <- FALSE
  
  # per-condition % kept
  tab_total <- table(cond)
  tab_keep  <- tapply(keep_idx, cond, sum, na.rm = TRUE)
  by_cond <- tibble::tibble(
    condition = names(tab_total),
    n_total   = as.integer(unname(tab_total)),
    n_keep    = as.integer(unname(tab_keep[names(tab_total)])),
    pct_keep  = 100 * n_keep / n_total
  )
  # legend labels like: "treated (87.3% kept)"
  lab_map <- setNames(sprintf("%s (%.1f%% kept)", by_cond$condition, by_cond$pct_keep),
                      by_cond$condition)
  
  per_feature <- do.call(rbind, lapply(feature_windows, function(w)
    data.frame(feature_id = w$id, antibody = w$label,
               median = w$median, mad = w$mad, low = w$low, high = w$high,
               stringsAsFactors = FALSE)
  ))
  
  # Plot with % kept in legend (and subtitle)
  plt <- NULL
  if (return_plot) {
    set.seed(seed)
    idx <- unlist(lapply(split(seq_along(cond), cond), function(ix){
      if (length(ix) > max_cells_per_condition) sample(ix, max_cells_per_condition) else ix
    }), use.names = FALSE)
    
    subtitle_txt <- paste(by_cond$condition, sprintf("(%.1f%% kept)", by_cond$pct_keep),
                          collapse = "   ")
    
    if (length(sel_ids) == 1 || combine == "first") {
      v <- as.numeric(Xtr[sel_ids[1], idx]); cc <- cond[idx]
      w <- feature_windows[[1]]
      df <- dplyr::tibble(value = v,
                          condition = factor(cc, levels = by_cond$condition),
                          condition_lab = factor(lab_map[cc], levels = lab_map[by_cond$condition]))
      plt <- ggplot2::ggplot(df, ggplot2::aes(x = value, color = condition_lab, fill = condition_lab)) +
        ggplot2::geom_density(alpha = 0.25, na.rm = TRUE) +
        ggplot2::annotate("rect", ymin = -Inf, ymax = Inf,
                          xmin = w$low, xmax = w$high, alpha = 0.12) +
        ggplot2::geom_vline(xintercept = c(w$low, w$high), linetype = "dashed") +
        ggplot2::labs(
          title = paste0("Keep within \u00B1", k, " MAD — ", w$label, " (", w$id, ")"),
          subtitle = subtitle_txt,
          x = paste0("Signal (", transform,
                     if (identical(transform, "arcsinh")) paste0(", cofactor=", cofactor) else "", ")"),
          y = "Density", color = "Condition", fill = "Condition"
        ) + theme_bw(base_size = 12)
    } else {
      long <- do.call(rbind, lapply(seq_along(sel_ids), function(i){
        dplyr::tibble(
          value = as.numeric(Xtr[sel_ids[i], idx]),
          condition = factor(cond[idx], levels = by_cond$condition),
          condition_lab = factor(lab_map[cond[idx]], levels = lab_map[by_cond$condition]),
          feature_id = sel_ids[i],
          antibody = sel_labs[i]
        )
      }))
      win_df <- do.call(rbind, lapply(feature_windows, function(w)
        data.frame(feature_id = w$id, low = w$low, high = w$high, antibody = w$label)
      ))
      plt <- ggplot2::ggplot(long, ggplot2::aes(x = value, color = condition_lab, fill = condition_lab)) +
        ggplot2::geom_density(alpha = 0.25, na.rm = TRUE) +
        ggplot2::geom_vline(data = win_df, ggplot2::aes(xintercept = low), linetype = "dashed", inherit.aes = FALSE) +
        ggplot2::geom_vline(data = win_df, ggplot2::aes(xintercept = high), linetype = "dashed", inherit.aes = FALSE) +
        ggplot2::facet_wrap(~ antibody + feature_id, scales = "free_y") +
        ggplot2::labs(
          title = paste0("Keep within \u00B1", k, " MAD (per feature)"),
          subtitle = subtitle_txt,
          x = paste0("Signal (", transform,
                     if (identical(transform, "arcsinh")) paste0(", cofactor=", cofactor) else "", ")"),
          y = "Density", color = "Condition", fill = "Condition"
        ) + theme_bw(base_size = 12)
    }
  }
  
  list(
    se_filtered = se[, keep_idx, drop = FALSE],
    keep_index  = stats::setNames(keep_idx, colnames(se)),
    stats_overall = by_cond,
    stats_windows = per_feature,
    plot = plt,
    params = list(channel = channel, combine = combine, k = k,
                  transform = transform, cofactor = cofactor, assay = assay_name)
  )
}

# Classify cell cycle from a DNA-content channel (e.g., FxCycle Violet)
# Strategy: per condition, find two modes (≈2N and 4N), estimate robust spread
# via MAD around nearest-mode assignments, build windows mu ± k*MAD, and
# call phases: subG1, G0/G1, S (between windows), G2/M, >4N.
# Robust DNA-based cell-cycle classifier with debris handling
# Classify cell cycle from DNA channel with shared/per-condition windows
classify_cell_cycle_by_dna <- function(
    se,
    channel = "Cell.Cycle",
    assay_name = NULL,
    condition_col = "condition",
    label_field = "antibody",
    transform = c("arcsinh","none","log10"),
    cofactor = 150,
    k = 2,
    sample_per_condition = 40000,
    seed = 1,
    return_plot = TRUE,
    # robustness
    exclude_below = c("auto","quantile","fixed","none"),
    exclude_q = 0.01,
    exclude_fixed = NA_real_,
    min_sep_iqr = 0.5,
    # window fitting
    window_scope = c("per_condition","global","reference"),
    reference_level = NULL,
    min_s_frac = 0.06,
    # peak alignment
    align_peaks = TRUE,
    align_scope = c("global","reference"),
    # widening G1 left tail
    g1_left_pad_mad = 0
){
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(condition_col %in% colnames(SummarizedExperiment::colData(se)))
  transform     <- match.arg(transform)
  exclude_below <- match.arg(exclude_below)
  window_scope  <- match.arg(window_scope)
  align_scope   <- match.arg(align_scope)
  
  if (is.null(assay_name)) assay_name <- SummarizedExperiment::assayNames(se)[1]
  X <- SummarizedExperiment::assay(se, assay_name)
  ids <- rownames(se)
  cond <- as.character(SummarizedExperiment::colData(se)[[condition_col]])
  
  # labels -> row
  rdat <- SummarizedExperiment::rowData(se)
  labs <- if (label_field %in% colnames(rdat)) as.character(rdat[, label_field]) else ids
  labs[is.na(labs) | labs==""] <- ids[is.na(labs) | labs==""]
  hit <- which(labs == channel); if (!length(hit)) hit <- which(tolower(labs) == tolower(channel))
  if (!length(hit)) hit <- which(ids == channel)
  if (!length(hit)) stop("Channel/antibody not found: ", channel)
  f <- hit[1]; feature_id <- ids[f]; feature_lab <- labs[f]
  
  # transform
  Xtr <- switch(transform,
                "arcsinh" = asinh(X / cofactor),
                "log10"   = { tmp <- X; tmp[tmp <= 0] <- NA_real_; log10(tmp) },
                "none"    = X)
  v_raw <- as.numeric(Xtr[f, ])
  
  # -------- helpers --------
  otsu_threshold <- function(x, nb=256){
    x <- x[is.finite(x)]; if (!length(x)) return(NA_real_)
    h <- graphics::hist(x, breaks=nb, plot=FALSE)
    cs <- cumsum(h$counts); mu <- cumsum(h$counts*h$mids)
    T  <- sum(h$counts); muT <- mu[length(mu)]
    den <- cs*(T-cs); den[den==0] <- NA_real_
    h$mids[which.max((muT*cs - mu)^2 / den)]
  }
  find_two_modes_robust <- function(x, floor_val, min_sep_iqr){
    xf <- x[is.finite(x)]; if (!length(xf)) return(c(NA_real_, NA_real_))
    if (is.finite(floor_val)) xf <- xf[xf >= floor_val]
    if (length(xf) < 50) return(sort(stats::kmeans(xf, centers=2)$centers))
    d <- stats::density(xf, n=2048); y <- d$y; locs <- which(diff(sign(diff(y)))==-2)+1
    centers <- if (length(locs) >= 2) d$x[locs[order(y[locs], decreasing=TRUE)][1:2]]
    else sort(stats::kmeans(xf, centers=2)$centers)
    centers <- sort(centers)
    iqr <- diff(stats::quantile(xf, c(0.25,0.75), na.rm=TRUE)); if (!is.finite(iqr) || iqr==0) iqr <- stats::sd(xf, na.rm=TRUE)
    if (is.finite(iqr) && (centers[2]-centers[1]) < min_sep_iqr*iqr) {
      q <- stats::quantile(xf, c(0.3,0.8), na.rm=TRUE)
      centers <- sort(stats::kmeans(xf, centers=as.numeric(q))$centers)
    }
    centers
  }
  build_peaks <- function(x){
    xf <- x[is.finite(x)]; if (!length(xf)) return(c(NA_real_, NA_real_, NA_real_))
    floor_val <- switch(exclude_below,
                        "quantile" = stats::quantile(xf, exclude_q, na.rm=TRUE),
                        "fixed"    = exclude_fixed,
                        "auto"     = { med <- stats::median(xf, na.rm=TRUE)
                        xl  <- xf[is.finite(med) & xf <= med]
                        thr <- otsu_threshold(xl)
                        if (!length(thr) || !is.finite(thr)) stats::quantile(xf, exclude_q, na.rm=TRUE) else thr },
                        "none"     = -Inf)
    centers <- find_two_modes_robust(x, floor_val, min_sep_iqr)
    c(centers[1], centers[2], floor_val)
  }
  build_windows_from_peaks <- function(x, mu1, mu2){
    grp <- ifelse(abs(x - mu1) <= abs(x - mu2), 1L, 2L)
    mad1 <- stats::mad(x[grp==1], center=mu1, constant=1, na.rm=TRUE)
    mad2 <- stats::mad(x[grp==2], center=mu2, constant=1, na.rm=TRUE)
    if (!is.finite(mad1) || mad1==0) { iqr <- diff(stats::quantile(x[grp==1], c(0.25,0.75), na.rm=TRUE)); mad1 <- ifelse(is.finite(iqr) && iqr>0, iqr/1.349, stats::sd(x[grp==1], na.rm=TRUE)) }
    if (!is.finite(mad2) || mad2==0) { iqr <- diff(stats::quantile(x[grp==2], c(0.25,0.75), na.rm=TRUE)); mad2 <- ifelse(is.finite(iqr) && iqr>0, iqr/1.349, stats::sd(x[grp==2], na.rm=TRUE)) }
    low1 <- mu1 - k*mad1; high1 <- mu1 + k*mad1
    low2 <- mu2 - k*mad2; high2 <- mu2 + k*mad2
    low1 <- low1 - g1_left_pad_mad * mad1
    # enforce minimum S-gap
    sep <- mu2 - mu1
    if (is.finite(sep) && sep > 0) {
      min_gap <- min_s_frac * sep
      if ((low2 - high1) < min_gap) {
        mid <- (mu1 + mu2)/2
        high1 <- min(mid - min_gap/2, high1)
        low2  <- max(mid + min_gap/2, low2)
      }
      if (high1 >= low2) { mid <- (mu1 + mu2)/2; high1 <- min(high1, mid); low2 <- max(low2, mid) }
    }
    list(mu1=mu1, mu2=mu2, low1=low1, high1=high1, low2=low2, high2=high2)
  }
  
  # -------- optional peak alignment --------
  cond_levels <- unique(cond)
  peaks_by_cond <- lapply(cond_levels, function(cc) build_peaks(v_raw[cond==cc]))
  names(peaks_by_cond) <- cond_levels
  
  if (align_peaks) {
    if (align_scope == "reference") {
      stopifnot(!is.null(reference_level), reference_level %in% cond_levels)
      target <- peaks_by_cond[[reference_level]][1:2]
    } else {
      mu <- build_peaks(v_raw)
      target <- mu[1:2]
    }
  }
  
  v <- v_raw
  if (align_peaks) {
    for (cc in cond_levels) {
      mu1c <- peaks_by_cond[[cc]][1]; mu2c <- peaks_by_cond[[cc]][2]
      if (any(!is.finite(c(mu1c, mu2c, target)))) next
      a <- (target[2] - target[1]) / (mu2c - mu1c)
      b <- target[1] - a * mu1c
      v[cond==cc] <- a * v_raw[cond==cc] + b
    }
  }
  
  # -------- fit windows --------
  fit_groups <- switch(window_scope,
                       "per_condition" = setNames(lapply(cond_levels, function(cc) which(cond==cc)), cond_levels),
                       "global"        = list(ALL = seq_along(cond)),
                       "reference"     = { stopifnot(!is.null(reference_level), reference_level %in% cond_levels)
                         setNames(list(which(cond==reference_level)), reference_level) })
  
  win_fit <- lapply(fit_groups, function(ix){
    pk <- build_peaks(v[ix])
    w  <- build_windows_from_peaks(v[ix], mu1 = pk[1], mu2 = pk[2])
    data.frame(mu1=w$mu1, mu2=w$mu2, low1=w$low1, high1=w$high1, low2=w$low2, high2=w$high2)
  })
  
  win_map <- list()
  if (window_scope == "per_condition") {
    for (cc in names(fit_groups)) win_map[[cc]] <- win_fit[[cc]]
  } else {
    ref_name <- names(win_fit)[1]
    for (cc in cond_levels) win_map[[cc]] <- win_fit[[ref_name]]
  }
  
  # -------- classify to 3 phases --------
  calls <- character(length(v))
  win_list <- list()
  for (cc in cond_levels) {
    w <- win_map[[cc]]
    ix <- which(cond==cc); x <- v[ix]
    call_cc <- ifelse(!is.finite(x), NA_character_,
                      ifelse(x <= w$high1, "G0/G1",
                             ifelse(x <  w$low2, "S", "G2/M")))
    calls[ix] <- call_cc
    w$condition <- cc; win_list[[cc]] <- w
  }
  
  SummarizedExperiment::colData(se)$cell_cycle <- factor(
    calls, levels = c("G0/G1","S","G2/M")
  )
  SummarizedExperiment::colData(se)$cell_cycle_channel <- feature_id
  
  stats_tbl <- dplyr::as_tibble(SummarizedExperiment::colData(se)[, condition_col, drop=FALSE]) |>
    dplyr::mutate(phase = SummarizedExperiment::colData(se)$cell_cycle) |>
    dplyr::count(.data[[condition_col]], phase, name="n") |>
    dplyr::group_by(.data[[condition_col]]) |>
    dplyr::mutate(pct = 100*n/sum(n)) |>
    dplyr::ungroup()
  
  win_df <- do.call(rbind, win_list)
  
  # -------- plot --------
  plt <- NULL
  if (return_plot) {
    set.seed(seed)
    idx_plot <- unlist(lapply(split(seq_along(cond), cond), function(ix){
      if (length(ix) > sample_per_condition) sample(ix, sample_per_condition) else ix
    }), use.names=FALSE)
    
    df_plot <- dplyr::tibble(value = v[idx_plot],
                             condition = factor(cond[idx_plot], levels = cond_levels))
    
    rects <- win_df |>
      dplyr::select(condition, low1, high1, low2, high2) |>
      dplyr::mutate(S_low = pmax(high1, low1), S_high = pmin(low2, high2))
    
    byc <- stats_tbl |>
      dplyr::group_by(.data[[condition_col]]) |>
      dplyr::summarise(G1 = sum(pct[phase=="G0/G1"], na.rm=TRUE),
                       S  = sum(pct[phase=="S"],    na.rm=TRUE),
                       G2 = sum(pct[phase=="G2/M"], na.rm=TRUE),
                       .groups="drop")
    subt <- paste(byc[[condition_col]],
                  sprintf("G1 %.1f%%  S %.1f%%  G2/M %.1f%%", byc$G1, byc$S, byc$G2),
                  collapse="   ")
    
    plt <- ggplot2::ggplot(df_plot, ggplot2::aes(x=value)) +
      ggplot2::geom_density(fill="grey80", color="grey40", alpha=0.4, na.rm=TRUE) +
      ggplot2::geom_rect(data=rects,
                         ggplot2::aes(xmin=low1, xmax=high1, ymin=-Inf, ymax=Inf),
                         inherit.aes=FALSE, alpha=0.10, fill="#56B4E9") +
      ggplot2::geom_rect(data=rects,
                         ggplot2::aes(xmin=pmax(high1,low1), xmax=pmin(low2,high2),
                                      ymin=-Inf, ymax=Inf),
                         inherit.aes=FALSE, alpha=0.10, fill="#F0E442") +
      ggplot2::geom_rect(data=rects,
                         ggplot2::aes(xmin=low2, xmax=high2, ymin=-Inf, ymax=Inf),
                         inherit.aes=FALSE, alpha=0.10, fill="#E69F00") +
      ggplot2::geom_vline(data=win_df, ggplot2::aes(xintercept=low1),  linetype="dashed") +
      ggplot2::geom_vline(data=win_df, ggplot2::aes(xintercept=high1), linetype="dashed") +
      ggplot2::geom_vline(data=win_df, ggplot2::aes(xintercept=low2),  linetype="dashed") +
      ggplot2::geom_vline(data=win_df, ggplot2::aes(xintercept=high2), linetype="dashed") +
      ggplot2::facet_wrap(~ condition, scales="free_y") +
      ggplot2::labs(
        title = paste0("Cell cycle windows — ", feature_lab, " (", feature_id, ")"),
        subtitle = paste0(subt, if (align_peaks) " | peaks aligned" else ""),
        x = paste0("DNA signal (", transform,
                   if (identical(transform,"arcsinh")) paste0(", cofactor=", cofactor) else "", ")"),
        y = "Density"
      ) + ggplot2::theme_bw(base_size=12)
  }
  
  list(
    se_with_phase = se,
    windows = win_df,
    stats = stats_tbl,
    plot = plt,
    params = list(channel=channel, assay=assay_name, k=k,
                  transform=transform, cofactor=cofactor,
                  exclude_below=exclude_below, exclude_q=exclude_q,
                  exclude_fixed=exclude_fixed, min_sep_iqr=min_sep_iqr,
                  window_scope=window_scope, reference_level=reference_level,
                  min_s_frac=min_s_frac, align_peaks=align_peaks,
                  align_scope=align_scope,
                  g1_left_pad_mad=g1_left_pad_mad)
  )
}


# Singlet gate using FSC-A vs FSC-H (preferred) or FSC-A vs FSC-W (fallback)
gate_singlets <- function(
    se,
    area_candidates   = c("FSC-A","FSC_A","FSC.A"),
    height_candidates = c("FSC-H","FSC_H","FSC.H"),
    width_candidates  = c("FSC-W","FSC_W","FSC.W"),
    method = c("auto","AH","AW"),        # AH = Area~Height band; AW = Area~Width one-sided
    assay_name = NULL,
    condition_col = "condition",
    label_field = "antibody",
    transform = c("arcsinh","none","log10"),
    cofactor = 150,
    k = 3,                               # band width in MAD units
    per_condition = TRUE,
    max_points_per_condition = 50000,    # for plotting only
    seed = 1
){
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(condition_col %in% colnames(SummarizedExperiment::colData(se)))
  transform <- match.arg(transform)
  method    <- match.arg(method)
  
  if (is.null(assay_name)) assay_name <- SummarizedExperiment::assayNames(se)[1]
  X    <- SummarizedExperiment::assay(se, assay_name)
  ids  <- rownames(se)
  cond <- as.character(SummarizedExperiment::colData(se)[[condition_col]])
  
  # map labels -> rownames
  rdat <- SummarizedExperiment::rowData(se)
  labs <- if (label_field %in% colnames(rdat)) as.character(rdat[, label_field]) else ids
  labs[is.na(labs) | labs==""] <- ids[is.na(labs) | labs==""]
  
  resolve_row <- function(cands){
    hit <- which(ids %in% cands | labs %in% cands)
    if (!length(hit)) hit <- grep(paste(cands, collapse="|"), ids, ignore.case=TRUE)
    if (!length(hit)) hit <- grep(paste(cands, collapse="|"), labs, ignore.case=TRUE)
    if (!length(hit)) return(NA_character_)
    ids[hit[1]]
  }
  id_A <- resolve_row(area_candidates)
  id_H <- resolve_row(height_candidates)
  id_W <- resolve_row(width_candidates)
  
  if (method == "auto") {
    method <- if (!is.na(id_H)) "AH" else if (!is.na(id_W)) "AW" else stop("No FSC-H or FSC-W channel found.")
  } else {
    if (method == "AH" && is.na(id_H)) stop("Requested AH but no height channel found.")
    if (method == "AW" && is.na(id_W)) stop("Requested AW but no width channel found.")
  }
  
  # transform
  Xtr <- switch(transform,
                arcsinh = asinh(X / cofactor),
                log10   = { tmp <- X; tmp[tmp <= 0] <- NA_real_; log10(tmp) },
                none    = X
  )
  A <- as.numeric(Xtr[id_A, ])
  H <- if (!is.na(id_H)) as.numeric(Xtr[id_H, ]) else NULL
  W <- if (!is.na(id_W)) as.numeric(Xtr[id_W, ]) else NULL
  
  # robust fit helpers
  robust_band_keep <- function(x, y, k, side = "both"){  # side = "both" for AH, "below" for AW
    ok <- stats::complete.cases(x, y)
    if (!any(ok)) return(rep(FALSE, length(x)))
    df <- data.frame(x = x[ok], y = y[ok])
    # robust line (M-estimator); fall back to lm if MASS not available
    fit <- tryCatch(MASS::rlm(y ~ x, data = df, maxit = 50), error = function(e) stats::lm(y ~ x, data = df))
    yhat <- stats::predict(fit, newdata = data.frame(x = x))
    resid <- y - yhat
    s <- stats::mad(resid[ok], constant = 1, na.rm = TRUE)
    if (!is.finite(s) || s == 0) s <- stats::sd(resid[ok], na.rm = TRUE)
    if (!is.finite(s) || s == 0) s <- 0
    if (side == "both") {
      keep <- abs(resid) <= k * s
    } else { # "below": keep narrow pulses (singlets) below the upper band
      keep <- resid <= k * s
    }
    keep[!ok] <- FALSE
    list(keep = keep, yhat = yhat, s = s, fit = fit)
  }
  
  keep <- rep(FALSE, length(A))
  stats_list <- list()
  plots <- list()
  
  set.seed(seed)
  groups <- if (per_condition) split(seq_along(cond), cond) else list(ALL = seq_along(cond))
  
  for (g in names(groups)) {
    ix <- groups[[g]]
    if (method == "AH") {
      res <- robust_band_keep(A[ix], H[ix], k, side = "both")
      keep[ix] <- res$keep
      # plotting sample
      samp <- if (length(ix) > max_points_per_condition) sample(ix, max_points_per_condition) else ix
      dfp <- data.frame(A = A[samp], H = H[samp], keep = res$keep[match(samp, ix)])
      # band lines
      ah <- data.frame(A = A[samp], yhat = res$yhat[match(samp, ix)])
      p <- ggplot2::ggplot(dfp, ggplot2::aes(A, H)) +
        ggplot2::geom_point(ggplot2::aes(color = keep), alpha = 0.3, size = 0.4) +
        ggplot2::geom_line(data = ah[order(ah$A), ], ggplot2::aes(A, yhat), color = "black") +
        ggplot2::geom_line(data = ah[order(ah$A), ], ggplot2::aes(A, yhat + k*res$s), linetype = "dashed") +
        ggplot2::geom_line(data = ah[order(ah$A), ], ggplot2::aes(A, yhat - k*res$s), linetype = "dashed") +
        ggplot2::labs(
          title = sprintf("Singlet gate (FSC-A vs FSC-H) — %s", if (per_condition) g else "all"),
          x = paste0(id_A, " (", transform, if (identical(transform,"arcsinh")) paste0(", c=",cofactor) else")"),
          y = paste0(id_H, " (", transform, if (identical(transform,"arcsinh")) paste0(", c=",cofactor) else")"),
          color = "Kept"
        ) + ggplot2::theme_bw(base_size = 12)
      plots[[g]] <- p
    } else { # AW
      res <- robust_band_keep(A[ix], W[ix], k, side = "below")
      keep[ix] <- res$keep
      samp <- if (length(ix) > max_points_per_condition) sample(ix, max_points_per_condition) else ix
      dfp <- data.frame(A = A[samp], W = W[samp], keep = res$keep[match(samp, ix)])
      aw <- data.frame(A = A[samp], yhat = res$yhat[match(samp, ix)])
      p <- ggplot2::ggplot(dfp, ggplot2::aes(A, W)) +
        ggplot2::geom_point(ggplot2::aes(color = keep), alpha = 0.3, size = 0.4) +
        ggplot2::geom_line(data = aw[order(aw$A), ], ggplot2::aes(A, yhat), color = "black") +
        ggplot2::geom_line(data = aw[order(aw$A), ], ggplot2::aes(A, yhat + k*res$s), linetype = "dashed") +
        ggplot2::labs(
          title = sprintf("Singlet gate (FSC-A vs FSC-W) — %s", if (per_condition) g else "all"),
          x = paste0(id_A, " (", transform, if (identical(transform,"arcsinh")) paste0(", c=",cofactor) else")"),
          y = paste0(id_W, " (", transform, if (identical(transform,"arcsinh")) paste0(", c=",cofactor) else")"),
          color = "Kept"
        ) + ggplot2::theme_bw(base_size = 12)
      plots[[g]] <- p
    }
    
    stats_list[[g]] <- data.frame(
      group = g, n_total = length(ix),
      n_keep = sum(keep[ix], na.rm = TRUE),
      pct_keep = 100 * sum(keep[ix], na.rm = TRUE) / length(ix),
      method = method, k = k,
      area = id_A, height = id_H, width = id_W,
      stringsAsFactors = FALSE
    )
  }
  
  stats_df <- do.call(rbind, stats_list)
  
  list(
    se_filtered = se[, keep, drop = FALSE],
    keep_index  = stats::setNames(keep, colnames(se)),
    stats = stats_df,
    plots = plots,
    params = list(method = method, k = k, per_condition = per_condition,
                  transform = transform, cofactor = cofactor,
                  area = id_A, height = id_H, width = id_W, assay = assay_name)
  )
}

# Debris gate with hard floors + optional ellipse; supports global or per-condition
gate_debris_fsc_ssc <- function(
    se,
    fsc = c("FSC-A","FSC_A","FSC.A"),
    ssc = c("SSC-A","SSC_A","SSC.A"),
    assay_name = NULL,
    condition_col = "condition",
    label_field = "antibody",
    transform = c("arcsinh","none","log10"),
    cofactor = 150,
    keep_fraction = 0.95,          # ellipse mass to keep (if used)
    per_condition = TRUE,
    # floors (hard lower bounds; anything below is always removed)
    floor_mode = c("none","quantile","fixed"),
    floor_q = 0.02,                # if floor_mode="quantile"
    fsc_floor_fixed = NA_real_,    # if floor_mode="fixed"
    ssc_floor_fixed = NA_real_,    # if floor_mode="fixed"
    # how to combine floors and ellipse
    keep_rule = c("floors_and_ellipse","floors_or_ellipse","floors_only","ellipse_only"),
    seed = 1
){
  stopifnot(inherits(se, "SummarizedExperiment"))
  transform  <- match.arg(transform)
  floor_mode <- match.arg(floor_mode)
  keep_rule  <- match.arg(keep_rule)
  if (!per_condition && !is.null(condition_col) && condition_col %in% colnames(SummarizedExperiment::colData(se))) {
    cond <- rep("ALL", ncol(se))
  } else {
    cond <- if (!is.null(condition_col) && condition_col %in% colnames(SummarizedExperiment::colData(se)))
      as.character(SummarizedExperiment::colData(se)[[condition_col]]) else rep("ALL", ncol(se))
  }
  if (is.null(assay_name)) assay_name <- SummarizedExperiment::assayNames(se)[1]
  
  ids  <- rownames(se)
  rdat <- SummarizedExperiment::rowData(se)
  labs <- if (label_field %in% colnames(rdat)) as.character(rdat[, label_field]) else ids
  labs[is.na(labs) | labs==""] <- ids[is.na(labs) | labs==""]
  resolve_row <- function(cands){
    hit <- which(ids %in% cands | labs %in% cands)
    if (!length(hit)) {
      pat <- paste0("^(", paste(gsub("\\.", "\\\\.", cands), collapse="|"), ")$")
      hit <- grep(pat, ids, ignore.case = TRUE); if (!length(hit)) hit <- grep(pat, labs, ignore.case = TRUE)
    }
    if (!length(hit)) stop("Could not find channel matching: ", paste(cands, collapse=", "))
    ids[hit[1]]
  }
  fsc_id <- resolve_row(fsc)
  ssc_id <- resolve_row(ssc)
  
  X  <- SummarizedExperiment::assay(se, assay_name)
  Xt <- switch(transform,
               arcsinh = asinh(X / cofactor),
               log10   = { tmp <- X; tmp[tmp <= 0] <- NA_real_; log10(tmp) },
               none    = X
  )
  FSC <- as.numeric(Xt[fsc_id, ]); SSC <- as.numeric(Xt[ssc_id, ])
  
  groups <- if (per_condition) split(seq_along(cond), cond) else list(ALL = seq_along(cond))
  ellipse_points <- function(center, covmat, level = keep_fraction, n = 200){
    ev <- eigen(covmat); r <- sqrt(stats::qchisq(level, df = 2))
    ang <- seq(0, 2*pi, length.out = n)
    sweep((r * cbind(cos(ang), sin(ang))) %*% t(ev$vectors %*% diag(sqrt(ev$values))), 2, center, `+`)
  }
  
  keep <- rep(FALSE, length(FSC))
  stats_list <- list(); plot_list <- list()
  
  set.seed(seed)
  for (g in names(groups)) {
    ix <- groups[[g]]
    
    # ----- floors (hard)
    if (floor_mode == "none") {
      f_floor <- -Inf; s_floor <- -Inf
    } else if (floor_mode == "quantile") {
      f_floor <- stats::quantile(FSC[ix], floor_q, na.rm = TRUE)
      s_floor <- stats::quantile(SSC[ix], floor_q, na.rm = TRUE)
    } else {
      f_floor <- fsc_floor_fixed; s_floor <- ssc_floor_fixed
    }
    pass_floor <- (FSC >= f_floor | !is.finite(f_floor)) & (SSC >= s_floor | !is.finite(s_floor))
    
    # ----- ellipse (fit only on points above floors)
    ok <- stats::complete.cases(FSC[ix], SSC[ix]) & pass_floor[ix]
    inside <- rep(FALSE, length(ix))
    ell_df <- NULL
    if (any(ok)) {
      xy <- cbind(FSC[ix][ok], SSC[ix][ok])
      center <- c(stats::median(xy[,1]), stats::median(xy[,2]))
      covmat <- tryCatch(MASS::cov.trob(xy)$cov, error = function(e) stats::cov(xy))
      d2 <- stats::mahalanobis(cbind(FSC[ix], SSC[ix]),
                               center = center, cov = covmat)
      inside <- is.finite(d2) & (d2 <= stats::qchisq(keep_fraction, df = 2))
      # for plotting
      ell_df <- as.data.frame(ellipse_points(center, covmat, keep_fraction))
    }
    
    # ----- combine rules
    kept_local <- switch(keep_rule,
                         floors_and_ellipse = pass_floor[ix] & inside,
                         floors_or_ellipse  = pass_floor[ix] | inside,
                         floors_only        = pass_floor[ix],
                         ellipse_only       = inside
    )
    keep[ix] <- kept_local
    
    stats_list[[g]] <- data.frame(
      group = g,
      n_total = length(ix),
      n_pass_floor = sum(pass_floor[ix], na.rm = TRUE),
      n_keep = sum(kept_local, na.rm = TRUE),
      pct_keep = 100 * sum(kept_local, na.rm = TRUE) / length(ix),
      fsc = fsc_id, ssc = ssc_id, keep_fraction = keep_fraction,
      floor_mode = floor_mode, floor_q = ifelse(floor_mode=="quantile", floor_q, NA),
      fsc_floor = f_floor, ssc_floor = s_floor,
      keep_rule = keep_rule,
      stringsAsFactors = FALSE
    )
    
    # ----- plot
    dfp <- data.frame(FSC = FSC[ix], SSC = SSC[ix], keep = kept_local)
    dfp <- dfp[stats::complete.cases(dfp), , drop = FALSE]
    p <- ggplot2::ggplot(dfp, ggplot2::aes(FSC, SSC)) +
      ggplot2::geom_point(ggplot2::aes(color = keep), alpha = 0.3, size = 0.4) +
      (if (!is.null(ell_df))
        ggplot2::geom_path(data = ell_df, ggplot2::aes(V1, V2), color = "black") else NULL) +
      (if (is.finite(f_floor)) ggplot2::geom_vline(xintercept = f_floor, linetype = "longdash", color = "red", alpha = 0.7) else NULL) +
      (if (is.finite(s_floor)) ggplot2::geom_hline(yintercept = s_floor, linetype = "longdash", color = "red", alpha = 0.7) else NULL) +
      ggplot2::labs(
        title = if (per_condition) sprintf("Debris gate: %s", g) else "Debris gate (all conditions)",
        subtitle = sprintf("Kept %.1f%% | floors: FSC%s, SSC%s | rule: %s",
                           stats_list[[g]]$pct_keep,
                           if (is.finite(f_floor)) sprintf(" ≥ %.3f", f_floor) else " —",
                           if (is.finite(s_floor)) sprintf(" ≥ %.3f", s_floor) else " —",
                           keep_rule),
        x = paste0(fsc_id, " (", transform, if (identical(transform,"arcsinh")) paste0(", c=",cofactor) else "", ")"),
        y = paste0(ssc_id, " (", transform, if (identical(transform,"arcsinh")) paste0(", c=",cofactor) else "", ")"),
        color = "Kept"
      ) + ggplot2::theme_bw(base_size = 12)
    plot_list[[g]] <- p
  }
  
  list(
    se_filtered = se[, keep, drop = FALSE],
    keep_index  = stats::setNames(keep, colnames(se)),
    stats = do.call(rbind, stats_list),
    plots = plot_list,
    params = list(fsc = fsc_id, ssc = ssc_id,
                  transform = transform, cofactor = cofactor,
                  keep_fraction = keep_fraction, per_condition = per_condition,
                  floor_mode = floor_mode, floor_q = floor_q,
                  fsc_floor_fixed = fsc_floor_fixed, ssc_floor_fixed = ssc_floor_fixed,
                  keep_rule = keep_rule)
  )
}

# Normalize selected channels by a reference channel (e.g., H3) and plot before/after
normalize_channels_by_reference <- function(
    se,
    targets        = c("H3K14ac","H3K23ac"),   # antibody labels or rownames (features)
    reference      = "H3",                     # antibody label or rowname
    assay_in       = NULL,                     # source assay (default: first)
    assay_out      = "mods_by_H3",             # name of new assay to create
    label_field    = "antibody",               # column in rowData with antibody names
    match          = c("exact","ignore_case","substring"),
    ref_floor      = 0,                        # denom <= ref_floor -> NA
    multiplier     = 1,                        # optional rescale of ratios
    keep_others    = TRUE,                     # keep non-target rows as-is in assay_out
    overwrite      = TRUE,                     # allow overwriting existing assay
    # ---- plotting ----
    condition_col  = "condition",
    make_plots     = TRUE,
    plot_transform_before = c("none","log10","arcsinh"),
    plot_transform_after  = c("none","log10","arcsinh"),
    cofactor       = 150,                      # only used if plot transform = arcsinh
    density_facet_ncol = 3,
    max_cells_per_condition = 20000,
    seed = 1
){
  stopifnot(inherits(se, "SummarizedExperiment"))
  match <- match.arg(match)
  plot_transform_before <- match.arg(plot_transform_before)
  plot_transform_after  <- match.arg(plot_transform_after)
  
  if (is.null(assay_in)) assay_in <- SummarizedExperiment::assayNames(se)[1]
  X <- SummarizedExperiment::assay(se, assay_in)   # features x cells
  
  ids  <- rownames(se)
  rdat <- SummarizedExperiment::rowData(se)
  labs <- if (label_field %in% colnames(rdat)) as.character(rdat[, label_field]) else ids
  labs[is.na(labs) | labs==""] <- ids[is.na(labs) | labs==""]
  
  # resolver
  resolve <- function(q) {
    hit <- switch(match,
                  exact       = which(labs == q | ids == q),
                  ignore_case = which(tolower(labs) == tolower(q) | tolower(ids) == tolower(q)),
                  substring   = unique(c(grep(q, labs, ignore.case=TRUE), grep(q, ids, ignore.case=TRUE)))
    )
    hit
  }
  
  # reference
  ref_idx <- resolve(reference)
  if (length(ref_idx) == 0) stop("Reference not found: ", reference)
  if (length(ref_idx) > 1) { warning("Multiple matches for reference; using first: ", ids[ref_idx[1]]); ref_idx <- ref_idx[1] }
  
  # targets
  tgt_idx <- unique(unlist(lapply(targets, resolve)))
  if (length(tgt_idx) == 0) stop("No targets matched: ", paste(targets, collapse=", "))
  
  # denominator with floor → NA
  den <- as.numeric(X[ref_idx, ])
  den_ok <- ifelse(is.finite(den) & den > ref_floor, den, NA_real_)
  
  # output assay
  Y <- if (keep_others) X else matrix(NA_real_, nrow=nrow(X), ncol=ncol(X), dimnames=dimnames(X))
  for (i in tgt_idx) Y[i, ] <- multiplier * (as.numeric(X[i, ]) / den_ok)
  
  if (!overwrite && assay_out %in% SummarizedExperiment::assayNames(se)) {
    stop("Assay '", assay_out, "' already exists. Set overwrite=TRUE to replace.")
  }
  SummarizedExperiment::assays(se)[[assay_out]] <- Y
  
  # ---- plots (densities per condition; before vs after) ----
  plots <- NULL
  if (make_plots) {
    # condition (graceful if missing)
    if (condition_col %in% colnames(SummarizedExperiment::colData(se))) {
      cond <- as.character(SummarizedExperiment::colData(se)[[condition_col]])
    } else {
      cond <- rep("ALL", ncol(se))
    }
    
    # helper: transform choice
    tfun <- function(v, how){
      switch(how,
             "none"   = v,
             "log10"  = { z <- v; z[!is.finite(z) | z <= 0] <- NA_real_; log10(z) },
             "arcsinh"= asinh(v / cofactor)
      )
    }
    
    # sample per condition to keep plotting light
    set.seed(seed)
    idx <- unlist(lapply(split(seq_along(cond), cond), function(ix){
      if (length(ix) > max_cells_per_condition) sample(ix, max_cells_per_condition) else ix
    }), use.names = FALSE)
    
    # BEFORE
    df_b <- as.data.frame(t(X[tgt_idx, idx, drop=FALSE]), check.names = FALSE)
    df_b[[condition_col]] <- cond[idx]
    long_b <- df_b |>
      tidyr::pivot_longer(cols = tidyselect::all_of(ids[tgt_idx]),
                          names_to = "feature_id", values_to = "value") |>
      dplyr::mutate(
        feature_label = labs[match(feature_id, ids)],
        value = tfun(value, plot_transform_before),
        panel = "Before (raw)"
      )
    
    p_before <- ggplot2::ggplot(long_b,
                                ggplot2::aes(x = value, fill = .data[[condition_col]])) +
      ggplot2::geom_density(alpha = 0.35, na.rm = TRUE) +
      ggplot2::facet_wrap(~ feature_label, scales = "free_y", ncol = density_facet_ncol) +
      ggplot2::labs(
        title = sprintf("Targets before normalization — assay '%s'", assay_in),
        x = sprintf("Value (%s%s)",
                    plot_transform_before,
                    if (plot_transform_before=="arcsinh") sprintf(", c=%s", cofactor) else ""),
        y = "Density", fill = "Condition"
      ) + ggplot2::theme_bw(base_size = 12)
    
    # AFTER
    df_a <- as.data.frame(t(Y[tgt_idx, idx, drop=FALSE]), check.names = FALSE)
    df_a[[condition_col]] <- cond[idx]
    long_a <- df_a |>
      tidyr::pivot_longer(cols = tidyselect::all_of(ids[tgt_idx]),
                          names_to = "feature_id", values_to = "value") |>
      dplyr::mutate(
        feature_label = paste0(labs[match(feature_id, ids)], " / ", labs[ref_idx]),
        value = tfun(value, plot_transform_after),
        panel = "After (ratio)"
      )
    
    p_after <- ggplot2::ggplot(long_a,
                               ggplot2::aes(x = value, fill = .data[[condition_col]])) +
      ggplot2::geom_density(alpha = 0.35, na.rm = TRUE) +
      ggplot2::facet_wrap(~ feature_label, scales = "free_y", ncol = density_facet_ncol) +
      ggplot2::labs(
        title = sprintf("Targets after normalization by %s — assay '%s'",
                        labs[ref_idx], assay_out),
        x = sprintf("Value (%s%s)%s",
                    plot_transform_after,
                    if (plot_transform_after=="arcsinh") sprintf(", c=%s", cofactor) else "",
                    if (multiplier != 1) sprintf(" × %.3g", multiplier) else ""),
        y = "Density", fill = "Condition"
      ) + ggplot2::theme_bw(base_size = 12)
    
    plots <- list(before = p_before, after = p_after)
  }
  
  # summary
  summ <- data.frame(
    target_id     = ids[tgt_idx],
    target_lab    = labs[tgt_idx],
    reference_id  = ids[ref_idx],
    reference_lab = labs[ref_idx],
    n_cells       = ncol(se),
    n_ref_NA      = sum(!is.finite(den_ok)),
    ref_floor     = ref_floor,
    multiplier    = multiplier,
    assay_in      = assay_in,
    assay_out     = assay_out,
    stringsAsFactors = FALSE
  )
  
  list(se = se, assay_out = assay_out, summary = summ, plots = plots)
}

plot_ptm_covariation <- function(
    se,
    x_feature, y_feature,
    assay_name = NULL,
    condition_col = "condition",
    label_field = "antibody",
    transform = c("none","arcsinh","log10"),
    cofactor = 150,
    cor_method = c("pearson","spearman"),
    sample_per_condition = 50000,
    colors = NULL,
    axis_label_prefix = '',
    outlier = c("none","quantile","mad"),
    outlier_q = c(0.01, 0.99),
    outlier_k = 5,
    compute_stats_on = c("filtered","raw"),
    plot_filtered_only = TRUE,
    point_alpha = 0.25,
    point_size  = 0.4,
    title = '',
    seed = 1
){
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(condition_col %in% colnames(SummarizedExperiment::colData(se)))
  if (is.null(assay_name)) assay_name <- SummarizedExperiment::assayNames(se)[1]
  transform        <- match.arg(transform)
  cor_method       <- match.arg(cor_method)
  outlier          <- match.arg(outlier)
  compute_stats_on <- match.arg(compute_stats_on)
  
  X   <- SummarizedExperiment::assay(se, assay_name)
  ids <- rownames(se)
  cd  <- as.data.frame(SummarizedExperiment::colData(se))
  cond <- as.character(cd[[condition_col]])
  cond_levels <- unique(cond)
  
  rdat <- SummarizedExperiment::rowData(se)
  labs <- if (label_field %in% colnames(rdat)) as.character(rdat[, label_field]) else ids
  labs[is.na(labs) | labs==""] <- ids[is.na(labs) | labs==""]
  
  resolve <- function(q){
    hit <- which(labs == q | ids == q)
    if (!length(hit)) hit <- which(tolower(labs) == tolower(q) | tolower(ids) == tolower(q))
    if (!length(hit)) hit <- unique(c(grep(q, labs, ignore.case=TRUE), grep(q, ids, ignore.case=TRUE)))
    if (!length(hit)) stop("Feature not found: ", q)
    hit[1]
  }
  ix <- resolve(x_feature); iy <- resolve(y_feature)
  
  Tfun <- switch(transform,
                 "none"    = function(v) v,
                 "log10"   = function(v){ v[!is.finite(v) | v <= 0] <- NA_real_; log10(v) },
                 "arcsinh" = function(v) asinh(v / cofactor)
  )
  x_all <- Tfun(as.numeric(X[ix, ])); y_all <- Tfun(as.numeric(X[iy, ]))
  
  set.seed(seed)
  keep_idx <- unlist(lapply(split(seq_along(cond), cond), function(ii){
    if (length(ii) > sample_per_condition) sample(ii, sample_per_condition) else ii
  }), use.names = FALSE)
  
  df <- data.frame(
    x = x_all[keep_idx],
    y = y_all[keep_idx],
    condition = factor(cond[keep_idx], levels = cond_levels)
  )
  
  # Outlier filtering
  filter_one <- function(d){
    d <- d[is.finite(d$x) & is.finite(d$y), , drop = FALSE]
    if (!nrow(d) || outlier == "none") { d$keep <- TRUE; return(d) }
    if (outlier == "quantile") {
      qx <- stats::quantile(d$x, outlier_q, na.rm = TRUE)
      qy <- stats::quantile(d$y, outlier_q, na.rm = TRUE)
      d$keep <- (d$x >= qx[1] & d$x <= qx[2] & d$y >= qy[1] & d$y <= qy[2])
    } else {
      medx <- stats::median(d$x, na.rm=TRUE); madx <- stats::mad(d$x, constant=1, na.rm=TRUE); if (!is.finite(madx) || madx==0) madx <- stats::sd(d$x, na.rm=TRUE)
      medy <- stats::median(d$y, na.rm=TRUE); mady <- stats::mad(d$y, constant=1, na.rm=TRUE); if (!is.finite(mady) || mady==0) mady <- stats::sd(d$y, na.rm=TRUE)
      d$keep <- (abs(d$x - medx) <= outlier_k*madx) & (abs(d$y - medy) <= outlier_k*mady)
    }
    d
  }
  df <- do.call(rbind, lapply(split(df, df$condition), filter_one))
  
  df_use   <- df[df$keep, , drop = FALSE]
  df_plot  <- if (plot_filtered_only) df_use else df
  stats_df <- if (compute_stats_on == "filtered") df_use else df
  
  x_limits <- range(df_plot$x, na.rm = TRUE)
  y_limits <- range(df_plot$y, na.rm = TRUE)
  
  stat_tbl <- do.call(rbind, lapply(levels(df$condition), function(cc){
    d <- stats_df[stats_df$condition == cc & is.finite(stats_df$x) & is.finite(stats_df$y), , drop = FALSE]
    if (nrow(d) < 3) return(data.frame(condition = cc, n = nrow(d), r = NA_real_, p = NA_real_))
    r <- suppressWarnings(stats::cor(d$x, d$y, method = cor_method, use = "complete.obs"))
    p <- tryCatch(stats::cor.test(d$x, d$y, method = cor_method)$p.value, error = function(e) NA_real_)
    data.frame(condition = cc, n = nrow(d), r = r, p = p)
  }))
  
  # Keep original color set
  if (is.null(colors)) {
    base_cols <- c("#317ec2", "#c03830", "#2ca25f", "#756bb1")
    colors <- setNames(rep(base_cols, length.out = length(cond_levels)), cond_levels)
  } else {
    colors <- setNames(rep(colors, length.out = length(cond_levels)), cond_levels)
  }
  
  if (is.null(axis_label_prefix)) {
    axis_label_prefix <- if (grepl("norm", assay_name, ignore.case = TRUE)) "Normalized"
    else if (transform=="none") "Raw" else "Transformed"
  }
  x_label <- sprintf("%s %s Signal", axis_label_prefix, labs[ix])
  y_label <- sprintf("%s %s Signal", axis_label_prefix, labs[iy])
  
  tp <- if (exists("theme_Publication")) theme_Publication() else ggplot2::theme_bw(base_size = 12)
  has_ggpubr <- requireNamespace("ggpubr", quietly = TRUE)
  
  make_panel <- function(cc){
    d_plot <- df_plot[df_plot$condition == cc, , drop = FALSE]
    d_fit  <- stats_df[stats_df$condition == cc & is.finite(stats_df$x) & is.finite(stats_df$y), , drop = FALSE]
    
    p <- ggplot2::ggplot(d_plot, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(color = colors[[cc]], alpha = point_alpha, size = point_size, na.rm = TRUE) +
      ggplot2::scale_x_continuous(limits = x_limits) +
      ggplot2::scale_y_continuous(limits = y_limits) +
      ggplot2::labs(subtitle = as.character(cc), x = x_label, y = y_label) +
      ggplot2::theme(legend.position = "none") + tp +
      ggplot2::theme(legend.margin = ggplot2::margin(0,0,0,0))
    
    if (nrow(d_fit) >= 3) {
      p <- p + ggplot2::geom_smooth(data = d_fit, method = "lm", color = colors[[cc]], se = FALSE)
      if (has_ggpubr) {
        p <- p + ggpubr::stat_cor(data = d_fit, method = cor_method,
                                  label.x.npc = "left", label.y.npc = "top",
                                  size = 3, color = colors[[cc]])
      }
    }
    p
  }
  
  plots <- lapply(cond_levels, make_panel)
  
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(plots, nrow = 1) +
      patchwork::plot_annotation(
        title = sprintf("%s vs %s — assay '%s'%s",
                        labs[ix], labs[iy], assay_name,
                        if (transform != "none")
                          sprintf(" (%s%s)", transform,
                                  if (transform=="arcsinh") sprintf(", cofactor=%s", cofactor) else "")
                        else "")
      )
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    combined <- do.call(gridExtra::grid.arrange, c(plots, nrow = 1))
  } else {
    combined <- plots
  }
  
  list(plot = combined, stats = stat_tbl,
       params = list(x_feature = labs[ix], y_feature = labs[iy],
                     assay = assay_name, transform = transform, cofactor = cofactor,
                     cor_method = cor_method, outlier = outlier,
                     outlier_q = outlier_q, outlier_k = outlier_k))
}

compare_ptm_by_condition_phase <- function(
    se,
    features,
    assay_name = NULL,
    condition_col = "condition",
    phase_col = "cell_cycle",
    label_field = "antibody",
    transform_plot = c("none","log10","arcsinh"),
    cofactor = 150,
    sample_per_group = 20000,
    facet_ncol = 3,
    colors = NULL,
    add_stat_labels = TRUE,
    drop_empty_facets = TRUE,
    min_n_per_group = 1,
    min_n_per_condition = 25,
    facet_label_size = 9,
    plot_title = NULL
    
){
  # --- helper: union-bind safely ---
  bind_rows_safe <- function(lst){
    lst <- Filter(Negate(is.null), lst)
    if (!length(lst)) return(data.frame())
    if (requireNamespace("dplyr", quietly = TRUE)) return(dplyr::bind_rows(lst))
    if (requireNamespace("data.table", quietly = TRUE)) return(as.data.frame(data.table::rbindlist(lst, fill = TRUE)))
    cols <- Reduce(union, lapply(lst, names))
    lst2 <- lapply(lst, function(d){ miss <- setdiff(cols, names(d)); if (length(miss)) d[miss] <- NA; d[, cols, drop = FALSE] })
    do.call(rbind, lst2)
  }
  
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(condition_col %in% colnames(SummarizedExperiment::colData(se)))
  transform_plot <- match.arg(transform_plot)
  if (is.null(assay_name)) assay_name <- SummarizedExperiment::assayNames(se)[1]
  
  X   <- SummarizedExperiment::assay(se, assay_name)
  ids <- rownames(se)
  cd  <- as.data.frame(SummarizedExperiment::colData(se))
  cond <- as.character(cd[[condition_col]])
  if (!phase_col %in% names(cd)) stop("`phase_col` not found in colData(se).")
  phase_raw <- as.character(cd[[phase_col]])
  cond_levels <- unique(cond)
  
  rdat <- SummarizedExperiment::rowData(se)
  labs <- if (label_field %in% colnames(rdat)) as.character(rdat[, label_field]) else ids
  labs[is.na(labs) | labs==""] <- ids[is.na(labs) | labs==""]
  
  # resolve features
  resolve <- function(q){
    hit <- which(labs == q | ids == q)
    if (!length(hit)) hit <- which(tolower(labs) == tolower(q) | tolower(ids) == tolower(q))
    if (!length(hit)) hit <- unique(c(grep(q, labs, ignore.case=TRUE), grep(q, ids, ignore.case=TRUE)))
    if (!length(hit)) stop("Feature not found: ", q)
    hit[1]
  }
  f_idx <- vapply(features, resolve, 1L)
  f_lab <- labs[f_idx]
  
  # plotting transform (stats on assay scale)
  Tfun <- switch(transform_plot,
                 "none"    = function(v) v,
                 "log10"   = function(v){ z <- v; z[!is.finite(z) | z <= 0] <- NA_real_; log10(z) },
                 "arcsinh" = function(v) asinh(v / cofactor)
  )
  
  # long data: ALL + per-phase
  mat <- X[f_idx, , drop = FALSE]
  base_df <- data.frame(
    feature_id = rep(ids[f_idx], each = ncol(se)),
    feature    = rep(f_lab, each = ncol(se)),
    value      = as.numeric(t(mat)),
    condition  = rep(cond, times = length(f_idx)),
    phase      = rep(phase_raw, times = length(f_idx)),
    stringsAsFactors = FALSE
  )
  all_df <- base_df; all_df$phase <- "ALL"
  long <- rbind(all_df, base_df)
  long$plot_value <- Tfun(long$value)
  long$condition  <- factor(long$condition, levels = cond_levels)
  
  # omit sparse facets
  if (requireNamespace("dplyr", quietly = TRUE) && requireNamespace("tidyr", quietly = TRUE)) {
    counts <- dplyr::as_tibble(long) |>
      dplyr::filter(is.finite(value)) |>
      dplyr::count(feature, phase, condition, name = "n") |>
      tidyr::complete(feature, phase, condition = cond_levels, fill = list(n = 0)) |>
      dplyr::group_by(feature, phase) |>
      dplyr::summarise(ok = all(n >= min_n_per_condition), min_n = min(n), .groups = "drop")
    valid_facets <- counts |> dplyr::filter(ok) |> dplyr::select(feature, phase)
    dropped_facets <- counts |> dplyr::filter(!ok)
    long <- dplyr::semi_join(long, valid_facets, by = c("feature","phase"))
  } else {
    tt <- long[is.finite(long$value), c("feature","phase","condition")]
    key <- paste(tt$feature, tt$phase, sep = "\r")
    tab <- split(tt$condition, key)
    ok_map <- vapply(names(tab), function(k){
      v <- tab[[k]]
      n <- sapply(cond_levels, function(cl) sum(v == cl))
      all(n >= min_n_per_condition)
    }, logical(1))
    kp <- names(tab)[ok_map]
    kdf <- data.frame(feature = sub("\r.*","", kp), phase = sub(".*\r","", kp), stringsAsFactors = FALSE)
    long <- merge(long, kdf, by = c("feature","phase"))
    dropped_facets <- data.frame()
  }
  
  # downsample
  set.seed(1)
  long$grp <- interaction(long$feature, long$condition, long$phase, drop = TRUE)
  keep_idx <- unlist(lapply(split(seq_len(nrow(long)), long$grp), function(ix){
    if (is.finite(sample_per_group) && length(ix) > sample_per_group) sample(ix, sample_per_group) else ix
  }), use.names = FALSE)
  plot_df <- long[keep_idx, , drop = FALSE]
  
  # colors
  if (is.null(colors)) {
    base_cols <- c("#317ec2","#c03830","#2ca25f","#756bb1","#e6ab02","#1b9e77")
    colors <- setNames(base_cols[seq_along(cond_levels)], cond_levels)
  }
  
  # stats
  feat_phase_split <- split(long, interaction(long$feature, long$phase, drop = TRUE))
  stat_rows <- lapply(feat_phase_split, function(df){
    df <- df[is.finite(df$value) & !is.na(df$condition), , drop = FALSE]
    if (nrow(df) < 2) return(NULL)
    tab <- tapply(df$value, df$condition, function(v) sum(is.finite(v)))
    nonempty <- names(tab)[tab >= max(min_n_per_group, 1)]
    if (length(nonempty) < 2) {
      return(data.frame(feature = df$feature[1], phase = df$phase[1],
                        test = "none", p = NA_real_, p_adj = NA_real_,
                        note = "too few per group", stringsAsFactors = FALSE))
    }
    if (length(nonempty) == 2) {
      g1 <- nonempty[1]; g2 <- nonempty[2]
      x <- df$value[df$condition == g1]; x <- x[is.finite(x)]
      y <- df$value[df$condition == g2]; y <- y[is.finite(y)]
      wt <- suppressWarnings(stats::wilcox.test(x, y, exact = FALSE))
      n1 <- length(x); n2 <- length(y)
      W  <- as.numeric(wt$statistic); U <- W - n1*(n1+1)/2
      AUC <- U / (n1*n2); delta <- 2*AUC - 1
      return(data.frame(feature = df$feature[1], phase = df$phase[1],
                        test = "wilcox", p = wt$p.value, p_adj = NA_real_,
                        n1 = n1, n2 = n2, group1 = g1, group2 = g2,
                        effect = delta, stringsAsFactors = FALSE))
    }
    df2 <- df[df$condition %in% nonempty, , drop = FALSE]
    kt <- suppressWarnings(stats::kruskal.test(value ~ condition, data = df2))
    data.frame(feature = df$feature[1], phase = df$phase[1],
               test = "kruskal", p = kt$p.value, p_adj = NA_real_,
               stringsAsFactors = FALSE)
  })
  stats_tbl <- bind_rows_safe(stat_rows)
  if (nrow(stats_tbl)) {
    idx <- interaction(stats_tbl$feature, stats_tbl$phase, drop = TRUE)
    stats_tbl$p_adj <- ave(stats_tbl$p, idx, FUN = function(v) stats::p.adjust(v, "BH"))
  }
  
  # facet labels (p-values)
  facet_labels <- do.call(rbind, lapply(split(plot_df, interaction(plot_df$feature, plot_df$phase, drop = TRUE)), function(dfp){
    f <- dfp$feature[1]; ph <- dfp$phase[1]
    st <- stats_tbl[stats_tbl$feature == f & stats_tbl$phase == ph, , drop = FALSE]
    lab <- if (!nrow(st)) "p: —"
    else if (!("p_adj" %in% names(st)) || !is.finite(st$p_adj[1]) || st$test[1] == "none")
      (if ("note" %in% names(st)) st$note[1] else "p: —")
    else sprintf("p(BH)=%.2g", st$p_adj[1])
    y_pos <- suppressWarnings(stats::quantile(dfp$plot_value, 0.98, na.rm = TRUE))
    data.frame(feature = f, phase = ph, label = lab, y = y_pos)
  }))
  
  tp <- if (exists("theme_Publication")) theme_Publication(base_size = 14) else ggplot2::theme_bw(base_size = 14)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = condition, y = plot_value, fill = condition)) +
    ggplot2::geom_violin(scale = "width", trim = FALSE, color = NA) +
    ggplot2::geom_boxplot(width = 0.18, outlier.size = 0.3, alpha = 0.8) +
    ggplot2::facet_grid(feature ~ phase, scales = "free_y", drop = drop_empty_facets) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      title = plot_title,   # <-- use string or NULL
      x = "Condition",
      y = sprintf("Signal (%s%s)", transform_plot,
                  if (transform_plot=="arcsinh") sprintf(", c=%s", cofactor) else "")
    ) + tp +
    ggplot2::theme(
      legend.position = "none",
      strip.text.x = ggplot2::element_text(size = facet_label_size),
      strip.text.y = ggplot2::element_text(size = facet_label_size)
    )
  
  list(plot = p, stats = stats_tbl,
       dropped = if (exists("dropped_facets")) dropped_facets else NULL,
       params = list(features = f_lab, assay = assay_name,
                     transform_plot = transform_plot, cofactor = cofactor,
                     sample_per_group = sample_per_group,
                     min_n_per_group = min_n_per_group,
                     min_n_per_condition = min_n_per_condition,
                     facet_label_size = facet_label_size,
                     plot_title = plot_title))
}
