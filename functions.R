# set protein colors
proteins <- c("H1.5", "H1.4", "H2A1A", "H2A1B", "H2AJ", "H2AX", "H3.1", "H4")

# Generate a Viridis color palette with the same number of colors as there are proteins
colors <- viridis_pal(end = 0.9)(length(proteins))

# Assign each protein a color
protein_colors <- setNames(colors, proteins)

categorize_peptide <- function(peptide) {
  if (grepl("\\[un\\]", peptide) && !grepl("\\[(me|ac|me2|me3|su|hib)\\]", peptide)) {
    return("Unmodified")
  } else if (grepl("\\[me\\];?", peptide) && !grepl("\\[(ac|me2|me3|su|hib)\\]", peptide)) {
    return("Mono-methylated")
  } else if (grepl("\\[me2\\]", peptide) && !grepl("\\[(me|ac|me3|su|hib)\\]", peptide)) {
    return("Di-methylated")
  } else if (grepl("\\[me3\\]", peptide) && !grepl("\\[(me|ac|me2|su|hib)\\]", peptide)) {
    return("Tri-methylated")
  } else if (grepl("\\[ac\\]", peptide) && !grepl("\\[(me|me2|me3|su|hib)\\]", peptide)) {
    return("Acetylated")
  } else if (grepl("\\[(su|hib)\\]", peptide)) {
    return("Other modifications")
  } else {
    return("Hybrid modifications")
  }
}

# The input format looks like this: K9unK14un
format_peptide_note <- function(peptide_note) {
  
  # Use regex to match patterns like K9me3 and K14ac
  # The pattern looks for a "K" followed by digits (\d+), and then any letters (\w+), optionally followed by digits (\d*)
  parts <- regmatches(peptide_note, gregexpr("K\\d+[a-z]+\\d*", peptide_note, perl = TRUE))[[1]]
  
  # Apply formatting to insert brackets around the modification part, including any numbers that are part of the modification
  formatted_parts <- sapply(parts, function(part) {
    sub("(K\\d+)([a-z]+\\d*)", "\\1[\\2]", part)
  })
  
  # Concatenate the formatted parts with semicolons
  formatted_note <- paste(formatted_parts, collapse = ";")
  
  return(formatted_note)
}

protein_readable_translation <- function(protein_names) {
  # Extract the simplified name
  simple_names <- gsub(".*\\|(.+?)_.*", "\\1", protein_names)
  
  # Format specific names
  simple_names <- gsub("^H14$", "H1.4", simple_names)
  simple_names <- gsub("^H15$", "H1.5", simple_names)
  simple_names <- gsub("^H31$", "H3.1", simple_names)
  
  return(simple_names)
}

ptm_readable_translation <- function(ptms_vector) {
  # Function to process each individual PTM string
  process_ptm <- function(ptm) {
    # Split the string into protein part and modification parts
    parts <- unlist(strsplit(ptm, "-", fixed = TRUE))
    protein <- parts[1]
    mods <- parts[2]
    
    # Special handling for "H31" to be renamed as "H3.1"
    if (protein == "H31") {
      protein <- "H3"
    }
    
    # Split multiple modifications
    mods_split <- unlist(strsplit(mods, ";", fixed = TRUE))
    
    # Initialize a variable to hold the valid modifications
    valid_mods <- vector("character")
    
    # Filter and clean modifications
    for (mod in mods_split) {
      if (!grepl("\\[un\\]", mod)) {
        mod_clean <- gsub("\\[(.*)\\]", "\\1", mod)
        mod_clean <- gsub("(.*)\\[.*\\]", "\\1", mod_clean)
        valid_mods <- c(valid_mods, mod_clean)
      }
    }
    
    # Combine protein with valid modifications
    if (length(valid_mods) > 0) {
      return(paste0(protein, paste0(valid_mods, collapse = "")))
    } else {
      return(protein)  # Return just the protein name if no valid modifications
    }
  }
  
  # Apply the processing function to each element in the vector
  sapply(ptms_vector, process_ptm)
}


# Filters out peptides based on skyline annotations
# or if there are peptides with multiple precursors
# IMPROVEMENT: Combine information from multiple precursors
filter_peptides_skyline <- function(df){
  
  # Remove peptides with "" or "NOT FOUND" peptide note (just the unmod and no modifications)
  if (nrow(df[df$Peptide.Note == "" | 
              grepl("NOT FOUND", df$Peptide.Note) |
              grepl("NOTFOUND", df$Peptide.Note), ]) > 0) {
    
    print("Peptides with no PTMs")
    print(as.character(unique(df[df$Peptide.Note == "" | 
                                   df$Peptide.Note == "NOT FOUND" |
                                   df$Peptide.Note == "NOTFOUND", "Peptide.Sequence"])))
    
    # Remove " NOT FOUND" and any other string after this from Peptide.Note
    df$Peptide.Note <- sub(" NOT FOUND.*", "", df$Peptide.Note)
    
    # Remove " NOT FOUND" and any other string after this from Peptide.Note
    df$Peptide.Note <- sub("(NOTFOUND).*", "", df$Peptide.Note)
    
    # clean peptide note names
    df$Peptide.Note <- gsub(" ", "", df$Peptide.Note)
    
    df.filter.1 <- df
    
  } else {
    
    print("All peptides have PTMs")
    df.filter.1 <- df
    
  }
  
  # For peptides with multiple precursor charges, take the one with higher signal
  if (nrow(df.filter.1 %>% 
           group_by(Peptide.Note) %>% 
           summarise(n_distinct_charges = n_distinct(Precursor.Charge)) %>% 
           filter(n_distinct_charges > 1)) > 0) {
    filter.2 <- df.filter.1 %>%
      group_by(Peptide.Note) %>%
      filter(n_distinct(Precursor.Charge) > 1) %>%
      ungroup() %>%
      group_by(Peptide.Note, Precursor.Charge) %>%
      summarise(Avg.Total.Area.MS1 = mean(Total.Area.MS1, na.rm = TRUE)) %>%
      group_by(Peptide.Note) %>%
      filter(Avg.Total.Area.MS1 == min(Avg.Total.Area.MS1)) %>% # min as these are the ones to filter out
      ungroup()
    
    df.filter.2 <- df.filter.1 %>%
      anti_join(filter.2, by = c("Peptide.Note", "Precursor.Charge"))
    
    print("Peptidoforms with multiple precursor charge")
    print(unique(semi_join(df.filter.1, filter.2, by = c("Peptide.Note", "Precursor.Charge"))$Peptide.Note))
    
  } else {
    print("All PTMs have a single precursor charge")
  }
  
  return(df.filter.2)
}

filter_peptides_diann <- function(df){
  
  # Remove peptides with "" or "NOT FOUND" peptide note (just the unmod and no modifications)
  if (nrow(df[df$Peptide.Note == "" | 
              grepl("NOT FOUND", df$Peptide.Note) |
              grepl("NOTFOUND", df$Peptide.Note), ]) > 0) {
    
    print("Peptides with no PTMs")
    print(as.character(unique(df[df$Peptide.Note == "" | 
                                   df$Peptide.Note == "NOT FOUND" |
                                   df$Peptide.Note == "NOTFOUND", "Peptide.Sequence"])))
    
    # Remove " NOT FOUND" and any other string after this from Peptide.Note
    df.filter.1$Peptide.Note <- sub(" NOT FOUND.*", "", df.filter.1$Peptide.Note)
    
    # Remove " NOT FOUND" and any other string after this from Peptide.Note
    df.filter.1$Peptide.Note <- sub("(NOTFOUND).*", "", df.filter.1$Peptide.Note)
    
    # clean peptide note names
    df.filter.1$Peptide.Note <- gsub(" ", "", df.filter.1$Peptide.Note)
    
  } else {
    
    print("All peptides have PTMs")
    df.filter.1 <- df
    
  }
  
  # For peptides with multiple precursor charges, take the one with higher signal
  if (nrow(df.filter.1 %>% 
           group_by(Protein.Names) %>% 
           summarise(n_distinct_charges = n_distinct(Precursor.Charge)) %>% 
           filter(n_distinct_charges > 1)) > 0) {
    filter.2 <- df.filter.1 %>%
      group_by(Protein.Names) %>%
      filter(n_distinct(Precursor.Charge) > 1) %>%
      ungroup() %>%
      group_by(Protein.Names, Precursor.Charge) %>%
      summarise(Ms1.Area = mean(Ms1.Area, na.rm = TRUE)) %>%
      group_by(Protein.Names) %>%
      filter(Ms1.Area == min(Ms1.Area)) %>% # min as these are the ones to filter out
      ungroup()
    
    df.filter.2 <- df.filter.1 %>%
      anti_join(filter.2, by = c("Protein.Names", "Precursor.Charge"))
    
    print("Peptidoforms with multiple precursor charge")
    print(unique(semi_join(df.filter.1, filter.2, by = c("Protein.Names", "Precursor.Charge"))$Protein.Names))
    
  } else {
    print("All PTMs have a single precursor charge")
  }
  
  return(df.filter.2)
}

filter_peptides_diann_old <- function(df){
  
  # Remove peptides with "" or "NOT FOUND" peptide note (just the unmod and no modifications)
  if (nrow(df[df$Peptide.Note == "" | 
              grepl("NOT FOUND", df$Peptide.Note) |
              grepl("NOTFOUND", df$Peptide.Note), ]) > 0) {
    
    print("Peptides with no PTMs")
    print(as.character(unique(df[df$Peptide.Note == "" | 
                                 df$Peptide.Note == "NOT FOUND" |
                                 df$Peptide.Note == "NOTFOUND", "Peptide.Sequence"])))
    
  } else {
    
    print("All peptides have PTMs")
    df.filter.1 <- df
    
  }
  
  # Remove " NOT FOUND" and any other string after this from Peptide.Note
  df.filter.1$Peptide.Note <- sub(" NOT FOUND.*", "", df.filter.1$Peptide.Note)
  
  # clean peptide note names
  df.filter.1$Peptide.Note <- gsub(" ", "", df.filter.1$Peptide.Note)
  
  # For peptides with multiple precursor charges, take the one with higher signal
  if (nrow(df.filter.1 %>% 
           group_by(Peptide.Note) %>% 
           summarise(n_distinct_charges = n_distinct(Precursor.Charge)) %>% 
           filter(n_distinct_charges > 1)) > 0) {
    filter.2 <- df.filter.1 %>%
      group_by(Peptide.Note) %>%
      filter(n_distinct(Precursor.Charge) > 1) %>%
      ungroup() %>%
      group_by(Peptide.Note, Precursor.Charge) %>%
      summarise(Ms1.Area = mean(Ms1.Area, na.rm = TRUE)) %>%
      group_by(Peptide.Note) %>%
      filter(Ms1.Area == min(Ms1.Area)) %>% # min as these are the ones to filter out
      ungroup()
    
    df.filter.2 <- df.filter.1 %>%
      anti_join(filter.2, by = c("Peptide.Note", "Precursor.Charge"))
    
    print("Peptidoforms with multiple precursor charge")
    print(unique(semi_join(df.filter.1, filter.2, by = c("Peptide.Note", "Precursor.Charge"))$Peptide.Note))
    
  } else {
    print("All PTMs have a single precursor charge")
  }
  
  return(df.filter.2)
}


deconvolute_h4 <- function(df.combine.filter){
  
  # created nested dataframe for h4 fragment ions
  h4.fragment <- df.combine.filter %>%
    subset(!Fragment.Ion %in% c("precursor", "precursor [M+1]", "precursor [M+2]")) %>%
    subset(grepl("H4-K5", Peptide.Note)) %>%
    select(Replicate.Name, Peptide.Note, Fragment.Ion, Transition.Area) %>%
    group_by(Peptide.Note) %>% 
    nest()
  
  # create new MS1 intensities for H4 peptides
  
  # mono
  
  # peak 1
  
  # get K5
  # get sum for all fragments
  # get sum for unique fragments - b2, b3, b4
  mono.k5 <- h4.fragment[h4.fragment$Peptide.Note == "H4-K5[ac];K8[un];K12[un];K16[un]", ]$data[[1]]
  mono.k5 <- mono.k5 %>%
    group_by(Replicate.Name) %>%
    summarize(total.Transition.Area.fragment = sum(Transition.Area, na.rm = TRUE),
              unique.Transition.Area.fragment = sum(Transition.Area[Fragment.Ion %in% c("b2", "b3", "b4")], na.rm = TRUE))
  
  # get K8
  # get sum for all fragments
  # get sum for unique combo fragments - b5, b6, b7, b8
  mono.k8 <- h4.fragment[h4.fragment$Peptide.Note == "H4-K5[un];K8[ac];K12[un];K16[un]", ]$data[[1]]
  mono.k8 <- mono.k8 %>%
    group_by(Replicate.Name) %>%
    summarize(total.Transition.Area.fragment = sum(Transition.Area, na.rm = TRUE),
              unique.Transition.Area.fragment = sum(Transition.Area[Fragment.Ion %in% c("b5", "b6", "b7", "b8")], na.rm = TRUE))
  
  # get K12
  # get sum for all fragments
  # get sum for unique combo fragments - y9, y8, y7, y6
  mono.k12 <- h4.fragment[h4.fragment$Peptide.Note == "H4-K5[un];K8[un];K12[ac];K16[un]", ]$data[[1]]
  mono.k12 <- mono.k12 %>%
    group_by(Replicate.Name) %>%
    summarize(total.Transition.Area.fragment = sum(Transition.Area, na.rm = TRUE),
              unique.Transition.Area.fragment = sum(Transition.Area[Fragment.Ion %in% c("y9", "y8", "y7", "y6")], na.rm = TRUE))
  
  # first get the proportion of unique fragments
  # K5 prop unique = K5 (unique) / all fragments
  mono.peak1 <- mono.k5 %>%
    mutate(k5.prop = unique.Transition.Area.fragment/total.Transition.Area.fragment) %>%
    select(Replicate.Name, k5.prop)
  
  # this is what we need to doeonvolute
  # K8 + K12 prop = 1 - K5 prop
  mono.peak1 <- mono.peak1 %>%
    mutate(k8.k12.prop = 1 - k5.prop)
  
  # second get proportion of combo unique frags
  # K8 prop combo = K8 (combo unique) / all fragments
  mono.k8 <- mono.k8 %>%
    mutate(prop.combo = unique.Transition.Area.fragment/total.Transition.Area.fragment)
  
  # K12 prop combo = K12 (combo unique) / all fragments
  mono.k12 <- mono.k12 %>%
    mutate(prop.combo = unique.Transition.Area.fragment/total.Transition.Area.fragment)
  
  # third get how their combo prop relates to the total proportion
  # K8 prop rel = K8 prop combo / (K8 prop combo + K12 prop combo)
  mono.k8 <- mono.k8 %>%
    mutate(k8.prop.rel = prop.combo / (prop.combo + mono.k12$prop.combo))
  
  # K12 prop rel = K12 prop combo / (K8 prop combo + K12 prop combo)
  mono.k12 <- mono.k12 %>%
    mutate(k12.prop.rel = prop.combo / (prop.combo + mono.k8$prop.combo))
  
  # fourth use this to deconvolute the combined proportion
  # K8 prop unique = K8 rel prop * K8+K12 prop
  mono.peak1 <- mono.peak1 %>%
    mutate(k8.prop = mono.k8$k8.prop.rel * k8.k12.prop)
  
  # K12 prop unique = K12 rel prop * K8+K12 prop
  mono.peak1 <- mono.peak1 %>%
    mutate(k12.prop = mono.k12$k12.prop.rel * k8.k12.prop)
  
  # fifth check if K5 prop + K8 prop + K12 prop = 1
  mono.peak1 <- mono.peak1 %>%
    select(-k8.k12.prop)
  
  rowSums(mono.peak1[,-1])
  
  # finally update MS1
  df.combine.filter.mono <- merge(df.combine.filter, mono.peak1, by = "Replicate.Name")
  
  df.combine.filter.mono <- df.combine.filter.mono %>%
    mutate(Total.Area.MS1 = case_when(
      Peptide.Note == "H4-K5[ac];K8[un];K12[un];K16[un]" ~ Total.Area.MS1 * k5.prop,
      Peptide.Note == "H4-K5[un];K8[ac];K12[un];K16[un]" ~ Total.Area.MS1 * k8.prop,
      Peptide.Note == "H4-K5[un];K8[un];K12[ac];K16[un]" ~ Total.Area.MS1 * k12.prop,
      TRUE ~ Total.Area.MS1 # Default case 
    ))
  
  # di
  
  # peak 1
  
  # K5K8
  # get sum for all fragments
  # get sum for unique fragments - b5, b6, b7, b8
  di.k5k8 <- h4.fragment[h4.fragment$Peptide.Note == "H4-K5[ac];K8[ac];K12[un];K16[un]", ]$data[[1]]
  di.k5k8 <- di.k5k8 %>%
    group_by(Replicate.Name) %>%
    summarize(total.Transition.Area.fragment = sum(Transition.Area, na.rm = TRUE),
              unique.Transition.Area.fragment = sum(Transition.Area[Fragment.Ion %in% c("b5", "b6", "b7", "b8")], na.rm = TRUE))
  
  # K5K12
  # get sum for all fragments
  # get sum for combo unique fragments - b2, b3, b4
  di.k5k12 <- h4.fragment[h4.fragment$Peptide.Note == "H4-K5[ac];K8[un];K12[ac];K16[un]", ]$data[[1]]
  di.k5k12 <- di.k5k12 %>%
    group_by(Replicate.Name) %>%
    summarize(total.Transition.Area.fragment = sum(Transition.Area, na.rm = TRUE),
              unique.Transition.Area.fragment = sum(Transition.Area[Fragment.Ion %in% c("b2", "b3", "b4")], na.rm = TRUE))
  
  # K8K12
  # get sum for all fragments
  # get sum for combo unique fragments - b2, b3, b4
  di.k8k12 <- h4.fragment[h4.fragment$Peptide.Note == "H4-K5[un];K8[ac];K12[ac];K16[un]", ]$data[[1]]
  di.k8k12 <- di.k8k12 %>%
    group_by(Replicate.Name) %>%
    summarize(total.Transition.Area.fragment = sum(Transition.Area, na.rm = TRUE),
              unique.Transition.Area.fragment = sum(Transition.Area[Fragment.Ion %in% c("b2", "b3", "b4")], na.rm = TRUE))
  
  # first get the proportion of unique fragments
  # K5K8 prop unique = K5K8 (unique) / all fragments
  di.peak1 <- di.k5k8 %>%
    mutate(k5k8.prop = unique.Transition.Area.fragment/total.Transition.Area.fragment) %>%
    select(Replicate.Name, k5k8.prop)
  
  # K5K12 + K8K12 prop = 1 - K5K8 props
  di.peak1 <- di.peak1 %>%
    mutate(k5k12.k8k12.prop = 1 - k5k8.prop)
  
  # second get proportion of combo unique frags
  # K5k12 prop combo = K5k12 (combo unique) / all fragments
  di.k5k12 <- di.k5k12 %>%
    mutate(prop.combo = unique.Transition.Area.fragment/total.Transition.Area.fragment)
  
  # K8k12 prop combo = K8k12 (combo unique) / all fragments
  di.k8k12 <- di.k8k12 %>%
    mutate(prop.combo = unique.Transition.Area.fragment/total.Transition.Area.fragment)
  
  # third get how their combo prop relates to the total proportion
  # K8k12 prop rel = K8k12 prop combo / (K8k12 prop combo + k8k12 prop combo)
  di.k5k12 <- di.k5k12 %>%
    mutate(k5k12.prop.rel = prop.combo / (prop.combo + di.k8k12$prop.combo))
  
  # k8k12 prop rel = k8k12 prop combo / (k8k12 prop combo + k5k12 prop combo)
  di.k8k12 <- di.k8k12 %>%
    mutate(k8k12.prop.rel = prop.combo / (prop.combo + di.k5k12$prop.combo))
  
  # fourth use this to deconvolute the combined proportion
  # k5k12 prop unique = k5k12 rel prop * K5K12 + K8K12 prop
  di.peak1 <- di.peak1 %>%
    mutate(k5k12.prop = di.k5k12$k5k12.prop.rel * k5k12.k8k12.prop)
  
  # k8k12 prop unique = k8k12 rel prop * K5K12 + K8K12 prop
  di.peak1 <- di.peak1 %>%
    mutate(k8k12.prop = di.k8k12$k8k12.prop.rel * k5k12.k8k12.prop)
  
  # fifth check if K5 prop + K8 prop + K12 prop = 1
  di.peak1 <- di.peak1 %>%
    select(-k5k12.k8k12.prop)
  
  rowSums(di.peak1[,-1])
  
  # finally update MS1
  df.combine.filter.mono.di1 <- merge(df.combine.filter.mono, di.peak1, by = "Replicate.Name")
  
  df.combine.filter.mono.di1 <- df.combine.filter.mono.di1 %>%
    mutate(Total.Area.MS1 = case_when(
      Peptide.Note == "H4-K5[ac];K8[ac];K12[un];K16[un]" ~ Total.Area.MS1 * k5k8.prop,
      Peptide.Note == "H4-K5[ac];K8[un];K12[ac];K16[un]" ~ Total.Area.MS1 * k5k12.prop,
      Peptide.Note == "H4-K5[un];K8[ac];K12[ac];K16[un]" ~ Total.Area.MS1 * k8k12.prop,
      TRUE ~ Total.Area.MS1 # Default case 
    ))
  
  # peak 2
  
  # k12k16
  # get sum for all fragments
  # get sum for unique fragments - y9, y8, y7, y6
  di.k12k16 <- h4.fragment[h4.fragment$Peptide.Note == "H4-K5[un];K8[un];K12[ac];K16[ac]", ]$data[[1]]
  di.k12k16 <- di.k12k16 %>%
    group_by(Replicate.Name) %>%
    summarize(total.Transition.Area.fragment = sum(Transition.Area, na.rm = TRUE),
              unique.Transition.Area.fragment = sum(Transition.Area[Fragment.Ion %in% c("y9", "y8", "y7", "y6")], na.rm = TRUE))
  
  # k8k16
  # get sum for all fragments
  # get sum for combo unique fragments - b2, b3, b4
  di.k8k16 <- h4.fragment[h4.fragment$Peptide.Note == "H4-K5[un];K8[ac];K12[un];K16[ac]", ]$data[[1]]
  di.k8k16 <- di.k8k16 %>%
    group_by(Replicate.Name) %>%
    summarize(total.Transition.Area.fragment = sum(Transition.Area, na.rm = TRUE),
              unique.Transition.Area.fragment = sum(Transition.Area[Fragment.Ion %in% c("b2", "b3", "b4")], na.rm = TRUE))
  
  # k5k16
  # get sum for all fragments
  # get sum for combo unique fragments - b2, b3, b4
  di.k5k16 <- h4.fragment[h4.fragment$Peptide.Note == "H4-K5[ac];K8[un];K12[un];K16[ac]", ]$data[[1]]
  di.k5k16 <- di.k5k16 %>%
    group_by(Replicate.Name) %>%
    summarize(total.Transition.Area.fragment = sum(Transition.Area, na.rm = TRUE),
              unique.Transition.Area.fragment = sum(Transition.Area[Fragment.Ion %in% c("b2", "b3", "b4")], na.rm = TRUE))
  
  # first get the proportion of unique fragments
  # k12k16 prop unique = k12k16 (unique) / all fragments
  di.peak2 <- di.k12k16 %>%
    mutate(k12k16.prop = unique.Transition.Area.fragment/total.Transition.Area.fragment) %>%
    select(Replicate.Name, k12k16.prop)
  
  # k8k16 + k5k16 prop = 1 - k12k16 props
  di.peak2 <- di.peak2 %>%
    mutate(k8k16.k5k16.prop = 1 - k12k16.prop)
  
  # second get proportion of combo unique frags
  # k8k16 prop combo = k8k16 (combo unique) / all fragments
  di.k8k16 <- di.k8k16 %>%
    mutate(prop.combo = unique.Transition.Area.fragment/total.Transition.Area.fragment)
  
  # k5k16 prop combo = k5k16 (combo unique) / all fragments
  di.k5k16 <- di.k5k16 %>%
    mutate(prop.combo = unique.Transition.Area.fragment/total.Transition.Area.fragment)
  
  # third get how their combo prop relates to the total proportion
  di.k8k16 <- di.k8k16 %>%
    mutate(k8k16.prop.rel = prop.combo / (prop.combo + di.k5k16$prop.combo))
  
  di.k5k16 <- di.k5k16 %>%
    mutate(k5k16.prop.rel = prop.combo / (prop.combo + di.k8k16$prop.combo))
  
  # fourth use this to deconvolute the combined proportion
  di.peak2 <- di.peak2 %>%
    mutate(k8k16.prop = di.k8k16$k8k16.prop.rel * k8k16.k5k16.prop)
  
  di.peak2 <- di.peak2 %>%
    mutate(k5k16.prop = di.k5k16$k5k16.prop.rel * k8k16.k5k16.prop)
  
  # fifth check if summed prop = 1
  di.peak2 <- di.peak2 %>%
    select(-k8k16.k5k16.prop)
  
  rowSums(di.peak2[,-1])
  
  # finally update MS1
  df.combine.filter.mono.di1.di2 <- merge(df.combine.filter.mono.di1, di.peak2, by = "Replicate.Name")
  
  df.combine.filter.mono.di1.di2 <- df.combine.filter.mono.di1.di2 %>%
    mutate(Total.Area.MS1 = case_when(
      Peptide.Note == "H4-K5[ac];K8[un];K12[un];K16[ac]" ~ Total.Area.MS1 * k5k16.prop,
      Peptide.Note == "H4-K5[un];K8[ac];K12[un];K16[ac]" ~ Total.Area.MS1 * k8k16.prop,
      Peptide.Note == "H4-K5[un];K8[un];K12[ac];K16[ac]" ~ Total.Area.MS1 * k12k16.prop,
      TRUE ~ Total.Area.MS1 # Default case 
    ))
  
  # tri
  
  # peak 1
  
  # k8k12k16
  # get sum for all fragments
  # get sum for unique fragments - b2, b3, b4
  tri.k8k12k16 <- h4.fragment[h4.fragment$Peptide.Note == "H4-K5[un];K8[ac];K12[ac];K16[ac]", ]$data[[1]]
  tri.k8k12k16 <- tri.k8k12k16 %>%
    group_by(Replicate.Name) %>%
    summarize(total.Transition.Area.fragment = sum(Transition.Area, na.rm = TRUE),
              unique.Transition.Area.fragment = sum(Transition.Area[Fragment.Ion %in% c("b2", "b3", "b4")], na.rm = TRUE))
  
  # k5k12k16
  # get sum for all fragments
  # get sum for combo unique fragments - b5, b6, b7, b8, y9, y8, y8, y6
  tri.k5k12k16 <- h4.fragment[h4.fragment$Peptide.Note == "H4-K5[ac];K8[un];K12[ac];K16[ac]", ]$data[[1]]
  tri.k5k12k16 <- tri.k5k12k16 %>%
    group_by(Replicate.Name) %>%
    summarize(total.Transition.Area.fragment = sum(Transition.Area, na.rm = TRUE),
              unique.Transition.Area.fragment = sum(Transition.Area[Fragment.Ion %in% c("b5", "b6", "b7", "b8", "y9", "y8", "y8", "y6")], na.rm = TRUE))
  
  # k5k8k16
  # get sum for all fragments
  # get sum for combo unique fragments - b5, b6, b7, b8, y9, y8, y8, y6
  tri.k5k8k16 <- h4.fragment[h4.fragment$Peptide.Note == "H4-K5[ac];K8[ac];K12[un];K16[ac]", ]$data[[1]]
  tri.k5k8k16 <- tri.k5k8k16 %>%
    group_by(Replicate.Name) %>%
    summarize(total.Transition.Area.fragment = sum(Transition.Area, na.rm = TRUE),
              unique.Transition.Area.fragment = sum(Transition.Area[Fragment.Ion %in% c("b5", "b6", "b7", "b8", "y9", "y8", "y8", "y6")], na.rm = TRUE))
  
  # first get the proportion of unique fragments
  # k8k12k16 prop unique = k8k12k16 (unique) / all fragments
  tri.peak1 <- tri.k8k12k16 %>%
    mutate(k8k12k16.prop = unique.Transition.Area.fragment/total.Transition.Area.fragment) %>%
    select(Replicate.Name, k8k12k16.prop)
  
  # k5k12k16 + K8K12 prop = 1 - k8k12k16 props
  tri.peak1 <- tri.peak1 %>%
    mutate(k5k12k16.k5k8k16.prop = 1 - k8k12k16.prop)
  
  # second get proportion of combo unique frags
  # k5k12k16 prop combo = k5k12k16 (combo unique) / all fragments
  tri.k5k12k16 <- tri.k5k12k16 %>%
    mutate(prop.combo = unique.Transition.Area.fragment/total.Transition.Area.fragment)
  
  # k5k8k16 prop combo = k5k8k16 (combo unique) / all fragments
  tri.k5k8k16 <- tri.k5k8k16 %>%
    mutate(prop.combo = unique.Transition.Area.fragment/total.Transition.Area.fragment)
  
  # third get how their combo prop relates to the total proportion
  tri.k5k12k16 <- tri.k5k12k16 %>%
    mutate(k5k12k16.prop.rel = prop.combo / (prop.combo + tri.k5k8k16$prop.combo))
  
  tri.k5k8k16 <- tri.k5k8k16 %>%
    mutate(k5k8k16.prop.rel = prop.combo / (prop.combo + tri.k5k12k16$prop.combo))
  
  # fourth use this to deconvolute the combined proportion
  tri.peak1 <- tri.peak1 %>%
    mutate(k5k12k16.prop = tri.k5k12k16$k5k12k16.prop.rel * k5k12k16.k5k8k16.prop)
  
  # K12 prop unique = K12 rel prop * K8+K12 prop
  tri.peak1 <- tri.peak1 %>%
    mutate(k5k8k16.prop = tri.k5k8k16$k5k8k16.prop.rel * k5k12k16.k5k8k16.prop)
  
  # fifth check if summed prop = 1
  tri.peak1 <- tri.peak1 %>%
    select(-k5k12k16.k5k8k16.prop)
  
  rowSums(tri.peak1[,-1])
  
  # finally update MS1
  df.combine.filter.mono.di1.di2.tri <- merge(df.combine.filter.mono.di1.di2, tri.peak1, by = "Replicate.Name")
  
  df.combine.filter.mono.di1.di2.tri <- df.combine.filter.mono.di1.di2.tri %>%
    mutate(Total.Area.MS1 = case_when(
      Peptide.Note == "H4-K5[un];K8[ac];K12[ac];K16[ac]" ~ Total.Area.MS1 * k8k12k16.prop,
      Peptide.Note == "H4-K5[ac];K8[un];K12[ac];K16[ac]" ~ Total.Area.MS1 * k5k12k16.prop,
      Peptide.Note == "H4-K5[ac];K8[ac];K12[un];K16[ac]" ~ Total.Area.MS1 * k5k8k16.prop,
      TRUE ~ Total.Area.MS1 # Default case 
    ))
  
  # remove all the proportions 
  # and return the same columns that were used in input
  return(df.combine.filter.mono.di1.di2.tri %>%
           select(colnames(df.combine.filter))
  )
}

# df <- df.combine.filter
# coldata <- coldata
load_summarized_experiment_diann <- function(df, coldata){
  
  # extract assay data
  
    # protein group MS2 based quantification
    df.quant.pg <- df %>%
      select(Peptide.Note, Run, PG.Quantity) %>%
      spread(key = Run, value = PG.Quantity)
    row.names(df.quant.pg) <- df.quant.pg$Peptide.Note
    df.quant.pg <- df.quant.pg[,-1]
    
    # precursor MS2 based quantification
    df.quant.pr <- df %>%
      select(Peptide.Note, Run, Precursor.Quantity) %>%
      spread(key = Run, value = Precursor.Quantity)
    row.names(df.quant.pr) <- df.quant.pr$Peptide.Note
    df.quant.pr <- df.quant.pr[,-1]
    
    # precursor MS1 total area
    df.area.ms1 <- df %>%
      select(Peptide.Note, Run, Ms1.Area) %>%
      spread(key = Run, value = Ms1.Area)
    row.names(df.area.ms1) <- df.area.ms1$Peptide.Note
    df.area.ms1 <- df.area.ms1[,-1]
    
    # MS2 quant split by fragments
    split_to_vector <- function(cell_value) {
      # Remove the trailing semicolon if it exists
      cell_value <- gsub(";$", "", cell_value)
      # Split the string into individual numeric values and return as a vector
      as.numeric(unlist(strsplit(cell_value, ";")))
    }
    df.quant.ms2 <- df %>%
      select(Peptide.Note, Run, Fragment.Quant.Raw) %>%
      mutate(Fragment.Quant.Raw.Split = sapply(Fragment.Quant.Raw, split_to_vector)) %>%
      select(Peptide.Note, Run, Fragment.Quant.Raw.Split) %>%
      spread(key = Run, value = Fragment.Quant.Raw.Split)
    row.names(df.quant.ms2) <- df.quant.ms2$Peptide.Note
    df.quant.ms2 <- df.quant.ms2[,-1]
    
    # MS2 quant sum
    sum_cell_values <- function(cell_value) {
      # Remove the trailing semicolon if it exists
      cell_value <- gsub(";$", "", cell_value)
      # Split the string into individual numeric values
      numbers <- as.numeric(unlist(strsplit(cell_value, ";")))
      # Sum the values
      sum(numbers, na.rm = TRUE)
    }
    df.quant.ms2.sum <- df %>%
      select(Peptide.Note, Run, Fragment.Quant.Raw) %>%
      mutate(Fragment.Quant.Raw.Sum = sapply(Fragment.Quant.Raw, sum_cell_values)) %>%
      select(Peptide.Note, Run, Fragment.Quant.Raw.Sum) %>%
      spread(key = Run, value = Fragment.Quant.Raw.Sum)
    row.names(df.quant.ms2.sum) <- df.quant.ms2.sum$Peptide.Note
    df.quant.ms2.sum <- df.quant.ms2.sum[,-1]
    
    # retention time
    df.rt <- df %>%
      select(Peptide.Note, Run, RT) %>%
      spread(key = Run, value = RT)
    row.names(df.rt) <- df.rt$Peptide.Note
    df.rt <- df.rt[,-1]
    
    # ion mobility
    df.im <- df %>%
      select(Peptide.Note, Run, IM) %>%
      spread(key = Run, value = IM)
    row.names(df.im) <- df.im$Peptide.Note
    df.im <- df.im[,-1]
  
    # q value
    df.qvalue <- df %>%
      select(Peptide.Note, Run, Q.Value) %>%
      spread(key = Run, value = Q.Value)
    row.names(df.qvalue) <- df.qvalue$Peptide.Note
    df.qvalue <- df.qvalue[,-1]
  
  # extract row data
  
    # fragment info
    # y4-unknown^1/517.3092651;y5-unknown^1/701.4304199;y3-unknown^1/446.2721558;b3-unknown^1/383.2288818;b5-unknown^1/654.3820801;
    # unsure what the numeric values here are
    extract_before_hyphen <- function(input_string) {
      # Remove any trailing semicolon
      input_string <- gsub(";$", "", input_string)
      # Split the string by semicolons
      split_string <- unlist(strsplit(input_string, ";"))
      # Extract the part before the hyphen using regular expressions
      result <- gsub("-.*", "", split_string)
      # Return unique values
      unique(result)
    }
    df.frag.info <- df %>%
      select(Peptide.Note, Fragment.Info) %>%
      distinct() %>%
      mutate(Fragments = sapply(Fragment.Info, extract_before_hyphen)) %>%
      select(Peptide.Note, Fragments)
    
    # precursor charge
    df.charge.pr <- df %>%
      select(Peptide.Note, Precursor.Charge) %>%
      distinct()
  
    # global qvalue
    df.qvalue.global <- df %>%
      select(Peptide.Note, Global.Q.Value) %>%
      distinct()
    
    # protein name
    
    # uniprot
    
    # modified sequence
    df.mod.seq <- df %>%
      select(Peptide.Note, Modified.Sequence) %>%
      distinct()
    
    # peptide sequence
    df.seq <- df %>%
      select(Peptide.Note, Stripped.Sequence) %>%
      distinct()
    
    # combine
    rowdata <- cbind(df.mod.seq,
                     df.seq,
                     df.charge.pr,
                     df.qvalue.global,
                     df.frag.info)
    row.names(rowdata) <- rowdata$Peptide.Note
    rowdata <- rowdata %>%
      select(-starts_with("Peptide.Note"))
    rowdata <- rowdata[match(rownames(df.quant.pg), rownames(rowdata)),]
    
  # filter coldata by unique run id
  # Extract the part after the last underscore from colnames of df.ms1
  extracted_run_id <- sapply(strsplit(colnames(df.quant.pg), "_", fixed = TRUE), function(x) tail(x, 1))
  coldata.filter <- coldata %>%
    filter(run_id %in% extracted_run_id)
  
  coldata.filter <- coldata.filter[match(extracted_run_id, coldata.filter$run_id),]
  
  # create summarized exp obj
  se <- SummarizedExperiment(assays = list(
    protein.quant = df.quant.pg,
    precursor.quant = df.quant.pr,
    ms1.area = df.area.ms1,
    ms2.quant = df.quant.ms2,
    ms2.quant.sum = df.quant.ms2.sum,
    retention.time = df.rt,
    ion.mobility = df.im,
    qvalue = df.qvalue),
    rowData = rowdata,
    colData = coldata.filter)
  
  return(se)
}

# Note that this does not change NA values to 0
# IMPROVEMENT: Be able to hold fragment ion data mapping 
# multiple fragment ions to a single peptide
load_summarized_experiment_skyline <- function(df, coldata){
  
  # remove spaces from strings
  df[] <- lapply(df, function(x) {
    if (is.character(x)) {
      gsub(" ", "", x) # Remove spaces
    } else {
      x # Leave numeric columns unchanged
    }
  })
  
  coldata[] <- lapply(coldata, function(x) {
    if (is.character(x)) {
      gsub(" ", "", x) # Remove spaces
    } else {
      x # Leave numeric columns unchanged
    }
  })
  
  # ms1 total area
  df.ms1 <- df %>%
    filter(Fragment.Ion == "precursor") %>%
    select(Peptide.Note, Replicate.Name, Total.Area.MS1) %>%
    spread(key = Replicate.Name, value = Total.Area.MS1)
  row.names(df.ms1) <- df.ms1$Peptide.Note
  df.ms1 <- df.ms1[,-1]
  
  # ms2 total area
  df.ms2 <- df %>%
    filter(Fragment.Ion == "precursor") %>%
    select(Peptide.Note, Replicate.Name, Total.Area.Fragment) %>%
    spread(key = Replicate.Name, value = Total.Area.Fragment)
  row.names(df.ms2) <- df.ms2$Peptide.Note
  df.ms2 <- df.ms2[,-1]
  
  # retention time
  df.rt <- df %>%
    subset(Fragment.Ion == "precursor") %>%
    select(Peptide.Note, Replicate.Name, Peptide.Retention.Time) %>%
    spread(key = Replicate.Name, value = Peptide.Retention.Time)
  row.names(df.rt) <- df.rt$Peptide.Note
  df.rt <- df.rt[,-1]
  
  # extract row data
  rowdata.pep <- df %>%
    subset(Fragment.Ion == "precursor") %>%
    select(Protein.Accession, 
           Protein.Gene, 
           Peptide.Sequence, 
           Peptide.Note, 
           Precursor.Charge,
           PrecursorMz) %>%
    distinct(Peptide.Note, .keep_all = TRUE) %>%
    as.data.frame()
  rownames(rowdata.pep) <- rowdata.pep$Peptide.Note
  rowdata.pep$Peptide.Note <- NULL
  
  rowdata.frag <- df %>%
    subset(!Fragment.Ion %in% c("precursor", "precursor [M+1]", "precursor [M+2]")) %>%
    select(Peptide.Note,
           Fragment.Ion,
           ProductMz,
           Product.Charge,
           Transition.Area,
           ProductMz) %>%
    group_by(Peptide.Note) %>%
    summarise(Fragment.Ion = list(Fragment.Ion),
              ProductMz = list(ProductMz),
              Product.Charge = list(Product.Charge),
              Transition.Area = list(Transition.Area)) %>%
    as.data.frame()
    rownames(rowdata.frag) <- rowdata.frag$Peptide.Note
    rowdata.frag$Peptide.Note <- NULL
    
    rowdata <- cbind(rowdata.pep, rowdata.frag)
    
    # cleaned histone peptidoform names
    rowdata$Peptide.Note.Clean <- gsub("\\[|\\]|;|-", "", rownames(rowdata))
    
    rowdata <- rowdata[match(rownames(df.ms1), rownames(rowdata)),]
  
  # filter coldata by unique run id
  # Extract the part after the last underscore from colnames of df.ms1
  extracted_run_id <- sapply(strsplit(colnames(df.ms1), "_", fixed = TRUE), function(x) tail(x, 1))
  
  # Filter colData
  # add sum of MS1 and MS2
  coldata.filter <- coldata %>%
    filter(run_id %in% extracted_run_id) %>%
    mutate(total.ms1 = colSums(df.ms1, na.rm = TRUE),
           total.ms2 = colSums(df.ms2, na.rm = TRUE))
  
  coldata.filter <- coldata.filter[match(extracted_run_id, coldata.filter$run_id),]
  
  # create summarized exp obj
  se <- SummarizedExperiment(assays = list(MS1 = df.ms1,
                                           MS2 = df.ms2,
                                           RT = df.rt),
                             rowData = rowdata,
                             colData = coldata.filter)
  
  return(se)
}

load_summarized_experiment_prop_eff_skyline <- function(df, coldata){
  
  # extract assay data
  
  # ms1 total area
  df.ms1 <- df %>%
    filter(Fragment.Ion == "precursor") %>%
    select(Peptide.Note, Replicate.Name, Total.Area.MS1) %>%
    spread(key = Replicate.Name, value = Total.Area.MS1)
  row.names(df.ms1) <- df.ms1$Peptide.Note
  df.ms1 <- df.ms1[,-1]
  
  # ms2 total area
  df.ms2 <- df %>%
    filter(Fragment.Ion == "precursor") %>%
    select(Peptide.Note, Replicate.Name, Total.Area.Fragment) %>%
    spread(key = Replicate.Name, value = Total.Area.Fragment)
  row.names(df.ms2) <- df.ms2$Peptide.Note
  df.ms2 <- df.ms2[,-1]
  
  # retention time
  df.rt <- df %>%
    filter(Fragment.Ion == "precursor") %>%
    select(Peptide.Note, Replicate.Name, Peptide.Retention.Time) %>%
    spread(key = Replicate.Name, value = Peptide.Retention.Time)
  row.names(df.rt) <- df.rt$Peptide.Note
  df.rt <- df.rt[,-1]
  
  # extract row data
  rowdata <- df %>%
    select(Protein.Name, Peptide.Modified.Sequence.Monoisotopic.Masses, Peptide.Note, Precursor.Charge) %>%
    unique() %>%
    group_by(Protein.Name, Peptide.Modified.Sequence.Monoisotopic.Masses, Peptide.Note, Precursor.Charge) %>%
    as.data.frame()
  rownames(rowdata) <- rowdata$Peptide.Note
  rowdata$Peptide.Note <- NULL
  rowdata <- rowdata[match(rownames(df.ms1), rownames(rowdata)),]
  
  # filter coldata by unique run id
  # Extract the part after the last underscore from colnames of df.ms1
  extracted_run_id <- sapply(strsplit(colnames(df.ms1), "_", fixed = TRUE), function(x) tail(x, 1))
  
  # Filter colData
  coldata.filter <- coldata %>%
    filter(run_id %in% extracted_run_id)
  
  coldata.filter <- coldata.filter[match(extracted_run_id, coldata.filter$run_id),]
  
  # create summarized exp obj
  se <- SummarizedExperiment(assays = list(MS1 = df.ms1,
                                           MS2 = df.ms2,
                                           RT = df.rt),
                             rowData = rowdata,
                             colData = coldata.filter)
  
  return(se)
}

load_summarized_experiment_ratio <- function(ratio.ms1, ratio.ms2, obj.filter){
  
  # extract assay data
  
  # ms1 total area
  df.ms1 <- ratio.ms1 %>%
    select(Peptide.Note, Replicate.Name, ratio) %>%
    spread(key = Replicate.Name, value = ratio)
  row.names(df.ms1) <- df.ms1$Peptide.Note
  df.ms1 <- df.ms1[,-1]
  
  df.ms1.log1p <- ratio.ms1 %>%
    select(Peptide.Note, Replicate.Name, ratio.log1p) %>%
    spread(key = Replicate.Name, value = ratio.log1p)
  row.names(df.ms1.log1p) <- df.ms1.log1p$Peptide.Note
  df.ms1.log1p <- df.ms1.log1p[,-1]
  
  # ms2 total area
  df.ms2 <- ratio.ms2 %>%
    select(Peptide.Note, Replicate.Name, ratio) %>%
    spread(key = Replicate.Name, value = ratio)
  row.names(df.ms2) <- df.ms2$Peptide.Note
  df.ms2 <- df.ms2[,-1]
  
  df.ms2.log1p <- ratio.ms2 %>%
    select(Peptide.Note, Replicate.Name, ratio.log1p) %>%
    spread(key = Replicate.Name, value = ratio.log1p)
  row.names(df.ms2.log1p) <- df.ms2.log1p$Peptide.Note
  df.ms2.log1p <- df.ms2.log1p[,-1]
  
  # retention time
  df.rt <- assay(obj.filter, "RT")
  df.rt <- df.rt[rownames(df.rt) %in% rownames(df.ms1),]
  df.rt <- df.rt[match(rownames(df.ms1), rownames(df.rt)),]
  df.rt <- df.rt[,colnames(df.rt) %in% colnames(df.ms1)]
  df.rt <- df.rt[,match(colnames(df.ms1), colnames(df.rt))]
  
  
  # add peptide information
  rowdata <- rowData(obj.filter)
  rowdata <- rowdata[rownames(rowdata) %in% rownames(df.ms1),]
  rowdata <- rowdata[match(rownames(df.ms1), rownames(rowdata)),]
  
  # add meta data
  # filter coldata by unique run id
  extracted_run_id <- sapply(strsplit(colnames(df.ms1), "_", fixed = TRUE), function(x) tail(x, 1))
  coldata <- colData(obj.filter)
  coldata <- coldata[coldata$run_id %in% extracted_run_id,]
  coldata <- coldata[match(extracted_run_id, coldata$run_id),]
  
  # create summarized exp obj.filter
  se <- SummarizedExperiment(assays = list(MS1.ratio = as.data.frame(df.ms1),
                                           MS1.ratio.log1p = as.data.frame(df.ms1.log1p),
                                           MS2.ratio = as.data.frame(df.ms2),
                                           MS2.ratio.log1p = as.data.frame(df.ms2.log1p),
                                           RT = as.data.frame(df.rt)),
                             rowData = rowdata,
                             colData = coldata)
  
  return(se)
}

file_path=list.files(pattern = "LINCS")[1]
columns=c("pr_gcp_base_peptide", "pr_gcp_histone_mark","A01", "A02", "A03")
load_LINCS_gct <- function(file_paths, columns) {
  data_list <- lapply(file_paths, function(file_path) {
    # Read the file
    data <- read.delim(file_path, skip=2)
    
    subset_columns_by_names <- function(data, search_terms) {
      # Combine all search terms into a single regular expression
      pattern <- paste(search_terms, collapse = "|")
      
      # Find columns that match the pattern
      matching_columns <- sapply(names(data), function(x) grepl(pattern, x))
      
      # Subset the dataframe to keep only matching columns
      subset_data <- data[, matching_columns]
      
      return(subset_data)
    }
    
    # Extract relevant columns and rows
    data <- subset_columns_by_names(data[-c(1:21),], columns)

    # Remove the normalizer peptides
    data <- data[!data$pr_gcp_histone_mark %in% c("H3NORM(41-49)", "H4(68-78)AltNorm"),]
    
    # Split modifications into histone and modification
    data$pr_gcp_histone_mark <- gsub("\\(.*?\\)", "", data$pr_gcp_histone_mark)
    data$pr_gcp_histone_mark <- gsub("H3\\.3", "H33", data$pr_gcp_histone_mark)
    # data$pr_gcp_histone_mark <-  gsub("\\bH3\\b(?!\\.3)", "H31", data$pr_gcp_histone_mark, perl=TRUE)

    res <- split_modifications(trimws(data$pr_gcp_histone_mark))
    data$Protein <- sapply(res, `[`, 1)  

    data$Protein <- gsub("\\bH3\\b", "H31", data$Protein)
    data$Protein[grep("H4", data$Protein)] <- "H4"
    data$Peptide.Sequence <- data$pr_gcp_base_peptide
    data$Peptide.Note <- sapply(res, function(x) {
      paste(x[-1], collapse = "")
    })
    
    # Format modifications
    format_modifications <- function(mods) {
      sapply(mods, function(x) {
        # Remove 1 from 'ac1' and 'me1'
        x <- gsub("ac1", "ac", x)
        x <- gsub("me1", "me", x)
        # Replace 'me0ac0' and 'ac0me0' with 'un'
        x <- gsub("me0ac0", "un", x)
        x <- gsub("ac0me0", "un", x)
        # Replace isolated 'me0' and 'ac0' with 'un'
        x <- gsub("me0", "un", x)
        x <- gsub("ac0", "un", x)
        # Remove any standalone 'me0' or 'ac0' at the end of the string
        x <- gsub("me0$", "", x)
        x <- gsub("ac0$", "", x)
        # Clean up any incorrect 'un' that follows a numerical modification like '1un', '2un', etc.
        gsub("([a-z]+\\d+)un", "\\1", x)
      })
    }
    data$Peptide.Note <- format_modifications(data$Peptide.Note)
    
    # Format H4 modifications
    format_h4_modifications <- function(mods) {
      # Define the base string with all sites as 'unmodified'
      base_string <- "K5unK8unK12unK16un"
      
      # Function to merge modifications into the base template
      sapply(mods, function(x) {
        # Start with the base string
        result <- base_string
        
        # Extract all existing modifications from the input
        # Matches 'K' followed by digits and any characters up to the next 'K' or end of string
        mods_in_string <- regmatches(x, gregexpr("K\\d+[a-z0-9]+", x))
        if(length(mods_in_string[[1]]) > 0) {
          for (mod in mods_in_string[[1]]) {
            # Extract the key (e.g., 'K5', 'K8', etc.)
            key <- substring(mod, 1, regexpr("[a-z]", mod) - 1)
            # Replace the corresponding 'un' entry in the result with the modification found
            result <- sub(paste0(key, "un"), mod, result)
          }
        }
        result
      })
    }
    data$Peptide.Note[data$Peptide.Sequence == "GKGGKGLGKGGAKR"] <- format_h4_modifications(data$Peptide.Note[data$Peptide.Sequence == "GKGGKGLGKGGAKR"])
    
    # Create a new Peptide.Note column
    data <- data %>%
      mutate(Peptide.Note.New = sapply(Peptide.Note, format_peptide_note),
             Peptide.Note.New = paste(Protein, Peptide.Note.New, sep = "-"))
    
    return(data)
  })
  
  return(data_list)
}


calculate_ptm_ratio <- function(obj, assay) {
  
  ratio <- assay(obj, assay) %>% 
    rownames_to_column("Peptide.Note") %>%
    pivot_longer(
      cols = -Peptide.Note,
      names_to = "Replicate.Name",
      values_to = "value") %>%
    merge(as.data.frame(colData(obj)), by.x = "Replicate.Name", by.y = "row.names") %>%
    merge(as.data.frame(rowData(obj)), by.x = "Peptide.Note", by.y = "row.names") %>%
    group_by(Peptide.Sequence, Replicate.Name) %>%
    filter(n_distinct(Peptide.Note) >= 2) %>%
    mutate(ratio = value / sum(value),
           ratio.log1p = log1p(ratio)) %>%
    as.data.frame()
  
  # print("Each peptide ratio sums to 1?")
  # print(all(round(df.verify$sum, digits = 2) == 1))
  
  return(as.data.frame(ratio))
}

calculate_single_ptm_ratio <- function(ptm.ratio, input) {
  
  if(input == "long"){
    df <- ptm.ratio %>%
      select(Peptide.Note, Replicate.Name, ratio) %>%
      spread(key = Replicate.Name, value = ratio)
    row.names(df) <- df$Peptide.Note
    df <- df[,-1]
  } else{
    df <- as.data.frame(ptm.ratio)
  }

  single_ptm_ratio <- df %>%
    rownames_to_column(var = "ptm") %>%
    mutate(protein = sub("-.*", "", ptm)) %>% # separate protein name
    mutate(ptm = gsub(".*?-", "", ptm)) %>% # remove the protein name from the ptm name
    mutate(sptm = strsplit(ptm, ";")) %>% # split the ptm into single ptms
    unnest(sptm) %>% # make each sptm into it's own row
    separate(sptm, into = c("aa", "mod"), sep = "\\[") %>% # separate into aa residue and type of mod
    group_by(protein, aa, mod) %>% # group 
    summarise(across(where(is.numeric), sum), .groups = 'drop') %>% # sum ratios for aa mods across peptides 
    mutate(sptm = paste0(protein, "-", aa, "[", mod)) %>%
    column_to_rownames(var = "sptm") %>%
    select(-c(protein, aa, mod))
  
  return(as.data.frame(single_ptm_ratio))
}

split_modifications <- function(strings) {
  # Use regmatches and regexpr with a regex pattern to capture two groups
  # The pattern captures the amino acid and its number separately from the modification
  sapply(strings, function(x) {
    matches <- regmatches(x, gregexpr("([A-Z][0-9]+)([a-z]+[0-9]*)?", x))
    unlist(matches)
  }, simplify = FALSE)
}

calculate_global_ptm_ratio <- function(sptm.ratio){
  
  me1 <- rownames(sptm.ratio)[grep("\\bme\\b", rownames(sptm.ratio))]
  print("me1 mods:")
  print(me1)
  
  me2 <- rownames(sptm.ratio)[grep("\\bme2\\b", rownames(sptm.ratio))]
  print("me2 mods:")
  print(me2)
  
  me3 <- rownames(sptm.ratio)[grep("\\bme3\\b", rownames(sptm.ratio))]
  print("me3 mods:")
  print(me3)
  
  ac <- rownames(sptm.ratio)[grep("\\bac\\b", rownames(sptm.ratio))]
  print("ac mods:")
  print(ac)
  
  global.ptm.ratio <- t(data.frame(me = colSums(sptm.ratio[rownames(sptm.ratio) %in% me1,], na.rm = TRUE),
                                 me2 = colSums(sptm.ratio[rownames(sptm.ratio) %in% me2,], na.rm = TRUE),
                                 me3 = colSums(sptm.ratio[rownames(sptm.ratio) %in% me3,], na.rm = TRUE),
                                 ac = colSums(sptm.ratio[rownames(sptm.ratio) %in% ac,], na.rm = TRUE)))
  
  return(as.data.frame(global.ptm.ratio))
  
}

calculate_h4_ac_summary_ratio <- function(ptm.ratio, input){
  
  if(input == "long"){
    df <- ptm.ratio %>%
      select(Peptide.Note, Replicate.Name, ratio) %>%
      spread(key = Replicate.Name, value = ratio)
    row.names(df) <- df$Peptide.Note
    df <- df[,-1]
  } else{
    df <- as.data.frame(ptm.ratio)
  }
  
  # Function to count occurrences of "ac"
  count_ac <- function(s) {
    matches <- gregexpr("ac", s)
    sapply(matches, function(x) sum(x > 0))
  }
  
  # Subset strings that start with "H4-"
  h4.peptides <- rownames(df)[grepl("^H4-", rownames(df))]
  
  # Count occurrences of "ac" in each string
  ac.counts <- count_ac(h4.peptides)
  
  
  h4.ac.ratio <- t(data.frame(zero_ac = colSums(df[h4.peptides[ac.counts == 0],], na.rm = TRUE),
                            one_ac = colSums(df[h4.peptides[ac.counts == 1],], na.rm = TRUE),
                            two_ac = colSums(df[h4.peptides[ac.counts == 2],], na.rm = TRUE),
                            three_ac = colSums(df[h4.peptides[ac.counts == 3],], na.rm = TRUE),
                            four_ac = colSums(df[h4.peptides[ac.counts == 4],], na.rm = TRUE)))
  
  return(as.data.frame(h4.ac.ratio))
  
}

calculate_cv <- function(x) {
  sd(abs(as.numeric(x)), na.rm = TRUE) / mean(abs(as.numeric(x)), na.rm = TRUE)
}

na.pad <- function(x,len){
  x[1:len]
}

makePaddedDataFrame <- function(l,...){
  maxlen <- max(sapply(l,length))
  data.frame(lapply(l,na.pad,len=maxlen),...)
}

filter_non_linear_peptides <- function(df.filter, obj, correlation, p_value){
  
  # get median R and p-value across replicates
  quant_filter_summary <- df.filter %>%
    mutate(id = gsub(" ", ";", id)) %>%
    filter(id %in% rownames(obj)) %>%
    group_by(id) %>%
    summarize(correlation = median(correlation),
              p_value = median(p_value)) %>%
    column_to_rownames(var = "id")
  
  # filter out peptides that don't meet criteria
  filter <- quant_filter_summary[quant_filter_summary$correlation <= correlation | 
                                 quant_filter_summary$p_value >= p_value,]
  
  if(nrow(filter) > 0){
    
    print("PTMs filtered due to low R across replicates")
    print(rownames(filter))
    
    obj.filter <- obj[!rownames(obj) %in% rownames(filter)]
    
    write.csv(filter, "filtered_non_quant_peptides.csv",
              row.names = TRUE)
    
    return(obj.filter)
    
  } else {
    print("No PTMs filtered due to low R")
    
    return(obj)
  }
  

}

filter_missing_peptides <- function(obj, threshold = 0.5){
  
  filter <- rownames(obj)[rowSums(is.na(assay(obj, "MS1"))) > sum(obj$multiplier == 1)*threshold]
  
  if(length(filter) > 0){
    
    print("Peptidoforms filtered due to missing values")
    print(filter)
    
    obj.filter <- obj[!rownames(obj) %in% filter,]
    
    write.csv(filter, "filtered_missing_peptides.csv",
              row.names = FALSE)
    
    return(obj.filter)
    
  } else {
    print("No Peptidoforms filtered due to missing values")
    
    return(obj)
  }
}

# IMPROVEMENT: If more than 50% of cells are <= the empty wells
# IMPROVEMENT: Do within each batch
filter_background_peptides <- function(obj, ratio = 1.1){
  
  filter <- assay(obj, "MS1") %>% 
    rownames_to_column("id") %>%
    pivot_longer(
      cols = -id,
      names_to = "variable",
      values_to = "value") %>%
    merge(as.data.frame(colData(obj)), by.x = "variable", by.y = "row.names") %>%
    subset(multiplier %in% c(0, 1)) %>%
    group_by(id, multiplier) %>%
    summarize(median = median(value)) %>%
    summarize(
      median_0cell = median[multiplier == "0"], # Calculate median for multiplier 0
      median_1cell = median[multiplier == "1"], # Calculate median for multiplier 1
      median_ratio = median_1cell / median_0cell, # Calculate ratio of medians (0/1)),
      .groups = 'drop') %>%
    filter(median_ratio < ratio)
  
  if(nrow(filter) > 0){
    
    print("Peptidoforms filtered due to low signal intensity")
    print(filter$id)
    
    obj.filter <- obj[!rownames(obj) %in% filter$id,]
    
    write.csv(filter, "filtered_background_peptides.csv",
              row.names = TRUE)
    
    return(obj.filter)
    
  } else {
    print("No Peptidoforms filtered due to low R")
    
    return(obj)
  }
}

# Need to fix the case when needing to filter both low and high
filter_prop_eff <- function(obj, prop.eff, dev = 2, filter_low = FALSE, filter_high = TRUE) {
  
  names(prop.eff) <- colnames(obj)
  
  # Calculate the overall median and median absolute deviation
  overall_median <- median(prop.eff, na.rm = TRUE)
  overall_mad <- mad(prop.eff, na.rm = TRUE)
  
  # Determine which samples to keep
  if(filter_low){
    samples_to_filter <- prop.eff < (overall_median - dev * overall_mad)
  } else if(filter_high){
    samples_to_filter <- prop.eff > (overall_median + dev * overall_mad)
  } else if(filter_low & filter_low){
    samples_to_keep <- prop.eff < (overall_median - dev * overall_mad) |
      prop.eff > (overall_median + dev * overall_mad)
  }
  
  samples_to_filter <- names(na.omit(samples_to_filter[samples_to_filter == TRUE]))
  
  if(length(samples_to_filter) > 0){
    
    print("Samples filtered due to exceeding MAD threshold")
    print(samples_to_filter)
    
    obj.filter <- obj[,!(colnames(obj) %in% samples_to_filter)]
    
    write.csv(samples_to_filter, "filtered_prop_eff_samples.csv")
    
    return(obj.filter)
    
  } else {
    print("No samples filtered due to exceeding MAD threshold")
    
    return(obj)
  }
  
}

filter_overall_intensity <- function(obj, assay, dev, filter_low = TRUE, filter_high = FALSE) {
  # Extract the assay data and set NA to 0
  df <- assay(obj, assay)
  df[is.na(df)] <- 0
  
  # Calculate the median for each sample
  sample_medians <- colMedians(as.matrix(df))
  
  # Calculate the overall median and median absolute deviation
  overall_median <- median(sample_medians)
  overall_mad <- mad(sample_medians)
  
  # Determine which samples to keep
  if(filter_low){
    samples_to_filter <- sample_medians < (overall_median - dev * overall_mad)
  } else if(filter_high){
    samples_to_filter <- sample_medians > (overall_median + dev * overall_mad)
  } else if(filter_low & filter_low){
    samples_to_keep <- sample_medians < (overall_median - dev * overall_mad) |
      sample_medians > (overall_median + dev * overall_mad)
  }
  
  if(length(samples_to_filter) > 0){
    
    print("Samples filtered due to exceeding MAD threshold")
    print(names(samples_to_filter[samples_to_filter == TRUE]))
    
    obj.filter <- obj[,!samples_to_filter]
    
    write.csv(names(samples_to_filter[samples_to_filter == TRUE]), "filtered_inten_distrib_samples.csv")
    
    return(obj.filter)
    
  } else {
    print("No samples filtered due to exceeding MAD threshold")
    
    return(obj)
  }
  
}

# Bootstrapping function for log2 fold change
bootstrap_log2fc <- function(data, indices, group_column, group1_name, group2_name) {
  sample_data <- data[indices, ]  # Resample with replacement
  g1_values <- sample_data[sample_data[[group_column]] == group1_name, "x"]
  g2_values <- sample_data[sample_data[[group_column]] == group2_name, "x"]
  
  # Exclude NAs from median calculations
  g1_median <- median(g1_values, na.rm = TRUE)
  g2_median <- median(g2_values, na.rm = TRUE)
  
  # Check for zero or NA in median values
  if (g1_median == 0 || g2_median == 0 || is.na(g1_median) || is.na(g2_median)) {
    return(NA)  # Return NA to avoid invalid calculations
  } else {
    return(log2(g2_median / g1_median))
  }
}

differential_abundance_confidence_intervals <- function(summarized_experiment, group_column, group1_name, group2_name, n_bootstrap = 1000) {
  # Extract assay data and metadata
  assay_data <- assay(summarized_experiment)
  col_data <- colData(summarized_experiment)
  
  # Convert group column to character if it's a factor
  if (is.factor(col_data[[group_column]])) {
    col_data[[group_column]] <- as.character(col_data[[group_column]])
  }
  
  # Check if the specified group column and group names exist
  if (!group_column %in% colnames(col_data)) {
    stop("The specified group column does not exist in colData.")
  }
  if (!(group1_name %in% col_data[[group_column]]) || !(group2_name %in% col_data[[group_column]])) {
    stop("The specified group names do not exist in the group column.")
  }
  
  # Perform Mann-Whitney U test and bootstrap for each analyte
  results <- apply(assay_data, 1, function(x) {
    g1 <- x[col_data[[group_column]] == group1_name]
    g2 <- x[col_data[[group_column]] == group2_name]
    p_value <- wilcox.test(g1, g2, exact = FALSE)$p.value  # Use approximation for ties
    
    # Bootstrap for log2 fold change
    bootstrap_results <- boot(data = data.frame(x, col_data[[group_column]]), 
                              statistic = bootstrap_log2fc, 
                              R = n_bootstrap, 
                              group_column = group_column, 
                              group1_name = group1_name, 
                              group2_name = group2_name)
    
    # Calculate confidence intervals, handling cases with NA
    if (any(!is.na(bootstrap_results$t))) {
      log2fc_conf_int <- boot.ci(bootstrap_results, type = "bca")[4:5]
    } else {
      log2fc_conf_int <- c(NA, NA)
    }
    c(log2fc = log2(median(g2) / median(g1)), log2fc_lower = log2fc_conf_int[1], log2fc_upper = log2fc_conf_int[2], p_value = p_value)
  })
  
  # Adjust p-values for multiple testing (Benjamini-Hochberg)
  p_values_adjusted <- p.adjust(as.numeric(sapply(results, function(x) x["p_value"])), method = "BH")
  
  # Format the results into a data frame and include row names
  results_df <- do.call(rbind, results)
  colnames(results_df) <- c("log2FoldChange", "CI_Lower", "CI_Upper", "pValue")
  results_df$pValue <- p_values_adjusted
  rownames(results_df) <- rownames(assay(summarized_experiment))
  
  return(results_df)
}

differential_abundance <- function(summarized_experiment, assay, test_type, group_column, group1, group2) {
  # Extract assay data and metadata
  assay_data <- assay(summarized_experiment, assay)
  col_data <- colData(summarized_experiment)
  
  # set 0 to NA
  assay_data[assay_data == 0] <- NA
  
  # remove NA
  assay_data <- na.omit(assay_data)
  
  # Check if the specified group column and group names exist
  if (!group_column %in% colnames(col_data)) {
    stop("The specified group column does not exist in colData.")
  }
  if (!(group1 %in% unique(col_data[[group_column]])) || !(group2 %in% unique(col_data[[group_column]]))) {
    stop("The specified group names do not exist in the group column.")
  }
  
  # Perform T-test for each analyte
  if(test_type == "t.test"){
    p_values <- apply(assay_data, 1, function(x) {
      g1 <- x[col_data[[group_column]] == group1]
      g2 <- x[col_data[[group_column]] == group2]
      
      # Check for sufficient number of observations
      if (length(g1) < 2 || length(g2) < 2 || all(is.na(g1)) || all(is.na(g2))) {
        return(NA)  # Return NA if there are not enough observations
      } else {
        return(t.test(g1, g2)$p.value)
      }
    })
    # Perform Mann-Whitney U test for each analyte
  } else if(test_type == "wilcox.test"){
    p_values <- apply(assay_data, 1, function(x) {
      g1 <- x[col_data[[group_column]] == group1]
      g2 <- x[col_data[[group_column]] == group2]
      
      # Check for sufficient number of observations
      if (length(g1) < 2 || length(g2) < 2 || all(is.na(g1)) || all(is.na(g2))) {
        return(NA)  # Return NA if there are not enough observations
      } else {
        return(wilcox.test(g1, g2)$p.value)
      }
    })
  }
  
  # adjust p-values
  p_values_adj <- p.adjust(p_values, method = "BH")
  
  # Calculate log2 fold changes (group2 / group1)
  log2_fold_changes <- apply(assay_data, 1, function(x) {
    g1_mean <- median(x[col_data[[group_column]] == group1], na.rm = TRUE)
    g2_mean <- median(x[col_data[[group_column]] == group2], na.rm = TRUE)
    log2(g2_mean / g1_mean)
  })
  
  # Combine results
  results <- data.frame(log2FC = log2_fold_changes, 
                        p.value = p_values,
                        p.adj = p_values_adj,
                        row.names = rownames(assay_data))
  
  results <- results[order(results$p.value, decreasing = FALSE),]
  
  return(results)
}

plot_umap <- function(obj, assay, n_neighbors, color, shape){
  
  set.seed(123)
  res.umap <- umap(t(na.omit(assay(obj, assay))), 
                   n_neighbors = n_neighbors)
  
  res.umap.plot <- res.umap$layout %>%
    merge(colData(obj), by = 'row.names') %>%
    column_to_rownames("Row.names") %>%
    as.data.frame()
  
  # Convert color to factor
  res.umap.plot[[color]] <- as.factor(res.umap.plot[[color]])
  
  # Base ggplot
  p <- ggplot(res.umap.plot, aes(x = V1, y = V2, color = !!sym(color)))
  
  # Add shape aesthetic only if shape is not NULL
  if (!is.null(shape)) {
    res.umap.plot[[shape]] <- as.factor(res.umap.plot[[shape]])
    p <- p + aes(shape = !!sym(shape))
  }
  
  p <- p + 
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
    geom_point(size = 3) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_Publication() +
    theme(legend.position = "right",
          legend.direction = "vertical",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.margin = margin()) +
    guides(color = guide_legend(ncol = 1), 
           shape = guide_legend(ncol = 1))
  
  print(p)
}

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
