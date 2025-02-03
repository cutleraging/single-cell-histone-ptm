# single-cell-histone-ptm (sc-hPTM)

Repository for processed data and scripts used in the manuscript: **Mass Spectrometry-based Profiling of Single-cell Histone Post-translational Modifications to Dissect Chromatin Heterogeneity**

## Data Structure

### 1. data - Processed bulk and single-cell data

#### Bulk Data
- **Dele-Oni et al 2021** - DOI: [10.1038/s41597-021-01008-4](https://doi.org/10.1038/s41597-021-01008-4). Global chromatin profiling assay level 3 processed Normalized and QC'ed Numerical Data (QCNORM) retrieved from: [PanoramaWeb](https://panoramaweb.org/project/LINCS/GCP/begin.view?).
- **Sidoli et al 2015** - DOIs: [10.1021/acs.analchem.5b00072](https://doi.org/10.1021/acs.analchem.5b00072), [10.1021/acs.analchem.5b03009](https://doi.org/10.1021/acs.analchem.5b03009). Processed global chromatin profiling assay data taken from supplementary tables.

#### Single-cell Data
- **Sample table** with metadata for injections.
- **Intensities**
  - `*_cells.csv` - Processed global chromatin profiling assay retention time and intensity values (MS1 and MS2 levels) for cell injections.
  - `*_hist_stand.csv` - Processed global chromatin profiling assay retention time and intensity values (MS1 and MS2 levels) for histone standard injections.
  - `*_prop_efficiency.csv` - Processed derivatization peptides retention time and intensity values (MS1 and MS2 levels).

---

### 2. derivatization-efficiency
- Script for calculating derivatization efficiency of single-cells from batches 1-5.

### 3. calibration-curves
#### Cells
- Script for plotting peptidoform titration curves and calculating correlation coefficients of cell injections.
#### Histone-standards
- Script for plotting peptidoform titration curves and calculating correlation coefficients of histone standard injections.

### 4. normalization
#### Auto
- Script for plotting quality control such as intensity distribution and background noise, filtering out low-quality peptidoforms and cells, peptide ratio normalization, batch correction, and data visualization.
  - Used to process cells for benchmarking experiments from batches 1-3.
#### Auto+NABUT
- Same as the `auto` script but used for processing batches 1-5.

### 5. sc-vs-bulk
- Script for comparing pseudobulk single-cell data (batches 1-3) with 100-cell bulk or previous bulk histone studies.

### 6. technical-variability
- Script for calculating variation across histone standard technical replicates and technical replicates from previous bulk studies.

### 7. biological-variability
#### Auto
- Script for calculating cell-to-cell variation at the peptidoform level while accounting for technical variation (measured by histone standards) using single-cell data from batches 1-3.
#### Auto+NABUT
- Script for calculating cell-to-cell variation at the peptidoform level, cell level, and performing hPTM covariance analysis using single-cell data from batches 1-5.

### 8. cell-classification
- Script for training a KNN classifier on labeled single-cell data and evaluating it on a test set.
- Re-trains the KNN classifier on all labeled single-cell data to classify unlabeled single-cell data (batches 1-5).

### 9. differential-abundance
#### With-mixed
- Differential abundance analysis of bulk and single-cell data, including unlabeled "mixed" cells classified using a KNN classifier.
- Results visualized using volcano plots and heatmaps.
- Control and treatment groups further sub-clustered and visualized.
#### Without-mixed
- Same differential abundance analysis as above but excludes unlabeled "mixed" cells.
- Includes correlation of results with the "with-mixed" analysis.

### 10. network-analysis
- Cytoscape analysis file and tables for control and treatment groups containing hPTM-hPTM correlation values and their significance (derived from `7-biological-variability`).

---

### Additional Files
#### cellenONE_runs
- Run files used with the cellenONE instrument for executing the sc-hPTM sample preparation method.

#### functions.R
- Contains all functions used throughout scripts in the analysis pipeline.
