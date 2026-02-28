# Dual-Task fMRI Brain-Behavior Association Study – Analytic Code Repository

## Overview

This repository contains analytic code and processed behavioral summary measures associated with the cognitive–motor dual-task fMRI study investigating age-related differences in motor–cognitive interference and functional network adaptations.

Neuroimaging data are publicly available via Dataverse:
https://doi.org/10.57782/BALF91

The associated manuscript is available on bioRxiv:
https://doi.org/10.1101/2025.09.28.678786

---

## Repository Contents

### 1. Behavioral Summary Data

Located in:

/behavioral/

- cognitive_behavior_summary.csv  
- motor_behavior_summary.csv  
- Readme.md

These files contain subject-level summary measures derived from task log files.

---

### 2. MATLAB / SPM Analysis Code

Located in:

/analysis_code/ROIs_extraction_from_SPM_activation_results/

Contains scripts used for:
- First-level modeling
- Second-level contrasts
- ROI extraction procedures

---

### 3. R Statistical Analysis (RGCCA)

Located in:

/analysis_code/RGCCA_analysis/

Contains scripts for Brain–behavior correlation analyses
- Cognitive task network output
- Motor task network output
- Behavioral output
- RGCCA modeling and visualization

---

## Software Requirements

MATLAB (R2023)  
SPM12  
R (version 4.4.3)  
Required R packages listed in scripts  

---

## Citation

If you use this repository, please cite:

Deng et al. (2025). Domain-Specific Functional Network Adaptations Supporting Dual-Task Performance in Older Adults. bioRxiv.
https://doi.org/10.1101/2025.09.28.678786

Neuroimaging dataset:
https://doi.org/10.57782/BALF91

---

## Contact

Yan Deng  
Carl von Ossietzky University of Oldenburg  
yan.deng1@uni-oldenburg.de
Supervisor Email: christiane.thiel@uni-oldenburg.de
