# Enhanced_inference_for_ITE_quantiles

This repository contains reproducibility materials for the manuscript:

> **Enhanced inference for distributions and quantiles of individual treatment effects in various experiments**  

It includes R code for:
- Simulation studies in Appendix A3
- Real data analyses in Sections 6.2 and 6.3
- Inference for average treatment effects on a binary outcome in Appendix A1

---

## 📂 Repository Structure

```
Enhanced_inference_for_ITE_quantiles/
├── README.md
├── helper_functions.R
├── Simulation Studies/
│     └── sim_section_A3.R
├── Real Data Analysis/
│     ├── education_experiment_section_6.2.R
│     ├── blood_cadmium_section_6.3.R
│     └── education_experiment_section_A1.R
├── output/
│     ├── simulation_figures/
│     └── real_data_figures/
│         ├── education experiment/
│         └── blood cadmium/

```
---

## Data Sources

All datasets used in the analyses are publicly available via R packages:

- Education Experiment Data

```r
data("electric_teachers", package = "RIQITE")
```

- Blood Cadmium Data

```r
data("cadmium", package = "QIoT")
```
---

## Script Execution

Source the helper functions before running other scripts:

```r
source("helper_functions.R")
```

Each script can be run independently after helper functions are loaded:
- **`helper_functions.R`**: Contains utility functions for statistical analysis used across all scripts.
- **`sim_section_A3`**: R code for simulation studies in Appendix A3.
- **`education_experiment_section_6.2`**:
  Real data analysis for evaluating the effectiveness of professional development in Section 6.2.
- **`blood_cadmium_section_6.3`**:
  Real data analysis for evaluating the effect of smoking on blood cadmium levels in Section 6.3.
- **`education_experiment_section_A1`**:
  Code for inference on average treatment effects on a binary outcome in Appendix A1.
---

## Outputs


| Script                      | Manuscript Section | Expected Outputs                | Output Directory                         |
|-----------------------------|--------------------------------------|----------------------------------------|----------------------------------------|
| `sim_section_A3.R`     | Appendix A3| Figure A2, A4       |  `output/simulation_figures/`   |
| `education_experiment_section_6.2.R`     | Section 6.2, Appendix A5.3| Figure 1(a)-(c), Figure 2(a)-(c), Figure A5(a)-(c)      |    `output/real_data_figures/education_experiment/` |
| `blood_cadmium_section_6.3.R`    | Section 6.3   | Figure 3(a)-(c), Figure 4(a)-(c)     |  `output/real_data_figures/blood_cadmium/` |
| `education_experiment_section_A1.R`    | Appendix A1 | Figure A1(a)-(b)       | `output/real_data_figures/education_experiment/`  |


---


