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
│   ├── simulation_figures/
│   └── real_data_figures/
│       ├── education experiment/
│       └── blood cadmium/

```
---

## Script Execution

### Source Helper Functions

Source the helper functions before running other scripts:

```r
source("helper_functions.R")
```

### Data Sources

Real data analyses use publicly available datasets:

#### Education Experiment Data

```r
data("electric_teachers", package = "RIQITE")
```

#### Blood Cadmium Data

```r
data("cadmium", package = "quantreg")
```

### Run individual scripts

Each script can be run independently after helper functions are loaded:
- **`helper_functions.R`**: Contains utility functions for statistical analysis used across all scripts.
- **`sim_section_A3`**: R code for simulation studies in Appendix A3.
- **`education_experiment_section_6.2`**:
  Real data analysis for evaluating the effectiveness of professional development in Section 6.2.
- **`blood_cadmium_section_6.3`**:
  Real data analysis for evaluating the effect of smoking on the blood cadmium level in Section 6.3.
   **Outputs:**  
    - Figure 3(a)-(c)  
    - Figure 4(a)-(c)  
  Saved in `output/real_data_figures/blood_cadmium/`

- **`education_experiment_section_A1`**:
  Code for inference on average treatment effects on a binary outcome in Appendix A1.
---


---

### ✅ Simulation Studies

1. **Appendix A1 - Binary Outcome Inference**

2. **Section 6.1 & Appendix A3 - Simulation Studies**

    ```
    scripts/simulations/A3_simulation_studies.R
    ```
### ✅ Real Data Analyses

3. **Section 6.2 - Evaluating the Effectiveness of Professional Development**

4. **Section 6.3 - Effect of Smoking on Blood Cadmium Levels**

   

---

## Expected Outputs


| Script                      | Outputs                              | Output Directory                         |
|-----------------------------|--------------------------------------|----------------------------------------|
| `sim_section_A3.R`     | Figure A2, Figure A4 for Appendix A3      |  `output/simulation_figures/`   |
| `education_experiment_section_6.2.R`     | Figure 1(a)-(c), Figure 2(a)-(c) for Section 6.2, Figure A5(a)-(c) for Section 6.2      |    `output/real_data_figures/education_experiment/` |
| `blood_cadmium_section_6.3.R`    | Figure 3(a)-(c), Figure 4(a)-(c) for Section 6.3      |  `output/real_data_figures/blood_cadmium/` |
| `education_experiment_section_A1.R`    | Figure A1(a)-(b) for Appendix A1      | `output/real_data_figures/education_experiment/`  |


---


