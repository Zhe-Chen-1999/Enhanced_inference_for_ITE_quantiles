# Enhanced_inference_for_ITE_quantiles

This repository contains reproducibility materials for the manuscript:

> **Enhanced inference for distributions and quantiles of individual treatment effects in various experiments**  
> **Manuscript ID:** JASA-T&M-2024-0594.R1

It includes R code for:
- Simulation studies in Section 6.1 and Appendix A3
- Real data analyses in Sections 6.2 and 6.3
- Inference for average treatment effects on a binary outcome in Appendix A1

---

## 📂 Repository Structure

```
Enhanced_inference_for_ITE_quantiles/
├── data/
├── scripts/
│   ├── helpers/
│   │   └── helper_functions.R
│   ├── simulations/
│   │   ├── A1_Inference_for_average_treatment_effects_on_a_binary_outcome.R
│   │   ├── A3_simulation_studies.R
│   │   ├── Sec6.2_Evaluating_the_effectiveness_of_professional_development.R
│   │   └── Sec6.3_Effect_of_smoking_on_the_blood_cadmium_level.R
├── output/
│   └── real_data_figures/
│       ├── education_experiment/
│       └── blood_cadmium/
├── README.md
└── LICENSE
```

---

## ⚙️ How to Set Up

### 1. Clone the Repository

```bash
git clone https://github.com/Zhe-Chen-1999/Enhanced_inference_for_ITE_quantiles.git
cd Enhanced_inference_for_ITE_quantiles
```


### 3. Load Required Data

These analyses use publicly available datasets:

#### Education Experiment Data

From the AER package:

```r
library(AER)
data("electric_teachers", package = "AER")
```

#### Blood Cadmium Data

From the quantreg package:

```r
library(quantreg)
data("cadmium", package = "quantreg")
```

If you store these datasets locally instead, place them in the `data/` folder.

---

### 4. Source Helper Functions

✅ Before running any scripts, always source the helper functions:

```r
source("scripts/helpers/helper_functions.R")
```

---

## ▶️ Script Execution Order

Run the scripts below **in this exact order** for full reproducibility:

---

### ✅ Simulation Studies

1. **Appendix A1 - Binary Outcome Inference**

    **Script:**

    ```
    scripts/simulations/A1_Inference_for_average_treatment_effects_on_a_binary_outcome.R
    ```

    **Outputs:**

    - Figure A1(a)-(b) saved in:
      ```
      output/real_data_figures/education_experiment/
      ```

2. **Section 6.1 & Appendix A3 - Simulation Studies**

    **Script:**

    ```
    scripts/simulations/A3_simulation_studies.R
    ```

---

### ✅ Real Data Analyses

3. **Section 6.2 - Evaluating the Effectiveness of Professional Development**

    **Script:**

    ```
    scripts/simulations/Sec6.2_Evaluating_the_effectiveness_of_professional_development.R
    ```

    **Outputs:**

    - Figure 1(a)-(c)
    - Figure 2(a)-(c)
    - Figure A5(a)-(c)

    Saved in:

    ```
    output/real_data_figures/education_experiment/
    ```

4. **Section 6.3 - Effect of Smoking on Blood Cadmium Levels**

    **Script:**

    ```
    scripts/simulations/Sec6.3_Effect_of_smoking_on_the_blood_cadmium_level.R
    ```

    **Outputs:**

    - Figure 3(a)-(c)
    - Figure 4(a)-(c)

    Saved in:

    ```
    output/real_data_figures/blood_cadmium/
    ```

---

## 📊 Expected Outputs

- **Appendix A1:**
  - Figure A1(a)-(b)

- **Section 6.2:**
  - Figure 1(a)-(c)
  - Figure 2(a)-(c)
  - Figure A5(a)-(c)

- **Section 6.3:**
  - Figure 3(a)-(c)
  - Figure 4(a)-(c)

---


