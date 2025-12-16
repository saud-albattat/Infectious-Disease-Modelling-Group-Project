# Scabies Outbreak Control in a Closed Population
**Course:** PUBH70060 - Infectious Disease Modelling  
**Date:** 13th December 2024  
**Institution:** Imperial College London

##  Project Overview
This project employs a deterministic compartmental model to simulate the transmission dynamics of *Sarcoptes scabiei* (Scabies) within a closed, isolated island population. The model evaluates the efficacy of different mass drug administration (MDA) strategies to control an outbreak centered around a single school.

##  Authors
* Saud Albahatt
* Carina Tang
* Ekwa Ameyaw
* Xinyi He

##  Scenario & Assumptions
The model simulates a population of **121 individuals** (101 children, 20 staff) on an isolated island.
* **Initial State:** 1 child initially exposed.
* **Transmission:** Occurs via Child-Child, Child-Adult, and Adult-Adult interactions.
* **Intervention:** One-time mass treatment applied to exposed children.

##  Model Dynamics
The simulation utilizes an extended SEIR model with a Treatment (T) compartment:
* **S:** Susceptible
* **E:** Exposed (Latent)
* **I:** Infectious
* **T:** Treated
* **R:** Recovered

##  Code Structure & Contributions

### 1. Core Model
Defines the differential equations governing the flow between compartments ($S_C, E_C, I_U, T, R, S_A, E_A, I_A$) based on contact rates ($\beta$) and latency periods ($\theta$).

### 2. Intervention Analysis (My Contribution)
This module focuses on **Optimization and Comparative Efficacy**. It answers the critical research questions regarding *when* to treat and *which* drug to use.

#### **A. Optimization of Timing and Coverage (Heatmap Analysis)**
* **Objective:** Determine the optimal treatment start day ($t_{start}$) and minimum necessary coverage ($Cov$) to minimize total infections.
* **Methodology:** Ran simulations varying $t_{start}$ (Days 0-100) and Coverage (0-100%).
* **Key Findings:**
    * **Optimal Timing:** The ideal treatment window is around **Day 24**.
    * **Coverage Threshold:** Coverage must be maintained above **30%** to be effective.
    * **Impact:** There is an observed ~25% reduction in infections between 0% and 100% coverage.

#### **B. Comparative Strategy: Cheap vs. Expensive Treatment**
* **Objective:** Compare the effectiveness of two distinct intervention profiles:
    1.  **Permethrin Cream (5%):** Lower efficacy ($\alpha=2$) but achieves **High Coverage** (>50%) due to low cost.
    2.  **Ivermectin (Oral):** High efficacy ($\alpha=4$) but achieves **Low Coverage** (<50%) due to high cost.
* **Methodology:** Plotted the "Peak Infected Children" against the start day of intervention for both drug profiles.
* **Key Findings:** * **Permethrin (High Coverage)** proved to be the more suitable strategy for reducing the peak of infection compared to the high-efficacy/low-coverage alternative.


##  Visualizations
The code generates the following key plots:
* **Total Infections Heatmap:** Visualizes the relationship between intervention timing and coverage.
* **Peak Infection Comparison:** Line graph comparing "Cheap/High Coverage" vs. "Expensive/Low Coverage" scenarios.



##  Dependencies
* R (v4.0+)
* `odin` (for generating the model)
* `ggplot2` (for visualization)
