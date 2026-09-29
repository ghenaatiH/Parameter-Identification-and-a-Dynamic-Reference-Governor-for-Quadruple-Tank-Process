# Dynamic Matrix Control for Quadruple Tank Systems Integrated with System Identification and Dynamic Reference Governor

This repository contains the scripts, optimization procedures, and numerical results associated with the paper:

**Dynamic Matrix Control for Quadruple Tank Systems Integrated with System Identification and Dynamic Reference Governor**

The repository is organized according to the main computational sections of the paper. The three main folders contain the files related to system identification, DMC control, and DMC integrated with a Dynamic Reference Governor (DRG).

## Repository Structure

### 1. `System Identification`

This folder contains the scripts and numerical results associated with **Section 2.2 (System Identification)** of the paper.

* **`System Identification.m`**
  Performs the parameter identification procedure using two software environments and two optimization methods: **Interior-Point** and **SQP**. The results obtained from the identification procedures and used in the paper are also provided in this folder.

* **`IdentificationResults.m`**
  Displays the identification results presented in **Section 2.2** of the paper.

* **`tau.m`**
  Calculates the time constants discussed in **Section 2.2** of the paper.

* **`QTP.m`**
  Contains the equations describing the height dynamics of the quadruple-tank system.

---

### 2. `DMC_Controler`

This folder contains the scripts and numerical results associated with **Section 5 (DMC Control)** of the paper.

* **`DMC.m`**
  Solves the DMC control problem presented in **Section 5** using three approaches: **Analytical**, **Interior-Point**, and **SQP**.

* **`StepResponc_QTP.m`**
  Calculates the step-response coefficients corresponding to a **12 V input-voltage step** for the quadruple-tank system with **two inputs and two outputs**.

* **`DMC_Cost.m`**
  Defines the objective function of the DMC optimization problem.

* **`NLC.m`**
  Defines the nonlinear constraint function used to model the applied-voltage constraints over the prediction horizon.

* **`results.m`**
  Displays the results and figures corresponding to **Figures 4–9** of the paper.

The numerical results used in the paper for the three DMC solution methods are stored in the following MATLAB files:

* **`Results_DMC_ANA.mat`** — Results obtained using the Analytical method.
* **`Results_DMC_INT.mat`** — Results obtained using the Interior-Point method.
* **`Results_DMC_SQP.mat`** — Results obtained using the SQP method.

---

### 3. `DMC_DRG_Controller`

This folder contains the scripts and numerical results associated with **Section 6 (DMC integrated with Dynamic Reference Governor)** of the paper.

* **`DMC_DRG.m`**
  Solves the DMC-DRG optimization problem presented in **Section 6** using the **Interior-Point** method.

* **`DMC_DRG_Cost.m`**
  Defines the objective function of the DMC-DRG optimization problem.

* **`NLC_Scenario1.m`**
  Defines the nonlinear constraint function for **Scenario 1**, including the applied-voltage constraints and the tank-height constraints over the prediction horizon.

* **`NLC_Scenario2.m`**
  Defines the nonlinear constraint function for **Scenario 2**, including the applied-voltage constraints and the tank-height constraints over the prediction horizon.

* **`Results_Section6.m`**
  Displays the results and figures corresponding to **Figures 11–16** of the paper.

The numerical results used in the paper for the two scenarios are stored in the following MATLAB files:

* **`Results_DMC_DRG_Scenario1.mat`** — Results obtained for Scenario 1.
* **`Results_DMC_DRG_Scenario2.mat`** — Results obtained for Scenario 2.

---

## Reproducibility

The files in this repository contain the computational procedures and numerical results used to generate the results presented in the paper.

The repository is organized so that each main computational section of the paper can be accessed independently:

* **Section 2.2:** System identification
* **Section 5:** DMC control
* **Section 6:** DMC integrated with Dynamic Reference Governor

The `.mat` files contain the numerical results used for the corresponding figures and analyses in the paper, while the MATLAB scripts provide the associated computational procedures.

## Software

The simulations, system identification procedures, and optimization calculations were implemented in **MATLAB** using the optimization methods specified in the corresponding sections of the paper.

## Paper

For the mathematical formulation, system model, system identification procedure, DMC formulation, and DMC-DRG methodology, please refer to the associated paper.
