## 🔋 Lithium-Ion Battery SoC Estimation and Modeling

This project focuses on **State of Charge (SoC) estimation** for lithium-ion batteries using **data-driven modeling and simulation techniques**.

### 🧠 Key Contributions

- Developed accurate battery models using **MATLAB**, **Simulink**, and **Simscape**, integrating:
  - **Coulomb Counting** method
  - **1-RC equivalent circuit** modeling
- Estimated key battery parameters (**Rs, Rp, Cp**) from **HPPC (Hybrid Pulse Power Characterization)** data using the **Recursive Least Squares with Variable Data Forgetting (RLS-VDF)** algorithm.
- Derived a **nonlinear OCV-SOC relationship** via **polynomial curve fitting** for enhanced modeling fidelity.
- Validated performance through metrics such as:
  - **Root Mean Square Error (RMSE)**
  - **Mean Absolute Error (MAE)**
  - **Mean Absolute Percentage Error (MAPE)**

### 📄 Publication

This work was **peer-reviewed and published at [IEEE INDISCON](https://ieeexplore.ieee.org/)**, highlighting significant contributions to **battery diagnostics and control in modern energy systems**.

### ⚙️ Tools & Technologies

- MATLAB
- Simulink / Simscape
- Curve Fitting Toolbox
- RLS Algorithm Implementation
- HPPC Data Processing

---

> 📌 This repository demonstrates a practical approach to battery SoC estimation for applications in electric vehicles and renewable energy systems.
