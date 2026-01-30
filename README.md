# Adaptive Input Shaper Design for Unknown Second-Order Systems with Real-Time Parameter Estimation

<p align="center">
  <b>Feedforward control scheme with parameter estimation and optimal input shaping for vibration suppression</b><br>
  <img src="Figures/ff_block.png" width="80%">
</p>

This work propose a feedforward control method that not only estimates system parameters online for a black-box second- order system, but also designs the input shaper adaptively in real-time. This repostitory contains the code, data, and methods used to design optimal input shaper and parameter estimator. 

**Note:** This work has been accepted for publication at **ACC 2026**.  
A preprint version is available on arXiv: https://arxiv.org/abs/2601.17210

---

## Project Structure

  - The `Libraries` folder contains custom modules for system modeling, parameter estimation, and optimal input shaper design used by the main scripts.
  - `run_ATDF.py` – Runs the feedforward simulation shown in the figure above, where the true system parameters $\zeta$ and $\omega_\mathrm{n}$ are estimated on the fly and the optimal input shaper is designed accordingly. Multiple simulations are supported with different sets of system parameters, allowing the algorithm to evaluate a diverse range of second-order system scenarios.

  - `ATDF_stepwise.py` – Simulates an adaptive input shaper with multiple step-wise changes after the estimation period, reflecting real-world scenarios such as 3D-printing motion and gantry-crane movements.

  - `parameter_sweep.py` – Generates 3D surface plots showing the dependency of $\mathcal{A}$ and $\mathcal{T}$ over ranges of $\omega_{\mathrm{n}}$ and $\zeta$.

---

## Requirements

- matplotlib==3.10.3
- numpy==2.3.0
- pandas==2.3.0
- pyDOE==0.3.8
- scipy==1.15.3
- seaborn==0.13.2

Install dependencies using:

```bash
pip install -r requirements.txt
```
---

## How to Run

### Step 1: Setup

Download or clone the repository:

```bash
git clone https://github.com/NyiNyi-14/A-TDF.git
```

Make sure all scripts are in the same directory.

 ### Step 2: Parameter Customization

Before running the code, adjust the system parameters to configure your simulation:

- **System parameters**: `omega_test`, `zeta_test`  
- **Time values**: `duration`, `identification_duration`, `dt`  

### Step 3: Simulate the System

Run `run_ATDF.py` to observe the adaptive input shaper’s single-step behavior over a range of estimated $\omega_{\mathrm{n}}$ and $\zeta$ values:

```bash
python run_ATDF.py
```

Run `ATDF_stepwise.py` to demonstrate the optimal input shaper performance under multiple step-wise reference changes:
```bash
python ATDF_stepwise.py
```

---

## Results

<p align="center">
  <b>Performance of the proposed method:</b><br>
  <img src="Figures/run_adaptive_1.png" width="100%"><br>
  (a–d) Varying &zeta; with &tau; = 2 s and &omega;<sub>n</sub> = &pi; rad/s <br>
  (e–h) Varying &omega;<sub>n</sub> with &tau; = 2 s and &zeta; = 0.707 <br>
  (i–l) Varying &tau; with &omega;<sub>n</sub> = 3&pi; rad/s and &zeta; = 0.707
</p>

<p align="center">
  <b>Feedforward control for step-wise reference tracking:</b><br>
  <img src="Figures/StepWise.png" width="70%"><br>
  (Top) &omega;<sub>n</sub> = 3&pi; rad/s  (b) &omega;<sub>n</sub> = 30&pi; rad/s
</p>

<p align="center">
  <b>Dependency of input shaper parameters 𝒜 and 𝒯 on &zeta; and &omega;<sub>n</sub> under parameter sweep:</b><br>
  <img src="Figures/3D_plot.png" width="70%"><br>
</p>


---

## Related Work

This project builds on developed control mechanisms, including:

- Design optimal input shaper online
- Real time parameter estimation

---

## Citation

If you use this work, please cite the related paper as follows:

---
<!-- 
## Author

**Nyi Nyi Aung** 

PhD Student, Mechanical and Industrial Engineering - LSU, USA

MSc, Sustainable Transportation and Electrical Power Systems - UniOvi, Spain

BE, Electrical Power - YTU, Myanmar

##

**Bradley Wight** 

B.S. Student, Mechanical Engineering

Louisiana State University

##

**Adrian Stein, PhD**

Assistant Professor

Department of Mechanical and Industrial Engineering

Louisiana State University

## -->
