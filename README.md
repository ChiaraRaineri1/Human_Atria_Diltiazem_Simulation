
# Human Atria Diltiazem Simulation

This repository contains a MATLAB-based simulation project developed to investigate the **bioelectric effects of the drug Diltiazem** on human atrial cardiomyocytes. It was developed by Chiara Raineri and Eleonora Di Consiglio as part of the *Computational Biology of the Heart* course during the academic year 2023-2024.

## 📘 Overview

The study models the effect of **Diltiazem**, a calcium channel blocker, on the **action potential (AP)** of atrial cells. The simulation is based on the **Maleckar et al. (2009)** human atrial cardiomyocyte model, which includes a detailed representation of 30+ ionic current states.

## 🧪 Objectives

- Simulate Diltiazem action as a **blockage of L-type Ca²⁺ channels (ICaL)**
- Assess changes in:
  - Resting potential (Vrest)
  - Maximum potential (Vmax)
  - Difference between Vmax and Vrest (Vdiff)
  - Duration of AP at 90% repolarization (APD90)
- Analyze 30 repeated stimuli, measuring parameters on the last AP
- Evaluate ionic current dynamics (IK1, IKr, IKs, IKur, It)

## 🧠 Theoretical Background

The action potential includes 5 phases from depolarization (phase 0) to full repolarization (phase 4). Diltiazem affects the **plateau phase (phase 2)** by blocking Ca²⁺ channels, thereby shortening AP duration and affecting current dynamics.

## ⚙️ Methods & Implementation

- **Model Used**: Maleckar 2009 model (extension of Nygren et al. 2001)
- **Platform**: MATLAB
- **Files**:
  - `maleckar_main.m`: simulation driver
  - `maleckarnorm.m`: core model definition
- **Diltiazem Blockade Simulation**:
  - `g_CaL` = 6.75 nS × `Const`
  - `Const = [1.0, 0.95, 0.9, ..., 0.05, 0.0]` (0 = full block)

## 📊 Results Summary

### Effects on Action Potential:

- **APD90**:
  - Reduced drastically for full ICaL block (`Const = 0`)
  - Non-linear variation as block intensity increases

- **Vrest**:
  - Generally decreases with Ca²⁺ channel block due to K⁺/Na⁺ redistribution
  - GHK equation used to relate ion permeability with resting voltage

- **Vmax & Vdiff**:
  - Moderate increase for mid-range block
  - Sharp drop in amplitude and peak potential at full block


---

## 📁 Repository Structure

- `CBH_project.m`: Main project file
- `maleckar_main.m`: Simulation driver
- `maleckarnorm.m`: Maleckar cell model file
- Other functions: `find_index.m`

## 👩‍💻 Authors

- **Chiara Raineri** 
- **Eleonora Di Consiglio**

## 📜 License

This project is intended for educational and academic purposes only.  
For licensing terms, see the [LICENSE](LICENSE) file.
