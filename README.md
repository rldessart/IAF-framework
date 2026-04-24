# IAF-framework
Inclination-based Attribute Framework (IAF) for magnetic data enhancement

# Magnetic Edge Detection Framework: ISA, p-Theta, and IAF

This repository contains GNU Octave / MATLAB-compatible scripts for generating synthetic magnetic models and computing enhancement attributes for magnetic data interpretation.

The implementation includes the Inclination of the Source Analytic Signal (ISA), the p-Theta operator, and the proposed Inclination-Based Angular Fusion (IAF) method.

---

## Description

Magnetic data interpretation often relies on derivative-based attributes to enhance geological structures. However, conventional methods may suffer from instability, noise amplification, or poor edge definition.

This repository provides implementations of the following methods:

- **ISA** – Inclination of the Source Analytic Signal  
- **p-Theta** – Theta-based edge detection operator  
- **IAF** – Inclination-Based Angular Fusion (proposed method)  

Synthetic magnetic models are generated for controlled scenarios, including:

- Two-layer depth model (shallow and deep sources)  
- Noise-contaminated datasets (15% and 25%)  

These models are used to evaluate the trade-off between resolution and stability in magnetic edge detection.

---

## Repository Structure
scripts/
compute_ISA.m
compute_pTheta.m
compute_IAF.m
compute_all_attributes.m

models/
synthetic_models.m

examples/
figure1_script.m
figure2_script.m
figure3_script.m

---


---

## Requirements

- GNU Octave (version 6.0 or later)  
  or  
- MATLAB-compatible environment  

No additional toolboxes are required.

---

## How to Run

1. Open GNU Octave or MATLAB  
2. Navigate to the repository folder  
3. Run:

```octave
figure1_script
figure2_script
figure3_script

These scripts will:

Generate synthetic magnetic data
Compute ISA, p-Theta, and IAF attributes
Reproduce the figures presented in the manuscript

## Reproducibility

All synthetic models, noise simulations, and attribute computations presented in the study can be reproduced using the scripts provided in this repository.

All processing steps, parameters, and noise levels are explicitly defined within the scripts.

## Author

Rafael Lima Dessart

## License

This code is provided for academic and research purposes. Users are free to use, modify, and distribute the scripts, provided that proper citation of the associated publication is given.
