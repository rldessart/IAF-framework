# IAF-framework  
**Inclination-Based Attribute Framework (IAF) for Magnetic Data Enhancement**

This repository contains GNU Octave / MATLAB-compatible scripts for generating synthetic magnetic models and computing enhancement attributes used in magnetic data interpretation.

The implementation includes:

- **ISA** – Inclination of the Source Analytic Signal  
- **p-Theta** – Theta-based edge detection operator  
- **IAF** – Inclination-Based Angular Fusion (proposed method)

---

## Description

Magnetic data interpretation commonly relies on derivative-based attributes to enhance geological structures. However, conventional methods may suffer from instability, noise amplification, or poor edge definition.

This repository provides a reproducible framework to evaluate and compare different enhancement techniques, with emphasis on the proposed **IAF method**, which combines angular information to improve structural coherence and stability.

Synthetic magnetic models are used under controlled conditions to isolate the intrinsic behavior of each method.

---

## Synthetic Models

The repository includes synthetic models designed to test performance under different conditions:

- **Two-layer depth model** (shallow and deep sources)  
- **Multi-source configurations**  
- **Noise-contaminated datasets** (e.g., 15% and 25%)

These scenarios allow evaluation of the trade-off between resolution, stability, and noise sensitivity.

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
quick_test.m

---

---

## Requirements

- GNU Octave (version 6.0 or later)  
  **or**  
- MATLAB-compatible environment  

No additional toolboxes are required.

---

## How to Run

1. Open GNU Octave or MATLAB  
2. Navigate to the repository folder  
3. Run:

```matlab
figure1_script
figure2_script
figure3_script

These scripts will:

Generate synthetic magnetic data
Compute ISA, p-Theta, and IAF attributes
Reproduce the figures presented in the study


Reproducibility

All synthetic models, noise simulations, and attribute computations presented in this repository are fully reproducible.

All processing steps, parameters, and noise levels are explicitly defined within the scripts.

## Author

Rafael Lima Dessart

## License

This project is licensed under the MIT License.

Users are free to use, modify, and distribute this code, provided that proper citation of the associated publication is given.
