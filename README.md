# PINNs_Thesis
Code in Pytorch for Physics-Informed Neural Networks (PINNs) for my diploma thesis

# Physics-Informed Neural Networks (PINNs) Thesis Code

This repository contains the code developed for my Master's thesis on Physics-Informed Neural Networks (PINNs) for solving wave equations. The thesis itself, as well as the comments in the code, are written in **Greek**.  

## Overview

The repository includes scripts for simulating various wave equations using PINNs, including:

- Linear Wave Equation
- Nonlinear Korteweg–de Vries (KdV) Equation
- Non-Linear Schrödinger (NLS) Equation
- Burgers' Equation
- Camassa-Holm Equation
- Matlab codes for

## Additional Material

The repository also includes MATLAB codes used for the numerical comparison of the results. In particular, a finite difference scheme combined with a 4th-order Runge–Kutta (RK4) method is implemented to validate the PINN approximations.

## Numerical Comparison Workflow

For the numerical comparison between PINNs and classical methods, the MATLAB implementation must be executed first.

The MATLAB scripts generate data files, which are then imported into the Python scripts for direct comparison with the PINN results.

Therefore, in order to reproduce the comparison results:

1. Run the corresponding MATLAB script.
2. Ensure that the generated data file is saved in the appropriate directory.
3. Execute the Python comparison script.

The Python code expects the MATLAB-generated file to be present in the specified path.


Each script is **well-commented in Greek** to explain the methodology and implementation details. 

## Requirements

- Python 3.x
- PyTorch
- NumPy, Matplotlib
- CUDA (optional, for GPU acceleration)

## Notes

- The code follows the methodology described in the thesis.
- The comments are in Greek to match the thesis documentation.
- Feel free to explore and adapt the code for further experiments.

## Contact

For questions or clarifications, you can contact me at [spyridonsakellaropoulos@gmail.com].
