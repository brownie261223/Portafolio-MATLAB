# MATLAB Numerical Methods Portfolio

Three academic projects from IF321 Numerical Methods in the Physical Engineering program at the National University of Engineering (UNI), Peru. The projects apply numerical methods to quantum mechanics, experimental data fitting, and linear systems.

## Projects

| Project | Method | Reference result |
| --- | --- | --- |
| **Hydrogen Stark Effect** | Assemble a finite Hamiltonian for the `m = 0` sector and compute its eigenvalues with Jacobi rotations. | The two `n = 2` levels are separated by **0.06351 eV** at `2 × 10⁸ V/m` in the `nmax = 4` model. |
| **Luminescence Decay Fit** | Fit an exponential to six coursework observations with Gauss–Newton and examine the residuals. | The fitted decay rate is **0.89425 ms⁻¹**, with an in-sample RMSE of **0.00909 a.u.** |
| **Jacobi and Gauss–Seidel** | Compare iteration counts, spectral radii, and errors against a direct solution. | In the baseline case, Gauss–Seidel met the same step tolerance in **135 iterations**, compared with **247** for Jacobi. |

## Run the examples

1. Download and extract [Portafolio-MATLAB.zip](./Portafolio-MATLAB.zip).
2. In MATLAB, set the current folder to one of the three extracted project folders.
3. Run the following command:

   ```matlab
   run_demo
   ```

Each project includes its source code, a figure, and an English README with the expected output and model assumptions. The reference results were checked with independent calculations. Execution of the packaged demos in MATLAB is pending.
