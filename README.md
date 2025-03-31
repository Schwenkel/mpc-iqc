# 	Output-feedback model predictive control under dynamic uncertainties using integral quadratic constraints 

## Description

This repository contains the MATLAB code that accompanies the paper 
> Lukas Schwenkel, Johannes Köhler, Matthias A. Müller, Frank Allgöwer "Output-feedback model predictive control under dynamic uncertainties using integral quadratic constraints", 2025, arXiv.


## Prerequisites

- [MATLAB](https://de.mathworks.com/products/matlab.html) (tested with version R2023b) 
- [YALMIP](https://yalmip.github.io/) (tested with version 22-June-2023)
- [Mosek](https://www.mosek.com/) (tested with version 10.1.21)
- [CasADi](https://web.casadi.org/) (tested with version v3.5.5)

## Usage

* Run `example1` to get the results of Examples 1 in the article. Specify in line 107 whether to use the initial state optimization and the observer or not.
* The controller K is synthesized using the code in the folder multi-objective-iqc-synthesis that is from [L. Schwenkel, J. Köhler, M.A. Müller, C.W. Scherer, F.Allgöwer. Multi-objective robust controller synthesis with integral quadratic constraints in discrete-time (2025) arxiv:2503.22429](https://arxiv.org/abs/2503.22429)
* For the observer synthesis the function in `obs_synthesis.m` is called.
* For the joint observer and controller analysis the function in `obs_K_analysis.m` is called.
* The MPC scheme is defined and the closed loop is simulated using `mpc_simulation.m` which calls `terminal_conditions.m` to compute the terminal cost and set.

## License

This project is licensed under the MIT License.

## Citation

If you use this code in your research, please cite our work:

```text
@article{Schwenkel2025,
  title={Output-feedback model predictive control under dynamic uncertainties using integral quadratic constraints},
  author={L. Schwenkel and J. K{\"o}hler and M. A. M{\"u}ller and F. Allg{\"o}wer}},
  year={2025}
  journal={arxiv}
}
```
  
## Contact

For any questions or issues related to this code, don't hesitate to get in touch with the author:

- Lukas Schwenkel schwenkel(at)ist.uni-stuttgart.de

## References
Some parts of the code implement methods from the following articles. These parts include a reference in the comments.

- V. Ionescu and M. Weiss, On computing the stabilizing solution of the discrete-time Riccati equation, Linear Algebra and its Applications 174 (1992), 229–238
- C. W. Scherer, Dissipativity and integral quadratic constraints: Tailored computational robustness tests for complex interconnections, IEEE Control Systems 42 (2022), no. 3, 115–139
- C. W. Scherer, J. Veenman, Stability analysis by dynamic dissipation inequalities: On merging frequency-domain techniques with time-domain conditions. Systems & Control Letters 121 (2018), 7–15
- C. W. Scherer, S. Weiland, Linear Matrix Inequalities in Control - Lecture Notes, Delft University of Technology (2005), Delft, Netherlands
- L. Schwenkel, J. Köhler, M. A. Müller, and F. Allgöwer, Robust peak-to-peak gain analysis using integral quadratic constraints, Proc. 22nd IFAC World Congress (2023), 11564–11569
- J. Veenman, C. W. Scherer, and H. Köroglu, Robust stability and performance analysis based on integral quadratic constraints, European J. Control 31 (2016), 1–32

