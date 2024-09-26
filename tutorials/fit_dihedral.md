# Fitting proper torsional constants for OPLS-AA force field

This tutorial is dedicated to run the `fit_dihedral.py` script. Before following the tutorial, be sure that scripts from the [DICEtools](https://github.com/hmcezar/dicetools) library are in your `$PATH` or that `fit_dihedral.py` is inside the DICEtools directory.

## Ryckaert-Bellemans dihedral

As developed by [William L. Jorgensen et al.](https://doi.org/10.1021%2Fja00214a001), the OPLS-AA force field describes the proper torsional interaction according to the Ryckaert-Bellemans dihedral. For each proper torsional angle, the 4-atom energy is given by,

$$E_{RB}=\sum_{n=1}^6C_n[\cos{(\phi+\delta_n)}]^n\equiv\frac{1}{2}\sum_{n=1}^4F_n[1+(-1)^{n+1}\cos{(n\phi+\delta_n)}] \ .$$

Both forms are equivalent, but the additional flexibility from the 6 coefficients ($C_n$) had a better performance in comparison to the 4 coefficients ($F_n$) in the Fourier form. Therefore, the `fit_dihedral.py` considers the full expression with 6 coefficients and all $\delta_n$ phases are set to zero as frequently done.

## Initial setup with DICEtools

