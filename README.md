# iGB5DOF
[![View iGB5DOF on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/107529-igb5dof)

Drop-in replacement for Bulatov, Reed, Kumar (BRK) function `GB5DOF.m`. We provide a higher fidelity and more versatile, dubbed `iGB5DOF.m` where `i` stands for "improved". Requires `mdl.mat` to be accessible on the path, or you can generate `mdl.mat` yourself on-the-fly. See [`iGB5DOF_example.m`](iGB5DOF_example.m) for details and usage. Note that this is a port of the [`interp5DOF`](https://github.com/sgbaird-5DOF/interp) package which provides a great deal more flexibility in using custom datasets, adjusting parameters, and using various machine learning models. This one is hardcoded to the Ni Olmsted GBs, similar to `GB5DOF.m`.

## Citing
> 1. Baird, S. G.; Homer, E. R.; Fullwood, D. T.; Johnson, O. K. Five Degree-of-Freedom Property Interpolation of Arbitrary Grain Boundaries via Voronoi Fundamental Zone Framework. Computational Materials Science 2021, 200, 110756. https://doi.org/10.1016/j.commatsci.2021.110756.
> 1. Chesser, I., Francis, T., De Graef, M., & Holm, E. A. (2020). Learning the Grain Boundary Manifold: Tools for Visualizing and Fitting Grain Boundary Properties. Acta Materialia. https://doi.org/10.2139/ssrn.3460311
> 1. Francis, T., Chesser, I., Singh, S., Holm, E. A., & De Graef, M. (2019). A geodesic octonion metric for grain boundaries. Acta Materialia, 166, 135–147. https://doi.org/10.1016/j.actamat.2018.12.034
