# iGB5DOF
Drop-in replacement for Bulatov, Reed, Kumar (BRK) function `GB5DOF.m`. Higher fidelity and more versatile, dubbed `iGB5DOF.m` where `i` stands for "improved". Requires `mdl.mat` to be accessible on the path, or you can generate `mdl.mat` yourself on-the-fly. See [`iGB5DOF_example.m`](iGB5DOF_example.m) for details and usage. Note that this is a port of the [`interp5DOF`](https://github.com/sgbaird-5DOF/interp) package which provides a great deal more flexibility in using custom datasets, adjusting parameters, and using various machine learning models. This one is hardcoded to the Ni Olmsted GBs, similar to `GB5DOF.m`.
