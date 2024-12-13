# RepulsiveScatterer

A Matlab library for iteratively solving the three-dimensional acoustic inverse obstacle scattering problem, where the tangent-point energy is used for regularization.

# Download/Installation
The library offers a mex-inclusion of the forward solver _BAEMM_ (https://github.com/HenrikSchumacher/BAEMM) and the tangent-point energy library _Repulsor_ (https://github.com/HenrikSchumacher/Repulsor ).
For installation, run
    git clone --recursive-submodules
or clone as usual and run
    git submodule update --init --recursive    
To use both, the .cpp-files in Forward_Problem/BAEMM_Matlab and Forward_Problem/Repulsor_Matlab need to be mex-compiled. Exemplary compilation files for compilation on MacOS are prvided. Beforehand, check the requirements for the submodule libraries which can be found in the respective README files (_Repulsor_ and _BAEMM_ require openblas, LAPACK and OpenCL to be installed). Also note that your Matlab version needs to support a compiler which itself at least supports C++20.

Alternatively, you can use your owh forward solver.

# Run examples

For running examples, just start the "main.m" function and specify the experiment (and the name of the file the reconstruction shall be saved in). Details of the experiments (parameters, noise-level, other obstacles etc.) can be set in "reconstruction_parameters.m".
