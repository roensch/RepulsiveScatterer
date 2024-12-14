# RepulsiveScatterer

A Matlab library for iteratively solving the three-dimensional acoustic inverse obstacle scattering problem, where the tangent-point energy is used for regularization.

We offer a pre-included forward solver for the _homogeneous Dirichlet boundary value problem_ and implementations of the boundary-to-far field map, its derivative, adjoint derivative and for the Gauss-Newton step.

# Download/Installation

The library offers a mex-inclusion of the C++-based forward solver _BÄMM_ (https://github.com/HenrikSchumacher/BAEMM) and the tangent-point energy library _Repulsor_ (https://github.com/HenrikSchumacher/Repulsor ).
For installation, run

    git clone --recursive-submodules
    
or clone as usual and run

    git submodule update --init --recursive 
    
To use both, the .cpp-files in "Forward_Problem/BAEMM_Matlab" and "Forward_Problem/Repulsor_Matlab" need to be mex-compiled. Exemplary compilation files for compilation on MacOS are provided. Beforehand, check the requirements for the submodule libraries which can be found in the respective README files (_Repulsor_ and _BÄMM_ require openblas, LAPACK and OpenCL to be installed). Also note that your Matlab version needs to support a compiler which itself at least supports C++20. For further details to get those libraries working we refer to their respective README-files.

Alternatively, you can use your own forward solver, just include it in "Forward_Problem/DirichletTri.m".

# Run examples

For running examples, just start the "main.m" function and specify the experiment (and the name of the file the reconstruction shall be saved in). Details of the experiments (parameters, noise-level, use of other obstacles etc.) can be set in "reconstruction_parameters.m".

The implementation of the library is focused on using multiple wavenumbers for incoming waves and the same set of incident directions for each wavenumber (see "reconstruction_parameters.m" for an exemplary setup). The numbers should always be powes of $2$ and the total number of waves (wavenumber_count x direction_count) should not exceed "64".
The directions of measurement can be customized in "Forward_Problem/DirichletTri.m"

# Credits

Some of the meshes were kindly made freely available by K. Crane (https://www.cs.cmu.edu/~kmcrane/Projects/ModelRepository/) and J. Miller (https://www.thingiverse.com/thing:11622/files).
