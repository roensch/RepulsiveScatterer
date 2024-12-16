#include "mex.h"
#include "matrix.h"

#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <pwd.h>
#include <filesystem>

#define LAPACK_DISABLE_NAN_CHECK

#define ACCELERATE_NEW_LAPACK /* Use new LAPACK version. */
// #define TENSORS_USE_ACCELERATE_OLD_LAPACK /* Use old LAPACK version. */

// #define REPULSOR_USE_AMD /* Use if suite-sparse is installed on AMD architecture. Then use linker flag ' -lamd'. */

#include "Helmholtz_OpenCL.hpp"

using namespace Tools;
using namespace Tensors;
using namespace Repulsor;

using Helmholtz_T = BAEMM::Helmholtz_OpenCL;

using LInt = mxInt64;
using Int = mxInt64;
using Real = double;
using Complex = std::complex<Real>;

void mexFunction(int output_count, mxArray * output[], int input_count, const mxArray * input[]) 
{
    cptr<Real> vertex_coordinates        = mxGetDoubles(input[0]);
    cptr<Int>  simplices                 = mxGetInt64s(input[1]);

    Int dim           = 3;
    Int vertex_count  = mxGetN(input[0]);
    Int simplex_count = mxGetN(input[1]);
    
    cptr<Real> kappa                     = mxGetDoubles(input[2]);
    cptr<Real> incident_directions       = mxGetDoubles(input[3]);
    cptr<Real> measurement_directions    = mxGetDoubles(input[4]);
    
    Int meas_count       = mxGetN(input[4]);
    Int wave_chunk_count = mxGetN(input[2]);
    Int wave_chunk_size  = mxGetN(input[3]);
    Int wave_count       = wave_chunk_count * wave_chunk_size;
    
    Real regpar = mxGetScalar(input[6]);

    const Real p     =  mxGetScalar(input[7]);
    const Real q     =  mxGetScalar(input[8]);
    const Real theta =  mxGetScalar(input[9]);
    
    Int thread_count   = mxGetScalar(input[10]);
    Int gpu_device_num = mxGetScalar(input[11]);
    
    const Int WC = 0;
    
    cptr<Real> h     = mxGetDoubles(input[5]);
    output[0]        = mxCreateDoubleMatrix( dim, vertex_count, mxREAL );
    
    Helmholtz_T H (
        vertex_coordinates,     vertex_count,
        simplices,              simplex_count,
        measurement_directions, meas_count,
        gpu_device_num,
        thread_count
    );
    
    H.UseDiagonal(true);
    H.SetBlockSize(32);

    Real cg_tol          = static_cast<Real>(0.00001);
    Real gmres_tol_inner = static_cast<Real>(0.005);
    Real gmres_tol_outer = static_cast<Real>(0.01);
    
    Complex * du_dn = nullptr;
    
    // set up the TP-metric
    
    using Mesh_T     = SimplicialMesh<2,3,Real,Int,LInt,Real,Real>;
    using Mesh_Ptr_T = std::shared_ptr<Mesh_T>;

    Mesh_Ptr_T M = std::make_shared<Mesh_T>(
        vertex_coordinates,    vertex_count,
        simplices, simplex_count,
        thread_count
        );

    M->cluster_tree_settings.split_threshold                        = 2;
    M->cluster_tree_settings.thread_count                           = thread_count; // take as many threads as there are used by SimplicialMesh M
    M->block_cluster_tree_settings.far_field_separation_parameter   = static_cast<Real>(theta);
    M->adaptivity_settings.theta                                    = static_cast<Real>(10.0);

    const Real s = (p - 2) / q;

    TangentPointMetric0<Mesh_T> tpm (q,p);
    
    auto A = [&]( cptr<Real> X, mptr<Real> Y )
    {
        tpm.MultiplyMetric( *M,
            regpar,            X, dim,
            Scalar::Zero<Real>, Y, dim,
            dim
        );
    };

    Real one_over_regpar = Inv<Real>(regpar);

    auto P = [&]( cptr<Real> X, mptr<Real> Y )
    {
        tpm.MultiplyPreconditioner( *M,
            one_over_regpar,    X, dim,
            Scalar::Zero<Real>, Y, dim,
            dim
        );
    };
    
    // perform the Gauss-Newton solve
    
    Int succeeded;
    
    succeeded = H.GaussNewtonSolve<WC>(
        kappa,                   wave_chunk_count,
        incident_directions,     wave_chunk_size,
        A,  P,
        Scalar::One<Real>,  h,                          dim,
        Scalar::Zero<Real>,  mxGetDoubles(output[0]),    dim, 
        du_dn, BAEMM::WaveType::Plane, 
        cg_tol, gmres_tol_inner, gmres_tol_outer
    );

    output[1] = mxCreateLogicalScalar(succeeded);
    
    if( succeeded != 1 )
    {
        mexErrMsgTxt("GMRES reached maximal iterations.");
    }
}
