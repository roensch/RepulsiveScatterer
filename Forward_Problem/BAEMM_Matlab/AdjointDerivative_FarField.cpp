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

    Int thread_count   = mxGetScalar(input[6]);
    Int gpu_device_num = mxGetScalar(input[7]);
    
    const Int WC = 0;

    cptr<Complex> h = reinterpret_cast<Complex * >( mxGetComplexDoubles(input[5]) );
    output[0]       = mxCreateDoubleMatrix( dim, vertex_count, mxREAL );

    Helmholtz_T H (
        vertex_coordinates,     vertex_count,
        simplices,              simplex_count,
        measurement_directions, meas_count,
        gpu_device_num,
        thread_count
    );

    H.UseDiagonal(true);
    H.SetBlockSize(32);

    const Real cg_tol    = static_cast<Real>(0.00001);
    const Real gmres_tol = static_cast<Real>(0.005);
    
    Complex * du_dn = nullptr;
    
    H.AdjointDerivative_FF<WC>(
        kappa,                   wave_chunk_count,
        incident_directions,     wave_chunk_size,
        h, reinterpret_cast<Real*>( mxGetDoubles( output[0] ) ), 
        du_dn, BAEMM::WaveType::Plane, cg_tol, gmres_tol
    );
}
