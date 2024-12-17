#include "mex.h"
#include "matrix.h"

#include <vector>
#include <map>
#include <algorithm>
#include <memory>
#include <string>
#include <sstream>


////////////////////////  BEGIN Step 1: Configuration  ////////////////////////

// Enable profiling

#define LAPACK_DISABLE_NAN_CHECK

// #define REPULSOR_USE_AMD /* Use if suite-sparse is installed on AMD architecture. Then use linker flag ' -lamd'. */

#ifdef __APPLE__
     #define ACCELERATE_NEW_LAPACK /* Use new LAPACK version. */
//      #define TENSORS_USE_ACCELERATE_OLD_LAPACK /* Use old LAPACK version. */

     #include "submodules/Tensors/Accelerate.hpp"
#else
    #include "submodules/Tensors/OpenBLAS.hpp"
#endif

#include "Repulsor.hpp"

using namespace Tools;
using namespace Tensors;
using namespace Repulsor;

using Real       = double;
using Int        = mxInt64;
using LInt       = mxInt64;
using Mesh_T     = SimplicialMesh<2,3,Real,Int,LInt,Real,Real>;

using Mesh_Ptr_T = std::shared_ptr<Mesh_T>;

using Handle_T   = int32_t;
using Map_T      = std::map<Handle_T, Mesh_Ptr_T>;
using Pair_T     = Map_T::value_type;

// List actions
enum class Action
{
    // create/destroy instance - REQUIRED
    New,
    Delete,
    // user-specified class functionality
    VertexCount,
    VertexCoordinates,
    AmbientDimension,
    DomainDimension,
    SimplexCount,
    Simplices,
    TangentPointEnergy,
    TangentPointEnergy_Differential,
    TangentPointMetric0,
    TangentPointMetric,
    Preconditioner,
    SemiStaticUpdate,
    MaximumSafeStepSize,
    Remesher,
    Refine
};

// Map string (first input argument to mexFunction) to an Action
const std::map<std::string, Action> actionTypeMap =
{   
    { "new",                                Action::New },
    { "delete",                             Action::Delete },
    { "VertexCount",                        Action::VertexCount },
    { "VertexCoordinates",                  Action::VertexCoordinates },
    { "AmbientDimension",                   Action::AmbientDimension },    
    { "DomainDimension",                    Action::DomainDimension },    
    { "SimplexCount",                       Action::SimplexCount },
    { "Simplices",                          Action::Simplices },
    { "TangentPointEnergy",                 Action::TangentPointEnergy },
    { "TangentPointEnergy_Differential",    Action::TangentPointEnergy_Differential },
    { "TangentPointMetric0",                Action::TangentPointMetric0 },
    { "TangentPointMetric",                 Action::TangentPointMetric },
    { "Preconditioner",                     Action::Preconditioner },
    { "SemiStaticUpdate" ,                  Action::SemiStaticUpdate },
    { "MaximumSafeStepSize" ,               Action::MaximumSafeStepSize },
    { "Remesher",                           Action::Remesher  },
    { "Refine",                             Action::Refine  }
};

/////////////////////////  END Step 1: Configuration  /////////////////////////


// getHandle pulls the integer handle out of input[1]
inline Handle_T getHandle( int input_count, const mxArray *input[] )
{
    if( input_count < 2 || mxGetNumberOfElements(input[1]) != 1 ) // mxIsScalar in R2015a+
    {
        mexErrMsgTxt("Specify an instance with an integer handle.");
    }
    return static_cast<Handle_T>(mxGetScalar(input[1]));
}

inline Map_T::const_iterator getIterator( const Map_T & m, Handle_T h )
{
    auto it = m.find(h);

    if( it == m.end() ) 
    {
        std::stringstream ss; 
        ss << "No instance corresponding to handle " << h << " found.";
        mexErrMsgTxt(ss.str().c_str());
    }

    return it;
}

// getInstance gets the instances that corresponds to the handle
inline Mesh_Ptr_T getInstance( const Map_T & m, Handle_T h )
{
    return getIterator( m, h )->second;
}


void mexFunction(int output_count, mxArray * output[], int input_count, const mxArray * input[]) 
{

    // static storage duration object for table mapping handles to instances
    static Map_T instances;

    if( input_count < 1 || !mxIsChar(input[0]) )
    {
        mexErrMsgTxt("First input must be an action string ('new', 'delete', or a method name).");
    }

    char * actionCstr = mxArrayToString(input[0]); // convert char16_t to char
    std::string actionStr (actionCstr); 
    mxFree(actionCstr);


    if( actionTypeMap.count(actionStr) == 0 )
    {
        mexErrMsgTxt(("Unrecognized action (not in actionTypeMap): " + actionStr).c_str());
    }

    Action action = actionTypeMap.at(actionStr);

    // If action is not "new" or "delete" try to locate an existing instance based on input handle


    Mesh_Ptr_T M;
    if( action != Action::New && action != Action::Delete ) 
    { 
        M = getInstance(instances, getHandle( input_count, input ) );
    }

	//////// Step 2: customize each action in the switch ////////
    switch( action )
    {
        case Action::New:
        {    
            Profiler::Clear("/Users/Jannik/github/Repulsor_Matlab");
            
            std::pair<Map_T::iterator, bool> result;

            if( input_count != 5 ) 
            {
                mexPrintf("Wrong input. Please input\n (i) a matrix of vertex coordinates,\n (ii) a matrix is simplex indices (zero-based),\n (iii) the separation parameter theta,\n and (iv) the number of threads to use.");
                return;
            }
            else
            {
                cptr<Real> vertex_coordinates = mxGetDoubles(input[1]);
                cptr<Int>  simplices          = mxGetInt64s(input[2]);
             
//                 Real q     = mxGetScalar(input[3]);
//                 Real p     = mxGetScalar(input[4]);
                Real theta = mxGetScalar(input[3]);

                Int thread_count = mxGetScalar(input[4]);
             
//                omp_set_num_threads(thread_count);

                Int amb_dim       = mxGetM(input[1]);
                Int vertex_count  = mxGetN(input[1]);

                Int dom_dim       = mxGetM(input[2])-1;
                Int simplex_count = mxGetN(input[2]);

        
//                 std::cout << "vertex_count  = " << vertex_count  << std::endl;
//                 std::cout << "amb_dim       = " << amb_dim       << std::endl;
//                 std::cout << "simplex_count = " << simplex_count << std::endl;
//                 std::cout << "dom_dim       = " << dom_dim       << std::endl;
//                 std::cout << "thread_count  = " << thread_count  << std::endl;                
//                 std::cout << "theta         = " << theta         << std::endl;

                M = std::make_shared<Mesh_T>(
                    vertex_coordinates,  vertex_count,
                    simplices,          simplex_count,
                    thread_count
                );

                M->cluster_tree_settings.split_threshold                        =  2;
                M->cluster_tree_settings.thread_count                           =  thread_count; // take as many threads as there are used by SimplicialMesh M
                M->block_cluster_tree_settings.far_field_separation_parameter   =  theta;
                M->adaptivity_settings.theta                                    = 10.0;


//                 std::cout << "primitive_count = " << M->GetClusterTree().PrimitiveCount() << std::endl;
                Handle_T H = instances.size() ? (instances.rbegin())->first + 1 : 1;

                // Keep M alive by inserting it into instances.
                result = instances.insert( Pair_T( H, M ) );

            }
    
            if( !result.second ) // sanity check
            {
                mexPrintf("Oh, bad news.  Tried to add an existing handle."); // shouldn't ever happen
            }
            else
            {
                mexLock(); // add to the lock count
            }

		    // return the handle
            output[0] = mxCreateDoubleScalar(result.first->first); // == handle
    
            break;
        }
            
        case Action::Delete:
        {
            instances.erase(
                getIterator(instances, getHandle( input_count, input ) )
            );
            
            mexUnlock();

            output[0] = mxCreateLogicalScalar(instances.empty()); // info

            break;
        }

        case Action::VertexCount:
        {
            mwSize dims[1] = {1};
            output[0] = mxCreateNumericArray(1,&dims[0],mxINT64_CLASS,mxComplexity::mxREAL);
            mxGetInt64s(output[0])[0] = M->VertexCount();
            break;
        }
            
        case Action::VertexCoordinates:
        {
            output[0] = mxCreateDoubleMatrix( M->AmbDim(), M->VertexCount(), mxREAL );
            M->VertexCoordinates().Write( mxGetDoubles(output[0]) );
            break;
        }

        case Action::SimplexCount:
        {
            mwSize dims[1] = {1};
            output[0] = mxCreateNumericArray(1,&dims[0],mxINT64_CLASS,mxComplexity::mxREAL);
            mxGetInt64s(output[0])[0] = M->SimplexCount();
            break;
        }
            
        case Action::Simplices:
            {
            mwSize dims[2] = { static_cast<mwSize>(M->DomDim()+1), static_cast<mwSize>(M->SimplexCount()) };
            output[0] = mxCreateNumericArray(2,&dims[0],mxINT64_CLASS,mxComplexity::mxREAL);
            M->Simplices().Write( mxGetInt64s(output[0]) );
            break;
        }

        case Action::AmbientDimension:
        {
            mwSize dims[1] = {1};
            output[0] = mxCreateNumericArray(1,&dims[0],mxINT64_CLASS,mxComplexity::mxREAL);
            mxGetInt64s(output[0])[0] = M->AmbDim();
            break;
        }

        case Action::DomainDimension:
        {
            mwSize dims[1] = {1};
            output[0] = mxCreateNumericArray(1,&dims[0],mxINT64_CLASS,mxComplexity::mxREAL);
            mxGetInt64s(output[0])[0] = M->DomDim();
            break;
        }

        case Action::TangentPointEnergy:
        {
            const Real q = mxGetScalar(input[2]);
            const Real p = mxGetScalar(input[3]);

            TangentPointEnergy<Mesh_T> energy ( q, p );

            output[0] = mxCreateDoubleScalar( energy.Value(*M) );

            break;
        }

        case Action::TangentPointEnergy_Differential:
        {
            const Real q = mxGetScalar(input[2]);
            const Real p = mxGetScalar(input[3]);

            TangentPointEnergy<Mesh_T> energy ( q, p );

            output[0] = mxCreateDoubleMatrix( M->AmbDim(), M->VertexCount(), mxREAL );

            energy.Differential( *M ).Write( mxGetDoubles(output[0]) );

            break;
        }

        case Action::TangentPointMetric0:
        {
            const Real q = mxGetScalar(input[2]);
            const Real p = mxGetScalar(input[3]);

            TangentPointMetric0<Mesh_T> metric ( q, p );

            output[0] = mxCreateDoubleMatrix( M->AmbDim(), M->VertexCount(), mxREAL );

            // alpha * A * x + beta * y
            metric.MultiplyMetric( 
                *M,
                1.0, mxGetDoubles( input[4]), M->AmbDim(),
                0.0, mxGetDoubles(output[0]), M->AmbDim(),
                M->AmbDim()
            );
            break;
        }

        case Action::Preconditioner:
        {
            const Real q = mxGetScalar(input[2]);
            const Real p = mxGetScalar(input[3]);

            TangentPointMetric0<Mesh_T> metric ( q, p );

            output[0] = mxCreateDoubleMatrix( M->AmbDim(), M->VertexCount(), mxREAL );

            metric.MultiplyPreconditioner( 
                *M,
                1.0, mxGetDoubles( input[4]), M->AmbDim(),
                0.0, mxGetDoubles(output[0]), M->AmbDim(),
                M->AmbDim()
            );
        
            break;
        }
            
        case Action::SemiStaticUpdate:
        {
            M->SemiStaticUpdate( mxGetDoubles(input[2]) );

            break;
        }

        case Action::MaximumSafeStepSize:
        {
            output[0] = mxCreateDoubleScalar( M->MaximumSafeStepSize( mxGetDoubles(input[2]), mxGetDoubles(input[3])[0] ) );

            break;
        }

        case Action::Remesher:
        {
            const Int thread_count = mxGetScalar(input[2]);
            const Int unify_iter = mxGetScalar(input[3]);
            const Int flip_iter = mxGetScalar(input[4]);
            const Int smooth_iter = mxGetScalar(input[5]);
            
            Real lower_bound;
            Real upper_bound;

//             using Remesher_T = SimplicialRemesher<2,3,Real,Int,Real,Int>
//
//             std::shared_ptr< Remesher_T > R = std::make_shared< Remesher_T >(
//                     M->VertexCoordinates().data(),  M->VertexCount(), false,
//                     M->Simplices().data(),          M->SimplexCount(), false,
//                     thread_count
//                 );
//             

            using RemesherBase_T = SimplicialRemesherBase<Real,Int,Real,Int>;
            
            SimplicialRemesher_Factory<RemesherBase_T,2,2,3,3> remesher_factory;
                        
            std::unique_ptr<RemesherBase_T> R = remesher_factory.Make(
                M->VertexCoordinates().data(), M->VertexCount(),  M->AmbDim(), false,
                M->Simplices().data(),         M->SimplexCount(), M->AmbDim(), false,
                nullptr,                                        0,          false,
                thread_count
            );

            if( input_count == 7 ) 
            {
                Real cut = mxGetScalar(input[6]);

                Tensor1<Real,Int> squared_edge_lengths = R->SquaredEdgeLengths();
    
                Sort( squared_edge_lengths.begin(), squared_edge_lengths.end() );
                
                lower_bound = std::sqrt( squared_edge_lengths[Int(Real(R->EdgeCount()-1) * cut)] );
                upper_bound = std::sqrt( squared_edge_lengths[Int(Real(R->EdgeCount()-1) * (1-cut))]  );
            }
            else if( input_count == 8 ) 
            {
                Real cut_lower = mxGetScalar(input[6]);
                Real cut_upper = mxGetScalar(input[7]);

                Tensor1<Real,Int> squared_edge_lengths = R->SquaredEdgeLengths();
    
                Sort( squared_edge_lengths.begin(), squared_edge_lengths.end() );
                
                lower_bound = std::sqrt( squared_edge_lengths[Int(Real(R->EdgeCount()-1) * cut_lower)] );
                upper_bound = std::sqrt( squared_edge_lengths[Int(Real(R->EdgeCount()-1) * (1-cut_upper))]  );
            }
            else
            {
                mexPrintf("Wrong input. ");
                return;
            }

            R->UnifyEdgeLengths( lower_bound, upper_bound, unify_iter );

            R->DelaunayFlip( flip_iter );

            R->TangentialSmoothing( smooth_iter );

            R->SelfCheck();
    
            R->Compress();
            
            R->SelfCheck();
            
            output[0] = mxCreateDoubleMatrix( R->AmbDim(), R->VertexCount(), mxREAL );
            size_t v_size = R->AmbDim() * R->VertexCount();

            mwSize dims[2] = { static_cast<mwSize>(R->DomDim()+1), static_cast<mwSize>(R->SimplexCount()) };
            output[1] = mxCreateNumericArray(2,&dims[0],mxINT64_CLASS,mxComplexity::mxREAL);
            size_t s_size = (R->DomDim()+1) * R->SimplexCount();

            copy_buffer(R->VertexCoordinates().data(), mxGetDoubles(output[0]), v_size, thread_count );
            copy_buffer(R->Simplices().data(), mxGetInt64s(output[1]), s_size, thread_count );

            break;
        }

        default:
        {
            mexErrMsgTxt(("Unhandled action: " + actionStr).c_str());
            break;
        }
    }
    ////////////////////////////////  DONE!  ////////////////////////////////
}
