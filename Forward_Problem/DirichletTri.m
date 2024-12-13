classdef DirichletTri < BaseOperatorTri
    % The Dirichlet operator for triangulated surfaces
    
    properties
        kappa;
        wave_count;

        device;
    end
    
    methods
        function F = DirichletTri(p)
            F@BaseOperatorTri(p);
            F.op_name = 'DirichletTri';
            F.kappa = p.kappa;
            celldata = load("/Users/Jannik/github/Repulsor_Matlab/Meshes/Sphere_00005120.mat").Expression1;
            data = celldata{1};
            F.meas_dir = data;
            F.wave_count = size(F.kappa,2) * size(F.inc_dir,1);
            F.Ydim = size(F.meas_dir,1) * F.wave_count; % C-dimension of Y
        end
        
        function [data,F] = create_synthetic_data(F,varargin)
           m = F.Sdag;
           
           if nargin == 1
               F.device = int64(0);
           elseif nargin == 2
               % set preferred GPU-device (defaults to F.device=0 if no argument is parsed)
               F.device = int64(varargin{1});
           else
               error('Wrong number of arguments.')
           end
    
               farfield_complex = FarField(transpose(m.Points), transpose(int64(m.ConnectivityList - 1)), F.kappa, transpose(F.inc_dir), transpose(F.meas_dir), int64(F.surface.thread_count), F.device);
               % alternatively insert custom implementation of F to set 'farfield'
    
               farfield(1,:,:) = real(farfield_complex);
               farfield(2,:,:) = imag(farfield_complex);

           [data,F] = addNoise(F,farfield);
        end
        
        function [farfield,F] = evaluate(F)
            m = F.surface;

            farfield_complex = FarField(transpose(m.triangulation.Points), transpose(int64(m.triangulation.ConnectivityList - 1)), F.kappa, transpose(F.inc_dir), transpose(F.meas_dir), int64(m.thread_count), F.device);

             % alternatively insert custom implementation of F to set 'farfield'

            farfield(1,:,:) = real(farfield_complex);
            farfield(2,:,:) = imag(farfield_complex);
        end
        
        function der = derivative(F,h)
            % Takes direction in the form of a 3xn matrix
            m = F.surface;

            if norm(h,"fro") ~= 0
                der_complex = Derivative_FarField(transpose(m.triangulation.Points), transpose(int64(m.triangulation.ConnectivityList - 1)), F.kappa, transpose(F.inc_dir), transpose(F.meas_dir), h, int64(m.thread_count), F.device);
                der(1,:,:) = real(der_complex);
                der(2,:,:) = imag(der_complex);
            else
                meas_count = size(F.meas_dir,1);

                der(1,:,:) = zeros(F.wave_count,meas_count);
                der(2,:,:) = zeros(F.wave_count,meas_count);
            end

            % alternatively insert custom implementation of DF to set 'der'
        end

        function der_adj = derivative_adjoint(F,y)
            % calculate DF*, i.e. the L^2-adjoint of DF, takes direction in the form of a wave_count x meas_count matrix
            m = F.surface;
            
            if norm(y,"fro") ~= 0
                y_c = complex(y(1,:,:) + 1i * y(2,:,:));
                der_adj = AdjointDerivative_FarField(transpose(m.triangulation.Points), transpose(int64(m.triangulation.ConnectivityList - 1)), F.kappa, transpose(F.inc_dir), transpose(F.meas_dir), y_c, int64(m.thread_count), F.device);
            else
                vertex_count = F.surface.VertexCount();
                der_adj = zeros(3,vertex_count);
            end

            % alternatively insert custom implementation of DF* to set 'der_adj'
            
            % In the case of using the tangent point energy, we need DF',
            % i.e. the L^2-waek form of DF*
            if F.metric == "Tangent_Point"
                [L,M] = m.cotan_Laplace_mass_matrix;
                der_adj = der_adj * M;
            end
        end

        function [step,succeeded] = GaussNewtonSolve_TPM0(F,v,alpha)
            % calculate Gauss-Newton type update step.
            m = F.surface;
            
            [step,succeeded] = GaussNewtonSolve(transpose(m.triangulation.Points), transpose(int64(m.triangulation.ConnectivityList - 1)), F.kappa, transpose(F.inc_dir), transpose(F.meas_dir), v, alpha, m.p, m.q, m.theta, int64(m.thread_count), F.device);
            % when using BAEMM, calculate the Gauss-Newton update directly
            % in a MEX-function, otherwise use
%             M = F.Make_Metric(alpha);
%             P = F.Make_Preconditioner(alpha);
% 
%             [step,flag] = F.GaussNewtonSolve(v,M,P);
%             succeeded = ~flag;
        end

        function [step,flag] = GaussNewtonSolve(F,h,M,P)
            % alternative implementation of GN-solver if BAEMM is not used
            % for the evaluation of boundary integral operators
            function y = A(x)
                X = reshape(x,3,[]);
                DF = F.derivative(X);
                Y = F.derivative_adjoint(DF);
                Y = Y + M(X);
                Y = P(Y);
                y = Y(:);
            end

            prec = P(h);

            [Step,flag] = gmres(@A,prec(:),30,1e-2,3);
            step = reshape(Step,3,[]);
        end
        
        function f = Make_Metric(F,alpha)
            m = F.surface.Make_TP_Multiplier();

            function y = M(x)
                y = alpha * m(x);
            end
            f = @M;
        end

        function f = Make_Preconditioner(F,alpha)
            p = F.surface.Make_TP_Preconditioner;

            function y = P(x)
                y = (1/alpha) * p(x);
            end
            f = @P;        
        end
    end
end
