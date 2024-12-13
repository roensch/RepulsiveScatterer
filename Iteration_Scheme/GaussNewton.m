classdef GaussNewton < regmethod
%% iteratively regularized Gauss-Newton type iteration with the tangent-point energy
    properties
        shrink_factor = 0.8;
        alpha_0;
        Armijo_constant;

        max_remesher_discrepancy;

        save_flag;
    end
    
    methods
        function R = GaussNewton(par)
            R = R@regmethod(par);

            R.alpha_0 = par.alpha_0;
            R.Armijo_constant = par.Armijo_constant;
            R.max_remesher_discrepancy = par.max_remesher_discrepancy;

            R.save_flag = par.save_flag;
        end

        function [x_k, stat,F] = solve(R,y_obs,F)
            % x_start: initial guess
            % y_obs  : observed data
            % F      : forward operator (see README.txt)

            alpha = R.alpha_0; % weight parameter for the Tangent-Point term
            tau_max = 10; %f actor for the search_direction
            ArmijoConstant = R.Armijo_constant;
            linesearch_aggression = 1.5;
            R.max_remesher_discrepancy = 2.5 * F.noiselevel * R.stoprule.tau;
            
            % plot initial guess together with true surface
            F.surface.Plot(F.Sdag);
            
            % get initial data
            x_start = F.surface.VertexCoordinates();
            y_k = F.evaluate();
            x_k = x_start;
            
            % check stopping-rule conditions
            [mustExit, R.stoprule] = R.stoprule.stop(F,x_k,y_k,y_obs,R);
            R.condPrint(1,'----------------------------------------------\n');
            [stat,F] = R.error_stat_output(F,[],x_k,y_obs,y_k,x_start);
            R.condPrint(1,'\n');
            
            % calc residual in the far-field space
            residual = y_k - y_obs;
            residual_norm = sqrt(F.dot_Y(residual,residual));

            R.it_step = 0;

            while ~mustExit
                R.it_step = R.it_step+1;
                % calculate the search direction of the gradient descent
                % algorithm

                DE = F.surface.TangentPointEnergy_Differential();

                DT = F.derivative_adjoint(residual) + alpha * DE;

                [search_direction,succeeded] = F.GaussNewtonSolve_TPM0(-DT,alpha);

                % make the Tikhonov functional
                T = R.Tikhonov_functional(F,y_obs,y_k,alpha);
                
                % Armijo line search in search direction
                [F.surface,y_k,tau_max,ArmijoConstant] = F.surface.Armijo_linesearch(T,DT,search_direction,ArmijoConstant,linesearch_aggression*tau_max);

                if residual_norm > R.max_remesher_discrepancy || R.it_step < 2 || succeeded == 0
                    % create remeshed new surface from linsearch result
                    tri_new = F.surface.Remesher(0,5,10,0.05,0.2);
                    F.surface = MeshSurface(tri_new);

                    y_k = F.evaluate();
                else
                    % create new surface from linsearch result
                    
                    tri_new = F.surface.triangulation;
                    F.surface = MeshSurface(tri_new);
                end
                
                % calc new alpha
                alpha = R.shrink_factor * alpha;
                
                x_k = F.surface.triangulation.Points;
                
                % plot updated surface together with true surface
                F.surface.Plot(F.Sdag);
                
                % check stopping-rule conditions
                [stat,F] = R.error_stat_output(F,stat,x_k,y_obs,y_k,x_start);
                condPrint(1,'\n');
                [mustExit, R.stoprule] = R.stoprule.stop(F,x_k,y_k,y_obs,R);

                % calculate new L^2 errors in the far-field
                residual = y_k - y_obs;
                residual_norm = sqrt(F.dot_Y(y_k - y_obs,y_k - y_obs));
            end

            if R.save_flag ~= ""
                filename = ['./Reconstructions/' R.save_flag '.stl'];
                stlwrite(F.surface.triangulation,filename);
            end
        end

        function f = Tikhonov_functional(R,F,y_data,y_eval,a)
            function [t,y] = functional(start)
                if start == true
                    y = y_eval;
                else
                    tic
                    y = F.evaluate();
                end
                t = (1/2) * F.dot_Y(y-y_data,y-y_data) + a * F.surface.TangentPointEnergy();
            end
            f = @functional;
        end
    end
end