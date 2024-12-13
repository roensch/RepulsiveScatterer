classdef MeshSurface < Repulsor
    % This is the class of general triangulated surfaces with several
    % implemented choices of energy

    properties
        prec_reg_fact = 0.001;

        triangulation;
    end

    methods
        function surface = MeshSurface(triang)
            theta = 1/8;
            threads = 4;
            surface@Repulsor(triang.Points, int64(triang.ConnectivityList), theta, threads);
            surface.triangulation = triang;
        end

        function grad = TangentPointMetric_inverse(surface,u,TOL)
                A = surface.Make_TP_Multiplier();
                P = surface.Make_TP_Preconditioner();

                grad = PCG(A,u,P,TOL,500);
        end
        
        function [surface,y,tau,sigma,varargout] = Armijo_linesearch(surface,F,DFx,search_direction,Armijo_constant,alpha_max)
            [E,~] = F(true);
            
            % calc derivative of Tikhonov functional in search_direction
            DE = sum(dot(DFx,search_direction),2);
            sigma = Armijo_constant;
            s = 1/(4*(1-sigma));
            
            % determine max step size to prevent collision
            tau0 = surface.MaximumSafeStepSize(search_direction,alpha_max);
            tau = tau0;          

            k = 0;
            maxIter = 10;
            
            coord = surface.triangulation.Points.';
            
            % update surface coordinates
            surface.Update(coord + tau*search_direction);
            
            [Enew0,y] = F(false);
            Enew = Enew0;

            done = false;
            inner_it = 0;

            % start linesearch iterations
            while ~done
                while (Enew > (E + tau * sigma * DE)) && (k <= maxIter) || isnan(Enew)
                    % update step-size
                    tau_interpolation = -(1/2) * tau^2 * DE / (Enew - E - tau * DE);
                    tau = max([s*tau,tau_interpolation]);
                    
                    % update surface coordinates
                    surface.Update(coord + tau * search_direction);

                    [Enew,y] = F(false);
                    k = k+1;
                end
                if k > maxIter && inner_it < 2
                    disp('Armijo constant gets shrunk.')
                    
                    % update surface coordinates
                    surface.Update(coord + tau*search_direction);

                    tau = tau0;
                    Enew = Enew0;

                    % update Armijo constant and shrink-parameter s
                    sigma = (1/2) * sigma;
                    s = 1/(4*(1-sigma));
                    k = 0;
                    inner_it = inner_it + 1;
                else
                    done = true;

                    % successively adapt Armijo constant, depending on
                    % linesearch performance
                    switch k
                        case {0}
                            sigma = min([2 * sigma,1/2]);
                        case {10,11}
                            sigma = (1/2) * sigma;
                        otherwise
                    end
                    fprintf('Linesearch iterations: %i\n',k);
                end
            end
            
            if nargout == 5
                varargout{1} = Enew;
            elseif nargout > 5
                disp('Wrong numbero of output arguments.');
            end
        end
        
        function f = Make_TP_Multiplier(surface)
            f = @surface.TangentPointMetric0;
        end

        function prec = Make_TP_Preconditioner(surface)            
            prec = @surface.TangentPointPreconditioner;
        end

        function surface = Update(surface, u)
            % update triangulation during linesearch. Do not use this to
            % update after complete GN-iteration (need to create new
            % MeshSurface in this case, see GaussNewton.m
            % u is supposed to be a 3 x n matrix

            Conn = surface.triangulation.ConnectivityList;
            surface.triangulation = triangulation(Conn, u.');
            surface.SemiStaticUpdate(u);
        end
        
        function u = ApplyMass(surface,v)
            [laplacian,mass] = surface.cotan_Laplace_mass_matrix();
            u = v * mass;
        end

        function u = ApplyMassInverse(surface,v)
            [laplacian,mass] = surface.cotan_Laplace_mass_matrix();
            sol = MakeCholeskySolver( mass );
            u = sol(v,1);
        end

        function Plot(surface,varargin)
            points = surface.triangulation.Points;
            h = trisurf(surface.triangulation.ConnectivityList, points(:,1), points(:,2), points(:,3),'LineStyle','none');
            if nargin == 2
                hold on
                trisurf(varargin{1},'FaceAlpha', 0,'EdgeColor', 'c','LineWidth',0.05,'LineStyle',':')
                hold off
            end
            drawnow
        end

  
        function FieldPlot(surface,u,scale)
            points = surface.triangulation.Points;
            quiver3( points(:,1), points(:,2), points(:,3), u(:,1),u(:,2),u(:,3), scale, 'r' )
        end
    end
end

