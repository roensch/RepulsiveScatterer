classdef BaseOperatorTri
    % Badse operator for triangulated surfaces
    
    properties
        op_name = 'Obstacle3dBaseOp';
        noiselevel = 0.01;
        surface;
        Sdag;
        Ydim;
        inc_dir;
        meas_dir;
        metric;
    end
    
    methods
        function F = BaseOperatorTri(p)
            F.Sdag = p.true_surface;
            F.inc_dir = p.inc_dir;
            F.surface = MeshSurface(p.init_guess);
            F.noiselevel = p.noiselevel;
            F.metric = p.metric;
        end

        function [data,F] = addNoise(F,mat)
            noise = randn(size(mat));
            F.noiselevel = F.noiselevel * sqrt(F.dot_Y(mat,mat));
            data = mat + F.noiselevel * noise/sqrt(F.dot_Y(noise,noise));
        end
        
        function v = apply_Gram_X(F,b)
            switch F.metric
                case 'Tangent_Point'
                    m = F.surface.Make_TP_Multiplier();
                    v = m(b);
                otherwise
                    disp('metric not implemented');
            end
        end


        function v = apply_Gram_X_inv(F,b)
            switch F.metric
                case 'Tangent_Point'
                    v = F.surface.TangentPointMetric_inverse(b,1e-4);
                otherwise
                    disp('metric not implemented');
            end
        end

        function s = dot_Y(F,A,B)
            s = sum(A.*F.apply_Gram_Y(B),'all');
        end

        function v = apply_Gram_Y(F,b)
            v = (4*pi/F.Ydim) * b;
        end

        function v = apply_Gram_Y_inv(F,b)
            v = (F.Ydim/(4*pi)) * b;
        end
    end
end

