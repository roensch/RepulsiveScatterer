classdef Repulsor < cppclass
    %REPULSOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % q and p are the two parameters of the tangent-point energy
        % q = 6 and p = 12 are quite a good choice for surfcaces.
        q =  6.0;
        p = 12.0;

        theta;

        thread_count;
    end
    
    methods

        function obj = Repulsor( vertex_coords, simplices, theta, threads )
            % vertex_coordinates -  a n x 3 matrix of doubles, ( n = no of vertices)
            % simplices          -  a n x 3 matrix of int64    ( m = no of triangles )
            % The indices have to go from 1 to n.
            % We subtract 1 to take this to indices from 0 to n-1 as used by C++.

            if ( size(vertex_coords,2) ~= 3 )
                error( "First input is expected to be a matrix with 3 columns.");
            end

            if ( size(simplices,2) ~= 3 )
                error( "Second input is expected to be a matrix with 3 columns.");
            end

            if ( min(simplices) < 1 )
                error( "Minimal entry of second input must be greater or equal to 1.");
            end

            if ( max(simplices) > size(vertex_coords,1) )
                error( "Maximal entry in second input must be smaller or equal to column count of first input.");
            end

            obj@cppclass( 'Repulsor_Wrapper', ...
                transpose(vertex_coords), ...
                transpose(simplices-1), ...
                theta, threads ...
            );

            obj.theta = theta;

            obj.thread_count = threads;
        end
        
        function result = TangentPointEnergy(obj)
            result = obj.cppmethod( 'TangentPointEnergy', obj.q, obj.p );
        end

        function result = TangentPointEnergy_Differential(obj)
            % returns derivative in the form of an 3 x n matrix.
            result = obj.cppmethod( 'TangentPointEnergy_Differential', obj.q, obj.p );
        end   

        function result = VertexCount(obj)
            result = obj.cppmethod( 'VertexCount' );
        end   

        function result = VertexCoordinates(obj)
            % returns VertexCount() x 3 matrix of doubles
            result = transpose( obj.cppmethod( 'VertexCoordinates' ) );
        end   

        function result = SimplexCount(obj)
            result = obj.cppmethod( 'SimplexCount' );
        end      

        function result = Simplices(obj)
            % returns SimplexCount() x 3 matrix of int64
            result = double(transpose( obj.cppmethod( 'Simplices' ) + 1 ) );
        end          

        function result = AmbientDimension(obj)
            result = obj.cppmethod( 'AmbientDimension' );
        end          

        function result = DomainDimension(obj)
            result = obj.cppmethod( 'DomainDimension' );
        end  

        function result = Triangulation(obj)
            V = obj.VertexCoordinates();
            S = obj.Simplices();
            result = triangulation(S,V);
        end
        
        function w = TangentPointMetric0(obj, v )
            % Evaluates w = v * A.
            % Expects v to be matrix of size 3 x n and 
            % always returns matrix of size 3 x n!
            if( size(v,1) ~= 3 )
                error("Expecting as input a matrix with 3 rows .")
            end
            if( size(v,2) ~= obj.VertexCount() )
                error("Expecting as input a matrix with VertexCount() columms.")
            end
            w = obj.cppmethod( 'TangentPointMetric0', obj.q, obj.p, v );
        end   

        function w = TangentPointPreconditioner(obj, v)
            % Evaluates w = v * A.
            % Expects v to be matrix of size 3 x n and 
            % always returns matrix of size 3 x n!
            if( size(v,1) ~= 3 )
                error("Expecting as input a matrix with 3 rows .")
            end
            if( size(v,2) ~= obj.VertexCount() )
                error("Expecting as input a matrix with VertexCount() columms.")
            end
            w = obj.cppmethod( 'Preconditioner', obj.q, obj.p, v);
        end

        function result = MaximumSafeStepSize(obj, u, T)
            % u - a 3 x n matrix of a moving direction
            % T > 0 a maximal step size.
            % The function returns T if no self-intersections occur on the
            % path from VertexCoordinates() to VertexCoordinates() + T * u.
            % Otherwise returns (almost) the greates t >= 0 for which this
            % holds.
            % Meant to determine an upper bound for the time step in line
            % search.
            result = obj.cppmethod( 'MaximumSafeStepSize', u, T );
            s = max(max(abs(obj.VertexCoordinates())));
            size = floor(log10(s));
            S = max(max(abs(result * u)));
            Size = floor(log10(S));
            factor = min(1,10^(size-Size));
            result = min([result, factor * result]);
        end

        function SemiStaticUpdate(obj, x)
            % x - a 3 x n matrix of a moving direction
            % Updates VertexCoordinates() to u,
            % but keeps combinatorial data as it was. 
            % Meant to be used for line search.
            obj.cppmethod( 'SemiStaticUpdate', x );
        end

        function T = Remesher(obj, unify_iter, flip_iter, smooth_iter, varargin)
            if nargin == 5
                [V,S] = obj.cppmethod( 'Remesher', obj.thread_count, unify_iter, flip_iter, smooth_iter, varargin{1});
            elseif nargin == 6
                [V,S] = obj.cppmethod( 'Remesher', obj.thread_count, unify_iter, flip_iter, smooth_iter, varargin{1}, varargin{2});
            end

            T = triangulation(  double(transpose(S))+1,transpose(V)    );
            if unify_iter > 0
                fprintf("Number of Vertices: %d \n",size(T.Points,1));
            end
        end

%         function T = Refine(obj, unify_iter, flip_iter, smooth_iter, lower_bound, upper_bound)
%             [V,S] = obj.cppmethod( 'Remesher', obj.thread_count, unify_iter, flip_iter, smooth_iter, lower_bound, upper_bound);
% 
%             T = triangulation(  double(transpose(S))+1,transpose(V)    );
%             if unify_iter > 0
%                 fprintf("Number of Vertices: %d \n",size(T.Points,1));
%             end
%         end

        function [L,m] = cotan_Laplace_mass_matrix(obj)
            C = obj.Simplices();
            V = obj.VertexCoordinates();
            nodes = length(V(:,1));
            L = sparse(nodes,nodes);
            m = sparse(nodes,nodes);
            for i = 1:3
                Cnew = C;
                Cnew(:,i) = [];
                edge1 = (V(Cnew(:,1),:)-V(C(:,i),:)).';
                edge2 = (V(Cnew(:,2),:)-V(C(:,i),:)).';
                cot = (1/2)*(dot(edge1,edge2)./vecnorm(cross(edge1,edge2)));
                L = L + sparse(Cnew(:,1),Cnew(:,2),cot,nodes,nodes);
                A = 0.5*vecnorm(cross(edge1,edge2));
                m = m + sparse(Cnew(:,1),Cnew(:,2),A,nodes,nodes);
            end
            L = L + L.';
            L = L - sparse(1:nodes,1:nodes,sum(L,2),nodes,nodes);
            m = m + m.';
            m = (1/12)*(m+sparse(1:nodes,1:nodes,sum(m,2),nodes,nodes));
        end
    end
end