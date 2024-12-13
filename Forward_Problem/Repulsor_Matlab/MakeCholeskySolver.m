function solver = MakeCholeskySolver(S)
% returns a function that acts as preconditioner;
% A is expected to be real, symmetric positive-definite

    % S =
    [U,~,Q] = chol(S);
    P = Q.';
    L = U.';

    function x = f( b, transp )
    % f captures L,U,P,Q from its enclosing scope and keeps a private copy.

        if ~exist('transp','var')
            transp = 0;
        end

        if( transp == 0 )
            x = Q * ( U \ ( L \ ( P * b ) ) ); 
        else
            x = ( ( ( b * Q ) / U ) / L ) * P; 
        end

    end

solver = @f;

end