function x = PCG( A, b, P, TOL, maxIter )
    % A       - a function handle that implements the matrix-vector multiplication
    % b       - a right hand side; can be vector or matrix
    % P       - a function handle that implementes the preconditioner
    % TOL     - tolerance; we try to approximate x = A \ b upto a relative
    %           error TOL in the preconditioner norm
    % maxIter - we abort after this may iterations

    x = zeros(size(b));
    r = b;
    z = P(r);
    p = z;
    rho_0 = dot(r(:),z(:));
    rho = rho_0;
    
    threshold = TOL * TOL * abs(rho_0);

    iter = 0;
    while ( abs(rho) > threshold ) && (iter <= maxIter)
        iter = iter+1;

        u = A(p);
        alpha = rho / dot(p(:),u(:));
        x = x + alpha * p;
        r = r - alpha * u;
        z = P(r);

        rho_old = rho;
        
        rho = dot(r(:),z(:));
        beta = rho / rho_old;
        p = z + beta * p;

    end

    if iter > maxIter
        warning("PCG aborted after requested number of iterations.");
        fprintf('Relative residual achieved     = '); % disp(sqrt(abs(rho/rho_0)));
        fprintf('Number of iterations performed = '); disp(iter);
    else
%         fprintf("PCG succeeded.\n");
%         fprintf('Relative residual achieved    = '); disp(sqrt(abs(rho/rho_0)));
%         fprintf('Number of iterations required = '); disp(iter);
    end
end