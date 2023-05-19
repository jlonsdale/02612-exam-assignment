function [H, g, A, b] = GenerateRandomEQP(n, density,beta)
    % GenerateRandomEQP generates a random equality-constrained quadratic programming (EQP) problem.
    
    % Inputs: 
    %   n: the size of the problem
    %   density: the density of the sparse matrices
    %   beta: value 0<beta<1 to determine the number of constraints
    
    %   It returns four output arguments: H, g, A, and b, representing the QP data.
    
    % Generate a random number of constraints
    m = round(beta * n);
    
    % Generate a sparse random matrix with given density
    M = sprandn(n, n, density);
    
    % Construct the Hessian matrix with a regularization factor alpha
    H = M' * M + 0.001 * eye(n);
    
    % Generate a sparse random matrix for the equality constraints
    A = sprandn(m, n, density);
    
    % Set the diagonal elements of A to 1
    A(1:m+1:end) = 1;
    
    % Generate random initial solution and Lagrange multipliers
    x_random = randn(n, 1);
    lambda_random = randn(m, 1);
    
    % Construct the KKT matrix
    KKT = [H, -A'; -A, zeros(m, m)];
    
    % use KKT system to obtain a g and b for the random initial solution
    solution = -KKT * [x_random; lambda_random];
    
    % Extract the gradient and Lagrange multipliers from the solution
    g = -solution(1:n);
    b = -solution(n+1:end);
    
end
