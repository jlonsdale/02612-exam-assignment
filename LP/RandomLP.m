function [g, A, b, C, dl, du, l, u] = RandomLP(n)
    % RandomLP generates random linear programming problem data

    m = round(n/2);

    % Create a full rank A matrix
    A = rand(n); 
    A = 10*(A+A') + n*eye(n);
    A = A(:, 1:m);  % Select only the first m columns of A
    b = randn(m, 1); 
    g = randn(n, 1);

    C = rand(n, m);  % Generate a random matrix C with dimensions n x m

    % Generate random upper and lower bounds
    l = -ones(n, 1) * 2.5;
    u = ones(n, 1) * 2.5;
    dl = -2.5 * rand(m, 1);
    du = 2.5 * rand(m, 1);
end

