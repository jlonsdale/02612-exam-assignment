function [g, A, b, C, dl, du, l, u] = RandomLP(n)   
    m=round(n);
    g = randn(n, 1); % Generate a random vector g  
    % Generate random matrices A and C
    % Create a full rank A matrix
    A = rand(n); 
    A = 10*(A+A')+n*eye(n);
    A = A(:,1:m); 
    b = randn(m,1); 
  
    % Generate random upper and lower bounds
    l = -ones(n,1)*2.5;
    u = ones(n,1)*2.5;

    C = rand(n,m);
    % Generate random upper and lower bounds
    dl = -2.5*rand(m,1);
    du = 2.5*rand(m,1);
end
