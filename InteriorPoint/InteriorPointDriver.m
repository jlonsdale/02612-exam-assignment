% 02612 Constrained Optimization
close all
clear
clc

counter=50;
n=100;
results = zeros(1, counter);
timeresults = zeros(1, counter);
options = optimoptions('quadprog', 'Algorithm', 'interior-point-convex', 'MaxIterations', 200, 'Display', 'iter', 'OptimalityTolerance', 1e-6);

i = 1;
while i <= counter
    [H, g, A, b, C, dl, du, l, u] = randomQPGenerator(n,0.75);
    Cbar = [C -C eye(n,n) -eye(n,n)];
    dbar = [-dl; du; -l; u];  
    [x_true, fval_true] = quadprog(H, g', Cbar',dbar,A',-b, [], [], [], options);
    tic;
    [x] = InteriorPointQP(H, g, A, -b, Cbar, dbar, zeros(n,1));
    timeresults(i) = toc;
    results(i) = (fval_true-0.5*x'*H*x-g'*x).^2;
    i = i + 1;  % Increment the loop variable
end

format short;

mseValue = mean(results)
avgtime = mean(timeresults)



