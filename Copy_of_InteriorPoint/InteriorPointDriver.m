% 02612 Constrained Optimization
close all
clear
clc

counter=100;
n=10;
results = zeros(1, counter);
timeresults = zeros(1, counter);
options = optimoptions('linprog', 'Algorithm', 'interior-point', 'MaxIterations', 200, 'Display', 'iter');

i = 1;
while i <= counter
    [H, g, A, b, C, dl, du, l, u] = randomQPGenerator(n,0.5);
    Cbar = [C -C eye(n,n) -eye(n,n)];
    dbar = [-dl; du; -l; u];  
    [x_true, fval_true] = linprog( g', Cbar',dbar,A',-b, [], [], options);
    tic;
    [x] = InteriorPointQP(g, A,  -b, Cbar, dbar, zeros(n,1));
    timeresults(i) = toc;
    results(i) = (abs(fval_true)-abs(g'*x)).^2;
    i = i + 1;  % Increment the loop variable
end
abs(fval_true)
abs(g'*x)
format short;

mseValue = sum(results)/counter
std = std(results)

avgtime = mean(timeresults)



