% 02612 Constrained Optimization
close all
clear
clc

%% Step 1 - prep the inputs
fAS_all = zeros(100,1);
fQPAS_all = zeros(100, 1);
fIP_all = zeros(100,1);
fQPIP_all = zeros(100, 1);
time_IP = zeros(100,1);
time_AS = zeros(100, 1);
time_quadprogIP = zeros(100,1);
time_quadprogAS = zeros(100,1);
broken_constraint_IP = zeros(100,1);
broken_constraint_AS = zeros(100,1);
IP_iter = zeros(100, 1);
quadprogIP_iter = zeros(100,1);
AS_iter = zeros(100,1);
quadprogAS_iter = zeros(100,1);
for i = 5:100
n=i;
[H, g, A, b, C, dl, du, l, u] = randomQPGenerator(i,0.5);
% Convert dl <= C' x <= cu, l <= x <= u to c(x) = Cbar' x + dbar >= 0
Cbar = [C -C eye(n,n) -eye(n,n)];
dbar = [-dl; du; -l; u];
[~,m] = size(Cbar);
[~,j] = size(A);

options = optimoptions('quadprog', 'Algorithm', 'interior-point-convex', 'MaxIterations', 200, 'Display', 'iter', 'OptimalityTolerance', 1e-6);
tic
[x, iter] = InteriorPointQP(H, g, A, b, Cbar, dbar, zeros(n,1));
ticelapsedTime1 = toc;
time_IP(i) = ticelapsedTime1;
IP_iter(i) = iter;
tic
[x2, ~, ~, output] = quadprog(H, g', Cbar',dbar,A',-b, [], [], [], options);
ticelapsedTime2 = toc;
time_quadprogIP(i) = ticelapsedTime2;
quadprogIP_iter(i) = output.iterations;


[x0] = linprog(zeros(1, n)', [A, Cbar]', [b;dbar]);
tic
[x3,lambdaopt,Wset,it] = PrimalActiveSet(H, g, A, -b, Cbar, dbar, x0)  ;
ticelapsedTime3 = toc;
time_AS(i) = ticelapsedTime3;
AS_iter(i) = it;
options = optimoptions('quadprog', 'Algorithm', 'active-set', 'MaxIterations', 200, 'Display', 'iter', 'OptimalityTolerance', 1e-6);
tic;
[x4, ~, ~, output] = quadprog(H, g', Cbar',dbar,A',-b, [], [],x0, options);
ticelapsedTime4 = toc;
time_quadprogAS(i) = ticelapsedTime4;
quadprogAS_iter(i) = output.iterations;


fIP = 0.5*x'*H*x+g'*x;
IP_EQ = (A'*x + b > 1e-6);
IP_EQ = any (IP_EQ == 1);
IP_IQ = (Cbar'*x + dbar > 1e-6);
IP_IQ = any (IP_IQ == 1);

if IP_EQ || IP_IQ
    broken_constraint_IP(i) = 1;
end
fQPIP = 0.5*x2'*H*x2+g'*x2;
fAS = 0.5*x3'*H*x2+g'*x3;
AS_EQ = (A'*x3 + b > 1e-6) ;
AS_EQ = any (AS_EQ == 1);
AS_IQ = (Cbar'*x3 + dbar > 1e-6);
AS_IQ = any (AS_IQ == 1);
if AS_EQ || AS_IQ
    broken_constraint_AS(i) = 1;
end
fQPAS =  0.5*x4'*H*x2+g'*x4;
fQPIP_all(i) = fQPIP;
fIP_all(i) = fIP;
fQPAS_all(i) = fQPAS;
fAS_all(i) = fAS;
end

figure;
plot((abs(fIP_all-fQPIP_all)),"-o");
hold on;
plot((abs(fAS_all-fQPAS_all)),"-o");
xlabel('N');
ylabel('fval');
legend('Interior Point Method', 'Active Set Method');
title('Comparison of absolute difference in function Values');
hold off;


figure;
plot(IP_iter, "-o");
hold on;
plot(quadprogIP_iter, "-o");
plot(AS_iter, "-o");
plot(quadprogAS_iter, "-o");

xlabel('N');
ylabel('Number of Iterations');
legend('Interior Point Method', 'quadprog Active Set', 'Active Set Method', 'quadprog Active Set');
title('Comparison of Number of Iterations');
hold off;



% Ns = 30:10:200;
% avg_time_IP = zeros(length(Ns),1);
% avg_time_QPIP = zeros(length(Ns),1);
% avg_time_AS = zeros(length(Ns),1);
% avg_time_QPAS = zeros(length(Ns),1);
% for k = 1:length(Ns)
%     n=Ns(k);
%     time_quadprogAS = 0;
%     time_quadprogIP = 0;
%     time_AS = 0;
%     time_IP = 0;
% for i = 1:10
% 
% [H, g, A, b, C, dl, du, l, u] = randomQPGenerator(n,0.5);
% % Convert dl <= C' x <= cu, l <= x <= u to c(x) = Cbar' x + dbar >= 0
% Cbar = [C -C eye(n,n) -eye(n,n)];
% dbar = [-dl; du; -l; u];
% [~,m] = size(Cbar);
% [~,j] = size(A);
% 
% options = optimoptions('quadprog', 'Algorithm', 'interior-point-convex', 'MaxIterations', 200, 'Display', 'iter', 'OptimalityTolerance', 1e-6);
% tic
% [x] = InteriorPointQP(H, g, A, b, Cbar, dbar, zeros(n,1));
% ticelapsedTime1 = toc;
% time_IP = time_IP +ticelapsedTime1;
% tic
% [x2] = quadprog(H, g', Cbar',dbar,A',-b, [],[], [], options);
% ticelapsedTime2 = toc;
% time_quadprogIP = time_quadprogIP + ticelapsedTime2;
% 
% 
% [x0] = linprog(zeros(1, n)', [A, Cbar]', [b;dbar]);
% %% Find the working sets, based on the initial guess.
% tic
% [x3,lambdaopt,Wset,it] = PrimalActiveSet(H, g, A, -b, Cbar, dbar, x0)  ;
% ticelapsedTime3 = toc;
% time_AS = time_AS + ticelapsedTime3;
% 
% options = optimoptions('quadprog', 'Algorithm', 'active-set', 'MaxIterations', 200, 'Display', 'iter', 'OptimalityTolerance', 1e-6);
% tic;
% [x4] = quadprog(H, g', Cbar',dbar,A',-b, [], [],x0, options);
% ticelapsedTime4 = toc;
% time_quadprogAS = time_quadprogAS + ticelapsedTime4;
% 
% 
% 
% end
% avg_time_IP(k) = time_IP/10;
% avg_time_QPIP(k) = time_quadprogIP/10;
% avg_time_AS(k) = time_AS/10;
% avg_time_QPAS(k) = time_quadprogAS/10;
% end
% 
% 
% figure;
% 
% plot(Ns, avg_time_IP, '-o');
% hold on;
% plot(Ns, avg_time_QPIP, '-o');
% plot(Ns, avg_time_QPIP, '-o');
% plot(Ns, avg_time_AS, '-o')
% plot(Ns, avg_time_QPAS, '-o')
% xlabel('N');
% ylabel('Time (seconds)');
% legend('Interior Point Method', 'quadprog Active Set', 'Active Set Method', 'quadprog Active Set');
% title('Comparison of Average Runtimes');
% hold off;