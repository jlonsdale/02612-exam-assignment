n_values = 5:10:200;

num_values = length(n_values);

diff_ldl_dense = zeros(1, num_values);
diff_ldl_sparse = zeros(1, num_values);
diff_lu_dense = zeros(1, num_values);
diff_lu_sparse = zeros(1, num_values);
diff_null_space = zeros(1, num_values);
diff_range_space = zeros(1, num_values);

time_quadprog = zeros(1, num_values);
time_ldl_dense = zeros(1, num_values);
time_ldl_sparse = zeros(1, num_values);
time_lu_dense = zeros(1, num_values);
time_lu_sparse = zeros(1, num_values);
time_null_space = zeros(1, num_values);
time_range_space = zeros(1, num_values);


for i = 1:length(n_values)
    n = n_values(i);
    [H, g, A, b] = GenerateRandomEQP(n, 0.15);
    
    tic;
    [x_true, fval_true] = quadprog(H, g', [], [], A, b);
    time_quadprog(i) = toc;
    
    tic;
    [x1] = EQPsolver(H, g, A, b, "LDLdense");
    fval1 = calculateFval(x1, H, g);
    time_ldl_dense(i) = toc;
    diff_ldl_dense(i) = norm(fval1 - fval_true);

    tic;
    [x2] = EQPsolver(H, g, A, b, "LDLsparse");
    fval2 = calculateFval(x2, H, g);
    time_ldl_sparse(i) = toc;
    diff_ldl_sparse(i) = norm(fval2 - fval_true);

    tic;
    [x3] = EQPsolver(H, g, A, b, "LUdense");
    fval3 = calculateFval(x3, H, g);
    time_lu_dense(i) = toc;
    diff_lu_dense(i) = norm(fval3 - fval_true);
    
    tic;
    [x4] = EQPsolver(H, g, A, b, 'LUsparse');
    fval4 = calculateFval(x4, H, g);
    time_lu_sparse(i) = toc;
    diff_lu_sparse(i) = norm(fval4 - fval_true);

    tic;
    [x5] = EQPsolver(H, g, A, b, "NullSpace");
    fval5 = calculateFval(x5, H, g);
    time_null_space(i) = toc;
    diff_null_space(i) = norm(fval5 - fval_true);

    tic;
    [x6] = EQPsolver(H, g, A, b, "RangeSpace");
    fval6 = calculateFval(x6, H, g);
    time_range_space(i) = toc;
    diff_range_space(i) = norm(fval6 - fval_true);
end

% Plot the differences for each method
figure;
plot(n_values, diff_ldl_dense, '-o', 'LineWidth', 2);
hold on;
plot(n_values, diff_ldl_sparse, '-o', 'LineWidth', 2);
plot(n_values, diff_lu_dense, '-o', 'LineWidth', 2);
plot(n_values, diff_lu_sparse, '-o', 'LineWidth', 2);
plot(n_values, diff_null_space, '-o', 'LineWidth', 2);
plot(n_values, diff_range_space, '-o', 'LineWidth', 2);

% Set plot labels and legend
xlabel('n');
ylabel('Difference');
title('Differences between Solutions fval and Quadprog fval');
legend('LDL Dense', 'LDL Sparse', 'LU Dense 1', 'LU Dense 2', 'Null Space', 'Range Space');
grid on;


% Plotting the computation times
figure;
plot(n_values, time_quadprog, '-o', 'DisplayName', 'Quadprog');
hold on;
plot(n_values, time_ldl_dense, '-o', 'DisplayName', 'LDL Dense');
plot(n_values, time_ldl_sparse, '-o', 'DisplayName', 'LDL Sparse');
plot(n_values, time_lu_dense, '-o', 'DisplayName', 'LU Dense');
plot(n_values, time_lu_sparse, '-o', 'DisplayName', 'LU Sparse');
plot(n_values, time_null_space, '-o', 'DisplayName', 'Null Space');
plot(n_values, time_range_space, '-o', 'DisplayName', 'Range Space');

xlabel('n');
ylabel('Computation Time');
title('Computation Times for Different Methods');
legend('Location', 'best');
grid on;

function fval = calculateFval(x, H, g)
    fval = 0.5 * x' * H * x + g' * x;
end

