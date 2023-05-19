n_values = 25:25:100;
num_runs=10;
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

    [H, g, A, b] = GenerateRandomEQP(n, 0.2, 0.2);   
    tic;
    [x_true, fval_true] = quadprog(H, g', [], [], A, b);
    time_quadprog(i) = toc;
    
    for run = 1:num_runs

        mse_ldl_dense = zeros(num_runs, 1);
        mse_ldl_sparse = zeros(num_runs, 1);
        mse_lu_dense = zeros(num_runs, 1);
        mse_lu_sparse = zeros(num_runs, 1);
        mse_null_space = zeros(num_runs, 1);
        mse_range_space = zeros(num_runs, 1);

        tic;
        [x1] = EQPsolver(H, g, A, b, "LDLdense");
        fval1 = calculateFval(x1, H, g);
        time = toc;

        if run == 1
           time_ldl_dense(i) = time;
        end

        mse_ldl_dense(run) = ((fval1 - fval_true).^2);
    
        tic;
        [x2] = EQPsolver(H, g, A, b, "LDLsparse");
        fval2 = calculateFval(x2, H, g);
        time = toc;

        if run == 1
            time_ldl_sparse(i) = time;
        end

        mse_ldl_sparse(run) = ((fval2 - fval_true).^2);
    
        tic;
        [x3] = EQPsolver(H, g, A, b, "LUdense");
        fval3 = calculateFval(x3, H, g);
        time = toc;
        
        if run == 1
            time_lu_dense(i) = time;
        end

        mse_lu_dense(run) = ((fval3 - fval_true).^2);
    
        tic;
        [x4] = EQPsolver(H, g, A, b, 'LUsparse');
        fval4 = calculateFval(x4, H, g);
        time = toc;
        
        if run == 1
            time_lu_sparse(i) = time;
        end

        mse_lu_sparse(run) = ((fval4 - fval_true).^2);
    
        tic;
        [x5] = EQPsolver(H, g, A, b, "NullSpace");
        fval5 = calculateFval(x5, H, g);
        time = toc;
        
        if run == 1
            time_null_space(i) = time;
        end

        mse_null_space(run) = ((fval5 - fval_true).^2);
    
        tic;
        [x6] = EQPsolver(H, g, A, b, "RangeSpace");
        fval6 = calculateFval(x6, H, g);
        time = toc;
        
        if run == 1
            time_range_space(i) = time;
        end
        mse_range_space(run) = ((fval6 - fval_true).^2);
    end

    diff_ldl_dense(i) = sum(mse_ldl_dense)/num_runs;
    diff_ldl_sparse(i) = sum(mse_ldl_sparse)/num_runs;
    diff_lu_dense(i) = sum(mse_lu_dense)/num_runs;
    diff_lu_sparse(i) = sum(mse_lu_sparse)/num_runs;
    diff_null_space(i) = sum(mse_null_space)/num_runs;
    diff_range_space(i) = sum(mse_range_space)/num_runs;

end

% Plot the differences for each method
figure;
plot(n_values, diff_ldl_dense, '-o');
hold on;
plot(n_values, diff_ldl_sparse, '-o');
plot(n_values, diff_lu_dense, '-o');
plot(n_values, diff_lu_sparse, '-o');
plot(n_values, diff_null_space, '-o');
plot(n_values, diff_range_space, '-o');

% Set plot labels and legend
xlabel('Problem size, n', 'FontSize', 16);
ylabel('MSE', 'FontSize', 16);
title('MSE of 20 runs per problem size "n"', 'FontSize', 16);
legend('LDL Dense', 'LDL Sparse', 'LU Dense', 'LU Sparse', 'Null Space', 'Range Space');
% Get the handle to the legend object
hLegend = legend;

% Increase the font size and other properties of the legend
set(hLegend, 'FontSize', 14)
grid on;


% Plotting the computation times
figure;
plot(n_values, time_quadprog, '-o');
hold on;
plot(n_values, time_ldl_dense, '-o');
plot(n_values, time_ldl_sparse, '-o');
plot(n_values, time_lu_dense, '-o');
plot(n_values, time_lu_sparse, '-o');
plot(n_values, time_null_space, '-o');
plot(n_values, time_range_space, '-o');

xlabel('Problem size, n', 'FontSize', 16);
ylabel('Computation Time, s', 'FontSize', 16);
title('Computation Times of a range of problem sizes', 'FontSize', 16);

legend('Quadprog', 'LDL Dense', 'LDL Sparse', 'LU Dense', 'LU Sparse', 'Null Space', 'Range Space');
% Get the handle to the legend object
hLegend = legend;

% Increase the font size and other properties of the legend
set(hLegend, 'FontSize', 14)
grid on;

% Adjust the figure properties for better readability
set(gcf, 'Position', [100, 100, 800, 500]);  % Set figure size and position
set(gca, 'FontName', 'Arial', 'FontSize', 12);  % Set axis font and font size
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');  % Add minor ticks to both axes
box on;  % Add a box around the plot area


function fval = calculateFval(x, H, g)
    fval = 0.5 * x' * H * x + g' * x;
end

