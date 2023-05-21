n=2;

% Original LP problem
[g, A, b, C, dl, du, l, u] = RandomLP(n);

% Convert inequality constraints to equality constraints
m = size(A, 1);
n = size(A, 2);

% Construct augmented matrices
Cbar = [C, -C, eye(n), -eye(n)];
dbar = [-dl; du; -l; u];

% Introduce slack variables
A_eq = [A, zeros(m, n)];
b_eq = -b;

% Solve the LP problem with equality constraints

% Slack variables
m_eq = size(A_eq, 1);
slack0 = zeros(m_eq, 1);
x0 = zeros(n, 1);

% Construct augmented LP problem
x_bar = [x0; slack0];
b_bar = [b_eq; dbar; zeros(m_eq, 1)]
A_bar = full([A_eq, Cbar])

A_bar'*x_bar + b_bar

% Solve the augmented LP problem
[x2] = linprog(g', [], [], A_bar', -b_bar);


