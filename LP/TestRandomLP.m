n=20;
% Original LP problem
[g, A, b, C, dl, du, l, u] = RandomLP(n);
linprog(g,A',b,C',dl-du,l,u)

