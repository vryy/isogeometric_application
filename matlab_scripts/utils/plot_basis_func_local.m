%% plot the NURBS basis function using Cox-de-Boor formula using the local knot vector
function plot_basis_func_local(Xi_local,p,num_points)

min_xi = min(Xi_local);
max_xi = max(Xi_local);
xi = linspace(min_xi, max_xi, num_points);
v = zeros(1,num_points);
for i = 1:num_points
    v(i) = CoxDeBoor(xi(i),1,p,Xi_local);
end

plot(xi, v);

