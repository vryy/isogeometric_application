%% plot the NURBS basis function using Cox-de-Boor formula using the local knot vector
function plot_basis_func_2d_local(Xi_local,Eta_local,p1,p2,num_points)

min_xi = min(Xi_local);
max_xi = max(Xi_local);
xi = linspace(min_xi, max_xi, num_points);
min_eta = min(Eta_local);
max_eta = max(Eta_local);
eta = linspace(min_eta, max_eta, num_points);
v1 = zeros(1,num_points);
v2 = zeros(1,num_points);
for i = 1:num_points
    v1(i) = CoxDeBoor(xi(i),1,p1,Xi_local);
    v2(i) = CoxDeBoor(eta(i),1,p2,Eta_local);
end

[X, Y] = meshgrid(xi, eta);
Z = v1'*v2;
surf(X, Y, Z)

