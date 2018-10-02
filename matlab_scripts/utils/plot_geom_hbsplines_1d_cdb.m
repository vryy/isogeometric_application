%% plot the hierarchical NURBS geometry using Cox-de-Boor formula
function plot_geom_hbsplines_1d_cdb(Xi,P,W,params)
min_xi = params.min_xi;
max_xi = params.max_xi;
num_points1 = params.num_points1;
x_xi = linspace(min_xi,max_xi,num_points1);
X = zeros(num_points1);
Y = zeros(num_points1);
Z = zeros(num_points1);
p1 = params.p1;
for i = 1:num_points1
     Denom = 0.0;
     for k = 1:length(W)
         N = W(k) * CoxDeBoor(x_xi(i),1,p1,Xi{k});
         Denom = Denom + N;
         X(i) = X(i) + N * P(k,1);
         Y(i) = Y(i) + N * P(k,2);
         Z(i) = Z(i) + N * P(k,3);
     end
     X(i) = X(i) / Denom;
     Y(i) = Y(i) / Denom;
     Z(i) = Z(i) / Denom;
end
plot3(X,Y,Z);

% params.num_points1 = 11;
% params.min_xi = 0.0;
% params.max_xi = 0.999999;
% params.p1 = 1;
% plot_geom_hnurbs_cdb(Xi,P,W,params);

