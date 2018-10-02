%% plot the hierarchical NURBS geometry using Cox-de-Boor formula
function plot_geom_hbsplines_2d_cdb(Xi,Eta,P,W,params)
min_xi = params.min_xi;
max_xi = params.max_xi;
min_eta = params.min_eta;
max_eta = params.max_eta;
num_points1 = params.num_points1;
num_points2 = params.num_points2;
x_xi = linspace(min_xi,max_xi,num_points1);
x_eta = linspace(min_eta,max_eta,num_points2);
X = zeros(num_points1,num_points2);
Y = zeros(num_points1,num_points2);
Z = zeros(num_points1,num_points2);
p1 = params.p1;
p2 = params.p2;
for i = 1:num_points1
     for j = 1:num_points2
         Denom = 0.0;
         for k = 1:length(W)
             N = W(k) * CoxDeBoor(x_xi(i),1,p1,Xi{k}) * CoxDeBoor(x_eta(j),1,p2,Eta{k});
             Denom = Denom + N;
             X(i,j) = X(i,j) + N * P(k,1);
             Y(i,j) = Y(i,j) + N * P(k,2);
             Z(i,j) = Z(i,j) + N * P(k,3);
         end
         X(i,j) = X(i,j) / Denom;
         Y(i,j) = Y(i,j) / Denom;
         Z(i,j) = Z(i,j) / Denom;
     end
end
surf(X,Y,Z);

% params.num_points1 = 11;
% params.num_points2 = 11;
% params.min_xi = 0.0;
% params.max_xi = 0.999999;
% params.min_eta = 0.0;
% params.max_eta = 0.999999;
% params.p1 = 1;
% params.p2 = 2;
% plot_geom_hnurbs_cdb(Xi,Eta,P,W,params);

