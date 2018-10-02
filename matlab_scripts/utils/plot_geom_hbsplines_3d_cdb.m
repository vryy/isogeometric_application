%% plot the hierarchical NURBS geometry using Cox-de-Boor formula
function plot_geom_hbsplines_3d_cdb(Xi,Eta,Zeta,P,W,params)
min_xi = params.min_xi;
max_xi = params.max_xi;
min_eta = params.min_eta;
max_eta = params.max_eta;
min_zeta = params.min_zeta;
max_zeta = params.max_zeta;
num_points1 = params.num_points1;
num_points2 = params.num_points2;
num_points3 = params.num_points3;
x_xi = linspace(min_xi,max_xi,num_points1);
x_eta = linspace(min_eta,max_eta,num_points2);
x_zeta = linspace(min_zeta,max_zeta,num_points3);
X = zeros(num_points1*num_points2*num_points3,1);
Y = zeros(num_points1*num_points2*num_points3,1);
Z = zeros(num_points1*num_points2*num_points3,1);
p1 = params.p1;
p2 = params.p2;
p3 = params.p3;
for i = 1:num_points1
    for j = 1:num_points2
        for l = 1:num_points3
            Denom = 0.0;
            num = (((i-1)*num_points2+j-1)*num_points3)+l;
            for k = 1:length(W)
                N = W(k) * CoxDeBoor(x_xi(i),1,p1,Xi{k}) * CoxDeBoor(x_eta(j),1,p2,Eta{k}) * CoxDeBoor(x_zeta(l),1,p3,Zeta{k});
                Denom = Denom + N;
                X(num) = X(num) + N * P(k,1);
                Y(num) = Y(num) + N * P(k,2);
                Z(num) = Z(num) + N * P(k,3);
            end
            X(num) = X(num) / Denom;
            Y(num) = Y(num) / Denom;
            Z(num) = Z(num) / Denom;
        end
    end
end
scatter3(X,Y,Z);

% params.num_points1 = 11;
% params.num_points2 = 11;
% params.num_points3 = 4;
% params.min_xi = 0.0;
% params.max_xi = 0.999999;
% params.min_eta = 0.0;
% params.max_eta = 0.999999;
% params.min_zeta = 0.0;
% params.max_zeta = 0.999999;
% params.p1 = 2;
% params.p2 = 2;
% params.p3 = 1;
% plot_geom_hnurbs_3d_cdb(Xi,Eta,Zeta,P,W,params);

