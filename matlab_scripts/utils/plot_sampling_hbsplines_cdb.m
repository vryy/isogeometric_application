%% plot the hierarchical B-Splines sampling points using Cox-de-Boor formula
function Pip=plot_sampling_hbsplines_cdb(Xi,Eta,P,W,params)
min_xi = params.min_xi;
max_xi = params.max_xi;
min_eta = params.min_eta;
max_eta = params.max_eta;
num_points1 = params.num_points1;
num_points2 = params.num_points2;
x_xi = linspace(min_xi,max_xi,num_points1);
x_eta = linspace(min_eta,max_eta,num_points2);
p1 = params.p1;
p2 = params.p2;
Pip = zeros(num_points1*num_points2,3);
n = length(W);
for i = 1:num_points1
     for j = 1:num_points2
         N = zeros(1,n);
         for k = 1:n
             N(k) = CoxDeBoor(x_xi(i),1,p1,Xi{k}) * CoxDeBoor(x_eta(j),1,p2,Eta{k});
         end
         R = W .* N / dot(W,N);
         num = (i-1)*num_points2 + j + 1;
         Pip(num,:) = R * P;
     end
end
scatter(Pip(:,1),Pip(:,2));

% params.num_points1 = 11;
% params.num_points2 = 11;
% params.min_xi = 0.0;
% params.max_xi = 0.999999;
% params.min_eta = 0.0;
% params.max_eta = 0.999999;
% params.p1 = 1;
% params.p2 = 2;
% plot_sampling_hnurbs_cdb(Xi,Eta,P,W,params);

