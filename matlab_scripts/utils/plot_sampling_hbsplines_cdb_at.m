%% plot the hierarchical B-Splines sampling points using Cox-de-Boor formula
function Pip=plot_sampling_hbsplines_cdb_at(Pts,Xi,Eta,P,W,params)
p1 = params.p1;
p2 = params.p2;
num_points = size(Pts,1);
Pip = zeros(num_points,3);
n = length(W);
N = zeros(1,n);
for i = 1:num_points
    for k = 1:n
        N(k) = CoxDeBoor(Pts(i,1),1,p1,Xi{k}) * CoxDeBoor(Pts(i,),1,p2,Eta{k});
    end
    R = W .* N / dot(W,N);
    Pip(i,:) = R * P;
end
scatter(Pip(:,1),Pip(:,2));

% params.p1 = 1;
% params.p2 = 2;
% plot_sampling_hnurbs_cdb(Xi,Eta,P,W,params);

