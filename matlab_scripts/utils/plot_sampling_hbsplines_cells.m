%% plot the hierarchical B-Splines sampling points using Cox-de-Boor formula
function plot_sampling_hbsplines_cells(Xi,Eta,P,W,Id,S,C,N,params)
Pts = zeros(length(S),2);
for c = 1:length(S)
    Pts(1,1) = S{c}(1,1);
    Pts(1,2) = S{c}(2,1);
    Pts(2,1) = S{c}(1,2);
    Pts(2,2) = S{c}(2,1);
    Pts(3,1) = S{c}(1,2);
    Pts(3,2) = S{c}(2,2);
    Pts(4,1) = S{c}(1,1);
    Pts(4,2) = S{c}(2,2);
%    for i = 1:size(Pts,1)
%        if abs(Pts(i,1) - params.max_xi) < 1.0e-6
%            Pts(i,1) = 0.9999*params.max_xi;
%        end
%        if abs(Pts(i,2) - params.max_eta) < 1.0e-6
%            Pts(i,2) = 0.9999*params.max_eta;
%        end
%    end
%    Pip = hnurbs_cdb_at(Pts,Xi,Eta,P,W,params);
    Pip = hnurbs_bezier_at(Pts,P,W,Id,S,C,N,params);
    line([Pip(1,1) Pip(2,1)],[Pip(1,2) Pip(2,2)]);
    line([Pip(2,1) Pip(3,1)],[Pip(2,2) Pip(3,2)]);
    line([Pip(3,1) Pip(4,1)],[Pip(3,2) Pip(4,2)]);
    line([Pip(4,1) Pip(1,1)],[Pip(4,2) Pip(1,2)]);
end

function Pip=hnurbs_cdb_at(Pts,Xi,Eta,P,W,params)
p1 = params.p1;
p2 = params.p2;
num_points = size(Pts,1);
Pip = zeros(num_points,3);
n = length(W);
N = zeros(1,n);
for i = 1:num_points
    for k = 1:n
        N(k) = CoxDeBoor(Pts(i,1),1,p1,Xi{k}) * CoxDeBoor(Pts(i,2),1,p2,Eta{k});
    end
    R = W .* N / dot(W,N);
    Pip(i,:) = R * P;
end

function Pip = hnurbs_bezier_at(Pts,P,W,Id,S,C,N,params)
% check if the point lie in which cell
num_points = size(Pts,1);
p1 = params.p1;
p2 = params.p2;
for i = 1:num_points
    xi = Pts(i,1);
    eta = Pts(i,2);
    found = 0;
    for c = 1:length(N)
        if xi>=S{c}(1,1) && xi<=S{c}(1,2) && eta>=S{c}(2,1) && eta<=S{c}(2,2)
            e = c;
            found = 1;
            break;
        end
    end
    if found == 0
        error(['the sampling point is not found']);
    end

    % extract the weight and nodal coordinates
    We = zeros(length(N{e}),1);
    Pe = zeros(length(N{e}),3);
    for n = 1:length(N{e})
        for b = 1:length(W)
            if Id(b) == N{e}(n)
                We(n) = W(b);
                Pe(n,:) = P(b,:);
            end
        end
        end

    % compute the local point
    local_xi = (xi-S{e}(1,1)) / (S{e}(1,2)-S{e}(1,1));
    local_eta = (eta-S{e}(2,1)) / (S{e}(2,2)-S{e}(2,1));
    
    % compute the Bezier shape function at local point
    B = zeros((p1+1)*(p2+1),1);
    for k = 0:p1
        y_xi = bezier(local_xi,k,p1);
        for l = 0:p2
            y_eta  = bezier(local_eta,l,p2);
            num = k*(p2+1)+l+1;
            B(num) = y_xi * y_eta;
        end
    end

    % compute the NURBS shape function
    Ne = C{e} * B;
    Re = We .* Ne / dot(We,Ne);
    Pip(i,:) = Re' * Pe;
end

% params.num_points1 = 11;
% params.num_points2 = 11;
% params.min_xi = 0.0;
% params.max_xi = 0.999999;
% params.min_eta = 0.0;
% params.max_eta = 0.999999;
% params.p1 = 1;
% params.p2 = 2;
% plot_sampling_hnurbs_cdb(Xi,Eta,P,W,params);

