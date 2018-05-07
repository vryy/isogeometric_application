%% plot the hierarchical B-Splines sampling points using Bezier discretization
function Pip=plot_sampling_hbsplines_bezier(P,W,Id,S,C,N,params)
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

% create a list of sampling points
x_xi = linspace(min_xi,max_xi,num_points1);
x_eta = linspace(min_eta,max_eta,num_points2);
Pip = zeros(length(x_xi)*length(x_eta),3);
for i = 1:length(x_xi)
    for j = 1:length(x_eta)
        % check if the point lie in which cell
        xi = x_xi(i);
        eta = x_eta(j);
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
        num = (i-1)*num_points2 + j + 1;
        Pip(num,:) = Re' * Pe;
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
% plot_sampling_hnurbs_bezier(P,W,Id,S,C,N,params);

