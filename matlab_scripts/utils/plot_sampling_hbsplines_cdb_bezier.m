%% plot the hierarchical B-Splines sampling points using Cox-de-Boor formula and Bezier discretization for comparison
function Pip=plot_sampling_hbsplines_cdb_bezier(Xi,Eta,P,W,Id,S,C,N,params)
clc
close all
hold on
axis equal
format long
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
Pip_cdb = zeros(num_points1*num_points2,3);
Pip_bezier = zeros(num_points1*num_points2,3);
for i = 1:num_points1
    for j = 1:num_points2
        xi = x_xi(i);
        eta = x_eta(j);

        % compute the basis function based on Cox-de-Boor formula
        Bs = zeros(1,length(W));
        for k = 1:length(W)
            Bs(k) = CoxDeBoor(xi,1,p1,Xi{k}) * CoxDeBoor(eta,1,p2,Eta{k});
        end
        R = W .* Bs / dot(W,Bs);
        num = (i-1)*num_points2 + j + 1;
        Pip_cdb(num,:) = R * P;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % compute the basis function based on Bezier discretization
        % check if the point lie in which cell
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
        Bz = zeros((p1+1)*(p2+1),1);
        for k = 0:p1
            y_xi = bezier(local_xi,k,p1);
            for l = 0:p2
                y_eta  = bezier(local_eta,l,p2);
                num = k*(p2+1)+l+1;
                Bz(num) = y_xi * y_eta;
            end
        end

        % compute the NURBS shape function
        Ne = C{e} * Bz;
        Re = We .* Ne / dot(We,Ne);
        num = (i-1)*num_points2 + j + 1;
        Pip_bezier(num,:) = Re' * Pe;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
        % compare and check why is the difference
        dif = norm(Pip_cdb(num,:) - Pip_bezier(num,:));
        if dif > 1.0e-6
            disp(['Warning: there are difference in computing basis function using Bezier discretization']);
            disp(['there are problems at Bezier element ' num2str(e)]);
            dif
            xi
            eta
            local_xi
            local_eta
            N{e}
            C{e}
            
            % compute all the nonzero basis function at this point using Cox-de-Boor formula
            cnt = 1;
            for k = 1:length(W)
                tmp = CoxDeBoor(xi,1,p1,Xi{k}) * CoxDeBoor(eta,1,p2,Eta{k});
                if abs(tmp) > 1.0e-6
                    nonzeroNid(cnt) = Id(k);
                    nonzeroNval(cnt) = tmp;
                    cnt = cnt + 1;
                end
            end
            nonzeroNid
            nonzeroNval' - Ne
            clear nonzeroNid nonzeroNval
            
            disp(['---------------------------------']);
        end
    end
end
scatter(Pip_cdb(:,1),Pip_cdb(:,2),'x');
scatter(Pip_bezier(:,1),Pip_bezier(:,2),'o');

% params.num_points1 = 11;
% params.num_points2 = 11;
% params.min_xi = 0.0;
% params.max_xi = 0.999999;
% params.min_eta = 0.0;
% params.max_eta = 0.999999;
% params.p1 = 1;
% params.p2 = 2;
% plot_sampling_hnurbs_cdb_bezier(Xi,Eta,P,W,Id,S,C,N,params);

