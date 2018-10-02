%% plot the hierarchical B-Splines sampling points using Cox-de-Boor formula
function plot_sampling_hbsplines_cells_2d(Xi,Eta,P,W,Id,S,C,N,params)
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

    % Pts
    if isfield(params,'adjust')
        for i = 1:size(Pts,1)

            if abs(Pts(i,1) - 0.0) < params.tol
                Pts(i,1) = params.tol;
            end
            if abs(Pts(i,2) - 0.0) < params.tol
                Pts(i,2) = params.tol;
            end

            if abs(Pts(i,1) - params.max_xi) < params.tol
                Pts(i,1) = 0.9999*params.max_xi;
            end
            if abs(Pts(i,2) - params.max_eta) < params.tol
                Pts(i,2) = 0.9999*params.max_eta;
            end
        end
    end
%    Pts

    if strcmp(params.method, 'cdb')
        Pip = hbsplines_cdb_at(Pts,Xi,Eta,P,W,params);
    elseif strcmp(params.method, 'bezier')
        Pip = hbsplines_bezier_at(Pts,P,W,Id,S,C,N,params);
    end
%    Pip

    line([Pip(1,1) Pip(2,1)],[Pip(1,2) Pip(2,2)]);
    line([Pip(2,1) Pip(3,1)],[Pip(2,2) Pip(3,2)]);
    line([Pip(3,1) Pip(4,1)],[Pip(3,2) Pip(4,2)]);
    line([Pip(4,1) Pip(1,1)],[Pip(4,2) Pip(1,2)]);
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

