%% plot the hierarchical B-Splines cell based on sampling points
function plot_sampling_hbsplines_cell_2d(cell_id,Xi,Eta,P,W,Id,S,C,N,CId,params)
nsampling1 = params.num_points1;
nsampling2 = params.num_points2;

for c = 1:length(S)
    if CId{c} == cell_id
        cnt = 1;
        for i = 0:nsampling1
            for j = 0:nsampling2
                Pts(cnt,1) = S{c}(1,1) + i*(S{c}(1,2) - S{c}(1,1)) / nsampling1;
                Pts(cnt,2) = S{c}(2,1) + j*(S{c}(2,2) - S{c}(2,1)) / nsampling2;
                if i == 0 && j == 0
                    c1 = cnt;
                elseif i == nsampling1 && j == 0
                    c2 = cnt;
                elseif i == nsampling1 && j == nsampling2
                    c3 = cnt;
                elseif i == 0 && j == nsampling2
                    c4 = cnt;
                end
                cnt = cnt + 1;
            end
        end

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
%        Pts(c1,:)
%        Pts(c2,:)
%        Pts(c3,:)
%        Pts(c4,:)
        if CId{c} == 69
            c
            Pts = [0.0694318442029737*0.03125, 0.0694318442029737*0.03125]
            N{c}
        end

        if strcmp(params.method, 'cdb')
            Pip = hbsplines_cdb_at(Pts,Xi,Eta,P,W,params);
        elseif strcmp(params.method, 'bezier')
            Pip = hbsplines_bezier_at(Pts,P,W,Id,S,C,N,params);
        end
%            Pip(c1,:)
%            Pip(c2,:)
%            Pip(c3,:)
%            Pip(c4,:)

        line([Pip(c1,1) Pip(c2,1)],[Pip(c1,2) Pip(c2,2)]);
        line([Pip(c2,1) Pip(c3,1)],[Pip(c2,2) Pip(c3,2)]);
        line([Pip(c3,1) Pip(c4,1)],[Pip(c3,2) Pip(c4,2)]);
        line([Pip(c4,1) Pip(c1,1)],[Pip(c4,2) Pip(c1,2)]);
        hold on

        scatter(Pip(:,1), Pip(:,2));
    end
end

