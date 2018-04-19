function nrbbasisfunsplot(nurbs,num_points)
if ~iscell(nurbs.knots) % curve
    nrbbasisfunsplot1d(nurbs,num_points);
else
    if size(nurbs.knots,2) == 2 % surface
        min_knot_1 = min(nurbs.knots{1});
        max_knot_1 = max(nurbs.knots{1});
        u = linspace(min_knot_1,max_knot_1,num_points);
        min_knot_2 = min(nurbs.knots{2});
        max_knot_2 = max(nurbs.knots{2});
        v = linspace(min_knot_2,max_knot_2,num_points);
        numbasis = nurbs.number(1)*nurbs.number(2);

        X = zeros(num_points,num_points);
        Y = zeros(num_points,num_points);
        Z = zeros(num_points,num_points,numbasis);
        [B,id] = nrbbasisfun({u,v},nurbs);
        for i = 1:num_points
            for j = 1:num_points
                X(i,j) = u(i);
                Y(i,j) = v(j);
            end
        end
%        size(B)
%        size(id)
        for i = 1:size(B,1) %num_points*num_points
            for j = 1:size(B,2) %(p+1)*(q+1)
                [row,col] = ind2sub([num_points,num_points],i);
                Z(row,col,id(i,j)) = B(i,j);
            end
        end

%        X = zeros(num_points,num_points);
%        Y = zeros(num_points,num_points);
%        Z = zeros(num_points,num_points,numbasis);
%        for i = 1:length(u)
%            for j = 1:length(v)
%                X(i,j) = u(i);
%                Y(i,j) = v(j);
%                [B,id] = nrbbasisfun([u(i);v(j)],nurbs);
%                for k = 1:length(B)
%                    Z(i,j,id(k)) = B(k);
%                end
%            end
%        end
        for k = 1:numbasis %note that the basis function is arranged in u-v direction
            surf(X,Y,Z(:,:,k));
        end
    elseif size(nurbs.knots,2) == 3 % volume
        % TODO
    else
        error('Invalid knot vector of NURBS');
    end
end

% wrong evaluation of nurbs basis functions (because of the id problem)
% U = [0 0 0 0 0.5 1 1 1 1];
% x = [0 1/3 0.5 2/3 1];
% y = [0 0 0 0 0];
% w = [1 1 1 1 1];
% nrb = nrbmak ([x;y;y;w], U);
% u = linspace(0, 1, 30);
% B = nrbbasisfun(u, nrb);
% plot(u, B)
% title('Cubic Bernstein polynomials')

% p = 2;   q = 3;   m = 4; n = 5;
% Lx  = 1; Ly  = 1; 
% nrb = nrb4surf   ([0 0], [1 0], [0 1], [1 1]);
% nrb = nrbdegelev (nrb, [p-1, q-1]);
% aux1 = linspace(0,1,m); aux2 = linspace(0,1,n);
% nrb = nrbkntins  (nrb, {aux1(2:end-1), aux2(2:end-1)});
% u = rand (1, 30); v = rand (1, 10);
% u = u - min (u); u = u / max (u);
% v = v - min (v); v = v / max (v);
% u
% v
% nrb
% [B, N] = nrbbasisfun ({0}, nrb);
