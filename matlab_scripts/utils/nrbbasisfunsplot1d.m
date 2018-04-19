function [u,Bn] = nrbbasisfunsplot1d(nurbs,num_points)
if ~iscell(nurbs.knots) % curve
    min_knot = min(nurbs.knots);
    max_knot = max(nurbs.knots);
    u = linspace(min_knot,max_knot,num_points);
    [B,id] = nrbbasisfun(u,nurbs);

    %rearrange to the correct sequence
    numbasis = nurbs.number;
    Bn = zeros(size(B,1),numbasis);
    for i = 1:size(B,1)
        for j = 1:size(B,2)
            Bn(i,id(i,j)) = B(i,j);
        end
    end

    plot(u,Bn);
    title(['B-splines basis function plot, n = ' num2str(nurbs.number) ', p = ' num2str(nurbs.order-1)]);
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

