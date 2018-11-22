% computing B-spline basis function based on Cox-de-Boor algorithm
function y = CoxDeBoor(u,i,p,knots)
% Input:
%   u       knot to be compute the function value
%   i       index of the basis function; for local knot vector, it must be 1
%   p       B-spline degree
%   knots   local knot vector
% Output: function value
y = 0.0;

if p == 0
%     disp(['knots i:' num2str(knots(i))]);
%     disp(['knots i+1:' num2str(knots(i+1))]);
%     disp(['u:' num2str(u)]);
%    if u >= knots(i) && u < knots(i+1)
    if u >= knots(i) && u < knots(i+1) + 1.0e-10 % at a very small value to enhance the accuracy at right boundary
%    if u >= knots(i) && u < knots(i+1) + 1.0e-10*(knots(i+1)-knots(i)) % at a very small value to enhance the accuracy at right boundary
        y = 1.0;
    else
        y = 0.0;
    end
    return;
end

if abs(knots(i+p) - knots(i)) > 1.0e-6
    y = y + (u - knots(i)) / (knots(i+p) - knots(i)) *...
                CoxDeBoor(u,i,p-1,knots);
end

if abs(knots(i+p+1) - knots(i+1)) > 1.0e-6
    y = y + (knots(i+p+1) - u) / (knots(i+p+1) - knots(i+1)) *...
                CoxDeBoor(u,i+1,p-1,knots);
end

