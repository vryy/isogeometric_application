
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def bspeval(d, c, k, u):

    # Local Variables: c, d, mc, i, k, nc, N, p, s, u, tmp1, nu
    # Function calls: bspeval, zeros, numel, basisfun, repmat, findspan, size
    #% BSPEVAL:  Evaluate B-Spline at parametric points.
    #% 
    #% Calling Sequence:
    #% 
    #%   p = bspeval(d,c,k,u)
    #% 
    #%    INPUT:
    #% 
    #%       d - Degree of the B-Spline.
    #%       c - Control Points, matrix of size (dim,nc).
    #%       k - Knot sequence, row vector of size nk.
    #%       u - Parametric evaluation points, row vector of size nu.
    #% 
    #%    OUTPUT:
    #%
    #%       p - Evaluated points, matrix of size (dim,nu)
    #% 
    #%    Copyright (C) 2000 Mark Spink, 2007 Daniel Claxton, 2010 C. de Falco
    #%
    #%    This program is free software: you can redistribute it and/or modify
    #%    it under the terms of the GNU General Public License as published by
    #%    the Free Software Foundation, either version 3 of the License, or
    #%    (at your option) any later version.
    #%    This program is distributed in the hope that it will be useful,
    #%    but WITHOUT ANY WARRANTY; without even the implied warranty of
    #%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #%    GNU General Public License for more details.
    #%
    #%    You should have received a copy of the GNU General Public License
    #%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    nu = numel(u)
    [mc, nc] = matcompat.size(c)
    #%   int bspeval(int d, double *c, int mc, int nc, double *k, int nk, double *u,int nu, double *p){
    #%   int ierr = 0;
    #%   int i, s, tmp1, row, col;
    #%   double tmp2;
    #%
    #%   // Construct the control points
    #%   double **ctrl = vec2mat(c,mc,nc);
    #%
    #%   // Contruct the evaluated points
    #%   double **pnt = vec2mat(p,mc,nu);
    #%
    #%   // space for the basis functions
    #%N = zeros(d+1,1);                               %   double *N = (double*) mxMalloc((d+1)*sizeof(double));
    #%
    #%   // for each parametric point i
    #%for col=1:nu                                    %   for (col = 0; col < nu; col++) {
    #%     // find the span of u[col]
    s = findspan((nc-1.), d, u.flatten(1), k)
    #%     s = findspan(nc-1, d, u[col], k);
    N = basisfun(s, u.flatten(1), d, k)
    #%     basisfun(s, u[col], d, k, N);
    #%
    tmp1 = s-d+1.
    #%     tmp1 = s - d;
    #%for row=1:mc                                %     for (row = 0; row < mc; row++)  {
    p = np.zeros(mc, nu)
    #%       tmp2 = 0.0;
    for i in np.arange(0., (d)+1):
        #%       for (i = 0; i <= d; i++)
        
    #%
    #%p(row,:) = tmp2;                      %       pnt[col][row] = tmp2;
    #%end                                         %     }
    #%end                                             %   }
    #%
    #%   mxFree(N);
    #%   freevec2mat(pnt);
    #%   freevec2mat(ctrl);
    #%
    #%   return ierr;
    #%   }
    return [p]