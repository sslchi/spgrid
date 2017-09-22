function [x, w, v] = laguerrepts(n, alpha)
%LAGURREPTS  Give zero points and weights for quadrature of gauss type.
%  
%  [x, w] = LAGURREPTS(n, alphabeta) give the zero points of p_n and weights
%  for quadrature of gauss type. The weight function is 
%  w =  x^alpha exp(-x), where alpha > -1. n must be smaller than 2000, for
%  n> 2000, you can use the matlab package Chebfun, chebfun/lagpts, www.chebfun.org.
%  Note that the weights in this routine is NOT THE SAME as that of glaguerre. 
% 
% See also glaguerre.

% $Date:2017/04/25$
% Copyright (c) G.Wang. 
% Email: wangguanjie0@126.com


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This program use the Jacobin matrix and matlab function eig to compute
%  the zero points and weigts. The cost is about O(n^3) due to eig,  
%  and can be improved to O(n^2). There are methods at cost about O(n),
%  see 
%  [1]: chebfun/lagpts, www.chebfun.org.

% The reccurenc coefficients can be found in
%  [2]: Dongbin Xiu, "Numerical methods for stochastic computations", 2010
%  [3]: Wikipedia, https://en.wikipedia.org/wiki/Laguerre_polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


a = 0;% set default value

% Parse input
if ( (nargin == 2) && isnumeric(alpha) )
 a = alpha;
end

if ( n > 2000 )
    error('Too many gauss points, use Chebfun to obtain.')
end

[x, w, v] = lagurrepts_gw(n, a);

end

%% ------------------------- Routines for GW ----------------------------
function [x, w, v] = lagurrepts_gw(n, alpha)

alph = 2*(1:n)-1 + alpha;  % 3-term recurrence coeffs
beta = sqrt( (1:n-1).*(alpha + (1:n-1) ) );
T = diag(beta,1) + diag(alph) + diag(beta,-1);  % Jacobi matrix
[V, D] = eig(T);                      % eigenvalue decomposition
[x, indx] = sort(diag(D));            % Laguerre points
w = gamma(alpha+1)*V(1,indx).^2;      % Quadrature weights
v = sqrt(x).*abs(V(1,indx)).';        % Barycentric weights
v = v./max(v);
v(2:2:n) = -v(2:2:n);

end