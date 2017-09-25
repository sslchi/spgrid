function [x, w, v] = glpts( n )
%GLPTS gives the Gauss–Lobatto points and quadrature/barycentric weights.
% [xk, wk, vk] = GLPTS( n )

% $Date: 24-Aug-2017$
% Copyright (c) Guanjie Wang, wangguanjie0@126.com

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Gauss-Lobatto popints are Gauss points for weights (1-x)(1+x), and can
% computed by JACOBIPTS. The barycentric weights can be computed by
% BARYWEIGHTS.The quadrature weights can be computed by 
% wi = 2/(n(n-1)*p_{n-1}^2(x_i)), where x_i neq +- 1. p_n(x_i) is computed by
% the recursion formula. 
%
% Note that for x_i ~ = 1 and x_i ~= -1,
% int (l_i(x)) dx = 1/(x_i^2-1) int (x^2-1) \hat{l}_i(x) dx = w_i/(1-x_i^2)
% \lambda_i = 1/(x_i^2-1)/prod_{k~=i}(x_i-x_k) = v_i/(x_i^2-1).
%
% Refernce:
% [1]:https://en.wikipedia.org/wiki/Gaussian_quadrature
% [2]:https://en.wikipedia.org/wiki/Legendre_polynomials
% [3]:Walter Gautschi, High-order Gauss–Lobatto formulae,Numerical Algorithms,
%     25: 213–222, 2000.
% 
% Historical note:
% 24-Aug-2017: compute by recursive procedure.
% 24-Aug-2017: by Nick Trefethen's methods.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note:
%   The approach used here is to observe that the Gauss-Lobatto points are
%   precisely the roots of (1-x^2)P'_{n-1}(x), and that the roots of P'_{n-1}(x)
%   are the same as the roots of P^(1,1)_{n-2}(x) [NIST, (18.9.15)], which can
%   be obtained for JACPTS. A similar identity [NIST, (18.9.16)] is used for the
%   computation of the quadrature weights from those of JACPTS, and the missing
%   barycentric weights are determined by enforcing the interpolation of f(x) =
%   x or x^2 at x = 0 in the even or odd case respectively.
%
%    x_j = roots of (1-x^2)P'_{n-1}(x)
%    w_j = { 2/(n*(n-1))                        : x_j = -1, 1
%          { 2/(n*(n-1)) * 1/[P_{n-1}(x_j)]^2   : otherwise
%
%   (Note that the weights for n-2 point Gauss-Jacobi with a = b = 1 satisfy 
%    u_j = C/(1-x_j^2)/[d/dx P^(1,1)_{n-2}(x_j)]^2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Trivial cases:
if ( n == 1 )
    error('CHEBFUN:lobpts:notSupported', 'N = 1 is not supported.');
elseif ( n == 2 )
    x = [-1 ; 1];
    w = [1, 1];
    v = [-1 ; 1];
    return
elseif ( n == 3 )
    x = [-1 ; 0 ; 1];
    w = [1, 4, 1]/3;
    v = [-.5 ; 1 ; -.5];
    return
end


%% Call JACPTS():
[x, w, v] = jacobipts(n - 2, 1, 1);

%% Nodes:
x = [-1 ; x ; 1];

%% Quadrature weights: 
w = [-1, w,  1];
w = w./(1-x.^2).';
w([1 end]) = 2/(n*(n - 1));

%% Barycentric weights:
v = v./(1 - x(2:n-1).^2);
v = v/max(abs(v));
if ( mod(n, 2) )
    v1 = -abs(sum(v.*x(2:end-1).^2)/2);
    sgn = 1;
else
    v1 = -abs(sum(v.*x(2:end-1))/2);
    sgn = -1;
end
v = -1 * [v1 ; v ; sgn*v1]; % make the first barycentric weight positive

end


% function [xk, wk, vk] = glpts( n )
% if n < 2
%     error('orthogonal:glpts',' input sholud be >= 2.')
% end
% 
% xk = jacobipts(n-2, 1,1); xk = [-1;xk;1]; 
% c = legpoly(n-1,xk); c(1) = 1; c(end) = 1; wk = 2./(n*(n-1).*c.^2); wk = wk.';
% 
% vk = baryweights(xk);
% 
% end
% 
% function c = legpoly(n,x)
% 
% a = 1; b = x;
% 
% switch n
%     case 0
%         c = 1;
%     case 1
%         c = x;
%     otherwise
%         for i = 2:n
%             c = (2*i-1).*x.*b - (i-1)*a;
%             c = c/i;
%             a = b; b= c;
%         end
% end
% 
% 
% end