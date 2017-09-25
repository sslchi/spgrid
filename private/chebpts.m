function [x,w,v] = chebpts(n)
%CHEBPTS compute the second kind chebyshev points.
%   [X,W,V] = CHEBPTS(N) retuns N Chebyshev points of the 1st kind in [-1, 1].
%   X is the second kind Chebshev points, W is the weights of Clenshaw-Curtis
%   quadrature, and V is the for barycentric polynomial interpolation in the
%   Chebyshev points X.
%

% This programe is adapted from Chebfun/chebtech2.chebpts, http://www.chebfun.org/


if n == 1
    x = 0;
    w = 2;
    v = 1;
else
    % Chebyshev points:
    m = n - 1;
    x = sin(pi*(-m:2:m)/(2*m)).';  % (Use of sine enforces symmetry.)
    
    % quadratrue weights
    c = 2./[1, 1-(2:2:(m)).^2];    % Exact integrals of T_k (even)
    c = [c, c(floor(n/2):-1:2)];   % Mirror for DCT via FFT
    w = ifft(c);                   % Interior weights
    w([1,n]) = w(1)/2;             % Boundary weights
    
    % barycentric interpolate weights
    v = [0.5;ones(n-1,1)];        % Note v(1) is positive.
    v(2:2:end) = -1;
    v(end) = .5*v(end);
end

end