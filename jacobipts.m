function [x, w, v] = jacobipts(n, a, b, int, meth)
%JACOBIPTS  Give zero points and weights for quadrature of gauss type.
%  
%  [x, w] = JACOBIPTS(n, a, b) give the zero points of p_n and weights
%  for quadrature of gauss type. The weight function is 
%  w =  (1-x)^a (1+x)^b, where a > -1, b >-1. 
%  Note that the weights in this routine is NOT THE SAME as that of gjacobi. 
% 
% See also gjacobi



% $Date:2017/04/25$
% Copyright (c) G.Wang. 
% Email: wangguanjie0@126.com


%  This program use some subroutines, and the copyright of those
%  sub-routines in this program belongs to Oxford and Chebfun Developedrs.

%  Copyright 2017 by The University of Oxford and The Chebfun Developers. See
%  http://www.chebfun.org/ for Chebfun information

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This program use the Jacobin matrix and matlab function eig to compute
%  the zero points and weigts. The cost is about O(n^3) due to eig,  
%  and can be improved to O(n^2). There are methods at cost about O(n), 
%  see
%  [1]: N. Hale and A. Townsend, "Fast and accurate computation of 
%       Gauss--Legendre and Gauss--Jacobi quadrature nodes and weights", 
%       SIAM Journal on Scientific Computing, 2013.
%  [2]: chebfun/jacpts, www.chebfun.org.  

%  'GW'  by Nick Trefethen, March 2009 - algorithm adapted from
%  [3]:  G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature
%  'REC' by Nick Hale, July 2011

%  see also 
%  [4]: John A. Gubner, "Gaussian Quadrature and the Eigenvalue Problem"
%  [5]: Dongbin Xiu, "Numerical methods for stochastic computations", 2010
%  [6]: Wikipedia, https://en.wikipedia.org/wiki/Jacobi_polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

interval = [-1, 1]; method = 'REC';
validStrings = {'GW', 'REC'};

% parse input
if ( nargin > 3 )
    if ( nargin == 5 && any(strcmpi(meth, validStrings)) )
        interval = int;
        method = meth;
    else
        if ( ischar(int) && any(strcmpi(int, validStrings)) )
            method = int;
        else
            interval = int;
        end
    end
end


if  ( (a <= -1) || (b <= -1) )
    error('alphabeta = [a, b] must satisfying a > -1, b >-1.');
elseif ( max(a, b) > 5 )
  warning('Orthogonal:jacpts',...
      'MAX(ALPHA, BETA) > 5. Results may not be accurate');
end

% Deal with trivial cases:
if ( n < 0 )
    error('CHEBFUN:jacpts:n', 'First input should be a positive number.');
elseif ( n == 0 )   % Return empty vectors if n == 0:
    x = []; 
    w = []; 
    v = []; 
    return
elseif ( n == 1 )
    x0 = (b-a)/(a+b+2);
    x = diff(interval)/2 * (x0+1) + interval(1); % map from [-1,1] to interval. 
    w = 2^(a+b+1)*beta(a+1, b+1) * diff(interval)/2;
    v = 1;
    return
end


if ( strcmpi(method, 'GW') )
    [x, w, v] = gw(n, a, b);
elseif ( strcmpi(method, 'REC') )
    [x, w, v] = rec(n, a, b);
end

% adjust quadrature weights on interval [-1,1]; 
C = 2^(a+b+1)*beta(a+1, b+1);
w = w/sum(w); w = w*C;


% Scale the nodes and quadrature weights:
[x, w] = rescale(x, w, interval, a, b);

% Scale the barycentric weights:
v = abs(v); 
v(2:2:end) = -v(2:2:end);
v = v./max(abs(v)); 


end

%% ========================= Sub- routines ===========================

function [x, w] = rescale(x, w, interval, a, b)
%RESCALE   Rescale nodes and weights to an arbitrary finite interval.
    if ( ~all(interval == [-1, 1]) )
        c1 = .5*sum(interval); 
        c2 = .5*diff(interval);
        w = c2^(a+b+1)*w;
        x = c1 + c2*x;    
    end
end

%% ------------------------- Routines for GW ----------------------------
    
function [x, w, v] = gw(n, a, b)
    ab = a + b;
    ii = (2:n-1)';
    abi = 2*ii + ab;
    aa = [(b - a)/(2 + ab)
          (b^2 - a^2)./((abi - 2).*abi)
          (b^2 - a^2)./((2*n - 2+ab).*(2*n+ab))];
    bb = [2*sqrt( (1 + a)*(1 + b)/(ab + 3))/(ab + 2) ; 
          2*sqrt(ii.*(ii + a).*(ii + b).*(ii + ab)./(abi.^2 - 1))./abi];
    TT = diag(bb,1) + diag(aa) + diag(bb,-1); % Jacobi matrix.
    [V, x] = eig( TT );                       % Eigenvalue decomposition.
    x = diag(x);                              % Jacobi points.
    % Quadrature weights:
    w = V(1,:).^2; 
    v = sqrt(1-x.^2).*abs(V(1,:))';           % Barycentric weights.
end

%% ------------------------- Routines for REC ---------------------------


function [x, w, v] = rec(n, a, b)
%REC   Compute nodes and weights using recurrrence relation.

   [x1, ders1] = rec_main(n, a, b, 1); % Nodes and P_n'(x)
   [x2, ders2] = rec_main(n, b, a, 0); % Nodes and P_n'(x)
   x = [-x2(end:-1:1) ; x1];
   ders = [ders2(end:-1:1) ; ders1];
   w = 1./((1-x.^2).*ders.^2)';        % Quadrature weights
   v = 1./ders;                        % Barycentric weights
end

function [x, PP] = rec_main(n, a, b, flag)
%REC_MAIN   Jacobi polynomial recurrence relation.

% Asymptotic formula (WKB) - only positive x.
if ( flag )
    r = ceil(n/2):-1:1;
else
    r = floor(n/2):-1:1;  
end
C = (2*r+a-.5)*pi/(2*n+a+b+1);
T = C + 1/(2*n+a+b+1)^2 * ((.25-a^2)*cot(.5*C) - (.25-b^2)*tan(.5*C));
x = cos(T).';

% Initialise:
dx = inf; 
l = 0;
% Loop until convergence:
while ( (norm(dx,inf) > sqrt(eps)/1000) && (l < 10) )
    l = l + 1;
    [P, PP] = eval_Jac(x, n, a, b);
    dx = -P./PP; 
    x = x + dx;
end
% Once more for derivatives:
[ignored, PP] = eval_Jac(x, n, a, b);

end

function [P, Pp] = eval_Jac(x, n, a, b)
%EVALJAC   Evaluate Jacobi polynomial and derivative via recurrence relation.

% Initialise:
ab = a + b;
P = .5*(a-b+(ab+2)*x);  
Pm1 = 1; 
Pp = .5*(ab+2);         
Ppm1 = 0; 

% n = 0 case:
if ( n == 0 )
    P = Pm1; 
    Pp = Ppm1; 
end

for k = 1:n-1
    % Useful values:
    A = 2*(k + 1)*(k + ab + 1)*(2*k + ab);
    B = (2*k + ab + 1)*(a^2 - b^2);
    C = prod(2*k + ab + (0:2)');
    D = 2*(k + a)*(k + b)*(2*k + ab + 2);

    % Recurrence:
    Pa1 = ( (B+C*x).*P - D*Pm1 ) / A;
    Ppa1 = ( (B+C*x).*Pp + C*P - D*Ppm1 ) / A;

    % Update:
    Pm1 = P; 
    P = Pa1;  
    Ppm1 =  Pp; 
    Pp = Ppa1;
end
end
