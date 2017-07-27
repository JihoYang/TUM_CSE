function [x w v t] = legpts(n,int,meth)
%LEGPTS  Legendre points and Gauss-Legendre Quadrature Weights.
%  LEGPTS(N) returns N Legendre points X in (-1,1).
%
%  [X,W] = LEGPTS(N) returns also a row vector W of weights for
%  Gauss-Legendre quadrature.
%
%  LEGPTS(N,D) scales the nodes and weights for the domain D. D can be
%  either a vector with two components or a domain object. If the interval
%  is infinite, the map is chosen to be the default 'unbounded map' with
%  mappref('parinf') = [1 0] and mappref('adaptinf') = 0.
%
%  [X,W,V] = LEGPTS(N) returns additionally a column vector V of weights in
%  the barycentric formula corresponding to the points X. The weights are
%  scaled so that max(abs(V)) = 1.
%
%  [X,W] = LEGPTS(N,METHOD) allows the user to select which method to use.
%    METHOD = 'REC' uses the recurrence relation for the Legendre 
%       polynomials and their derivatives to perform Newton iteration 
%       on the WKB approximation to the roots. Default for N < 100.
%    METHOD = 'ASY' uses the Hale-Townsend fast algorithm based up
%       asymptotic formulae, which is fast and accurate. Default for 
%       N >= 100.
%    METHOD = 'GLR' uses the Glaser-Liu-Rokhlin fast algorithm [2], which
%       is fast and can give better relative accuracy for the -.5<x<.5
%       than 'ASY' (although the accuracy of the weights is usually worse).
%    METHOD = 'GW' will use the traditional Golub-Welsch eigenvalue method, 
%       which is maintained mostly for historical reasons.
%
%  See also chebpts, jacpts, legpoly.

%  Copyright 2011 by The University of Oxford and The Chebfun Developers. 
%  See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

%  'GW' by Nick Trefethen, March 2009 - algorithm adapted from [1].
%  'GLR' by Nick Hale, April 2009 - algorithm adapted from [2].
%  'REC' by Nick Hale, July 2011
%  'ASY' by Nick Hale & Alex Townsend, May 2012 - see [3].
%
%  References:
%   [1] G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature
%       rules", Math. Comp. 23:221-230, 1969, 
%   [2] A. Glaser, X. Liu and V. Rokhlin, "A fast algorithm for the 
%       calculation of the roots of special functions", SIAM Journal  
%       on Scientific Computing", 29(4):1420-1438:, 2007.
%   [3] N. Hale and A. Townsend, "Fast computation of Gauss-Jacobi 
%       quadrature nodes and weights",In preparation, 2012.

% Defaults
interval = [-1,1];
method = 'default';
method_set = nargin == 3;
t = [];

% Deal with trivial cases
if n < 0
    error('CHEBFUN:legpts:n', 'First input should be a positive number.');
elseif n == 0   % Return empty vector if n == 0
    x = []; w = []; v = []; return
elseif n == 1   % n == 1    
    x = 0;  w = 2;  v = 1; t = 1; return
elseif n == 2   % n == 2
    x = [-1 ; 1]/sqrt(3); w = [1 1]; v = [1 ; -1]; t = acos(x); return
end

% Check the inputs
if nargin > 1
    if nargin == 3
        interval = int; method = meth;
    elseif nargin == 2
        if ischar(int), method = int; method_set = true; else interval = int; end
    end
    if ~any(strcmpi(method,{'default','GW','fast','fastsmall','GLR','ASY','REC'}))
        error('CHEBFUN:legpts:inputs',['Unrecognised input string.', method]); 
    end
    if isa(interval,'domain')
        interval = interval.endsandbreaks;
    end
    if numel(interval) > 2,
        warning('CHEBFUN:legpts:domain',...
            'Piecewise intervals not supported and will be ignored.');
        interval = interval([1 end]);
    end
end

% Choose the method
if (n < 100 && ~method_set) || any(strcmpi(method,{'fastsmall','rec'}))
   [x w v] = rec(n);        % REC ('fastsmall' is for backward compatibiilty)
elseif strcmpi(method,'GW')
   [x w v] = gw(n);         % GW  see [1]
elseif strcmpi(method,'GLR')
   [x w v] = alg0_Leg(n);   % GLR see [2]
else
   [x w v t] = asy(n);        % HT  see [3]
end
v = abs(v); v = v./max(v); v(2:2:end) = -v(2:2:end);

if nargout == 4 && isempty(t), t = acos(x); end

% Rescale to arbitrary interval
if ~all(interval == [-1 1])
    if ~any(isinf(interval))    % Finite interval
        dab = diff(interval);
        x = (x+1)/2*dab + interval(1);
        w = dab*w/2;
    else                        % Infinite interval
        m = maps(fun,{'unbounded'},interval); % use default map
        if nargout > 1, w = w.*m.der(x.'); end
        x = m.for(x);
        x([1 end]) = interval([1 end]);
    end
end

end

%% -------------------- Routines for GW algorithm ------------------

function [x w v] = gw(n)
   beta = .5./sqrt(1-(2*(1:n-1)).^(-2)); % 3-term recurrence coeffs
   T = diag(beta,1) + diag(beta,-1);     % Jacobi matrix
   [V,D] = eig(T);                       % Eigenvalue decomposition
   x = diag(D); [x,i] = sort(x);         % Legendre points
   w = 2*V(1,i).^2;                      % Quadrature weights
   v = sqrt(1-x.^2).*abs(V(1,i))';       % Barycentric weights
   % Enforce symmetry
   ii = 1:floor(n/2);  x = x(ii);  w = w(ii); 
   vmid = v(floor(n/2)+1); v = v(ii);
   if mod(n,2)
        x = [x ; 0 ; -x(end:-1:1)];  w = [w  2-sum(2*w) w(end:-1:1)];
        v = [v ; vmid ; v(end:-1:1)];
   else
        x = [x ; -x(end:-1:1)];      w = [w w(end:-1:1)];      
        v = [v ; v(end:-1:1)];
   end
end

%% -------------------- Routines for REC algorithm ------------------

function [x, w, v] = rec(n)

    % Asymptotic formula (WKB) - only positive x.
    if mod(n,2), s = 1; else s = 0; end 
    k = (n+s)/2:-1:1; theta = pi*(4*k-1)/(4*n+2);
    x = (1-(n-1)/(8*n^3)-1/(384*n^4)*(39-28./sin(theta).^2)).*cos(theta);

    % Initialise
    Pm2 = 1; Pm1 = x;  PPm2 = 0; PPm1 = 1;
    dx = inf; l = 0;

    % Loop until convergence
    while norm(dx,inf) > eps && l < 10
        l = l + 1;
        for k = 1:n-1, 
            P = ((2*k+1)*Pm1.*x-k*Pm2)/(k+1);           Pm2 = Pm1; Pm1 = P; 
            PP = ((2*k+1)*(Pm2+x.*PPm1)-k*PPm2)/(k+1);  PPm2 = PPm1; PPm1 = PP;  
        end
        dx = -P./PP; x = x + dx;    
        Pm2 = 1; Pm1 = x; PPm2 = 0; PPm1 = 1;
    end

    % Once more for derivatives
    for k = 1:n-1, 
        P = ((2*k+1)*Pm1.*x-k*Pm2)/(k+1);           Pm2 = Pm1; Pm1 = P; 
        PP = ((2*k+1)*(Pm2+x.*PPm1)-k*PPm2)/(k+1);  PPm2 = PPm1; PPm1 = PP;  
    end    
%     PP = -n*(x.*P-Pm2)./(1-x.^2);

    % Reflect for negative values
    x = [-x(end:-1:1+s) x].';
    ders = [PP(end:-1:1+s) PP].';

    w = 2./((1-x.^2).*ders.^2)';          % Quadrature weights
    v = 1./ders;                          % Barycentric weights  

end

%% -------------------- Routines for GLR algorithm ------------------------

function [x, w, v] = alg0_Leg(n) % Driver for GLR

    % Compute coefficients of P_m(0), m = 0,..,N via recurrence relation.
    Pm2 = 0; Pm1 = 1; 
    for k = 0:n-1, P = -k*Pm2/(k+1); Pm2 = Pm1; Pm1 = P; end

    % Get the first roots and derivative values to initialise.
    x = zeros(n,1); ders = zeros(n,1);               % Allocate storage
    if mod(n,2)                                      % n is odd
        x((n-1)/2) = 0;                              % Zero is a root
        ders((n+1)/2) = n*Pm2;                       % P'(0)    
    else                                             % n is even
        [x(n/2+1), ders(n/2+1)] = alg2_Leg(P,n);      % Find first root
    end       
    [x, ders] = alg1_Leg(x,ders);          % Other roots and derivatives

    w = 2./((1-x.^2).*ders.^2)';          % Quadrature weights
    v = 1./ders;                          % Barycentric weights

end

% ---------------------

function [roots, ders] = alg1_Leg(roots,ders)  % Main algorithm for GLR

    n = length(roots);
    if mod(n,2), N = (n-1)/2; s = 1; else N = n/2; s = 0; end   

    % Approximate roots via asymptotic formula.
    k = (n-2+s)/2:-1:1; theta = pi*(4*k-1)/(4*n+2);
    roots(((n+4-s)/2):end) = (1-(n-1)/(8*n^3)-1/(384*n^4)*(39-28./sin(theta).^2)).*cos(theta);
    x = roots(N+1);

    % Number of terms in Taylor expansion.
    m = 30;

    % Storage
    hh1 = ones(m+1,1); zz = zeros(m,1); u = zeros(1,m+1); up = zeros(1,m+1);

    % Loop over all the roots we want to find (using symmetry).
    for j = N+1:n-1
        % Distance to initial approx for next root
        h = roots(j+1) - x;

        % Recurrence Taylor coefficients (scaled & incl factorial terms).
        M = 1/h;                           % Scaling
        c1 = 2*x/M; c2 = 1./(1-x^2);       % Some constants
        % Note, terms are flipped for more accuracy in inner product
        u([m+1 m]) = [0 ders(j)/M];  up(m+1) = u(m);
        for k = 0:m-2
            up(m-k) = (c1*(k+1)*u(m-k)+(k-n*(n+1)/(k+1))*u(m-k+1)/M^2)*c2;
            u(m-(k+1)) = up(m-k)/(k+2);
        end
        up(1) = 0;  

        % Newton iteration
        hh = hh1; step = inf;  l = 0; 
        while (abs(step) > eps) && (l < 10)
            l = l + 1;
            step = (u*hh)/(up*hh)/M;
            h = h - step;        
            Mhzz = (M*h)+zz;
            hh = [1;cumprod(Mhzz)];     % Powers of h (This is the fastest way!)
            hh = hh(end:-1:1);          % Flip for more accuracy in inner product 
        end

        % Update
        x = x + h;
        roots(j+1) = x;
        ders(j+1) = M*(up*hh);  

    end

    % Nodes are symmetric.
    roots(1:N+s) = -roots(n:-1:N+1);
    ders(1:N+s) = ders(n:-1:N+1);

end

% ---------------------

function [x1, d1] = alg2_Leg(Pn0,n) % Find the first root (note P_n'(0)==0)

    % Approximate first root via asymptotic formula
    k = ceil(n/2); theta = pi*(4*k-1)/(4*n+2);
    x1 = (1-(n-1)/(8*n^3)-1/(384*n^4)*(39-28./sin(theta).^2)).*cos(theta);

    m = 30; % Number of terms in Taylor expansion.

    % Recurrence Taylor coefficients (scaled & incl factorial terms).
    M = 1/x1; % Scaling
    zz = zeros(m,1); u = [Pn0 zeros(1,m)]; up = zeros(1,m+1); % Storage
    for k = 0:2:m-2
        up(k+2) = (k-n*(n+1)/(k+1))*u(k+1)/M^2;
        u(k+3) = up(k+2)/(k+2);
    end
    % Flip for more accuracy in inner product calculation.
    u = u(m+1:-1:1); up = up(m+1:-1:1);

    % Newton iteration
    x1k = ones(m+1,1); step = inf; l = 0;
    while (abs(step) > eps) && (l < 10)
        l = l + 1;
        step = (u*x1k)/(up*x1k)/M;
        x1 = x1 - step;
        x1k = [1;cumprod(M*x1+zz)]; % Powers of h (This is the fastest way!)
        x1k = x1k(end:-1:1);        % Flip for more accuracy in inner product
    end

    % Get the derivative at this root, i.e. P'(x1).
    d1 = M*(up*x1k);

end

%% -------------------- Routines for ASY algorithm ------------------------

function [x w v t] = asy(n)

    % Determine switch between interior and boundary regions
    nbdy = min(10,floor(n/2));

    % Interior
    [x w v t] = asy1(n,nbdy);   

    % Boundary
    [xbdy wbdy vbdy tbdy] = asy2(n,nbdy); 
    
    % Combine
    bdyidx1 = n-(nbdy-1):n; bdyidx2 = nbdy:-1:1;
    x(bdyidx1) = xbdy;  w(bdyidx1) = wbdy; v(bdyidx1) = vbdy; t(bdyidx1) = tbdy;
    x(bdyidx2) = -xbdy; w(bdyidx2) = wbdy; v(bdyidx2) = vbdy; t(bdyidx2) = -tbdy;
    
end

function [x w v t] = asy1(n,nbdy)
    % Interior method
  
    % Approximate roots via asymptotic formula. (Tricomi)
    s = mod(n,2);
    k = (n-2+s)/2+1:-1:1; theta = pi*(4*k-1)/(4*n+2);
    x = (1-(n-1)/(8*n^3)-1/(384*n^4)*(39-28./sin(theta).^2)).*cos(theta);
    t = acos(x);
    
    if n < 666
        % Approximation for Legendre roots (See Olver 1974)
        idx = (x>.5); npts = sum(idx);
        % Roots pf the Bessel function J_0 (Precomputed in Mathematica)
        jk = [2.404825557695773     5.520078110286310    8.653727912911012 ...
             11.791534439014281    14.930917708487785   18.071063967910922 ...
             21.211636629879258    24.352471530749302   27.493479132040254 ... 
             30.634606468431975    33.775820213573568].';
        if npts > 11
            % Esimate the larger Bessel roots (See Branders et al., JCP 1981).
            p = ((length(jk)+1:npts).'-.25)*pi;  pp = p.*p;
            num = 0.0682894897349453 + pp.*(0.131420807470708 + ...
                  pp.*(0.0245988241803681 + pp.*0.000813005721543268));
            den = p.*(1.0 + pp.*(1.16837242570470 + pp.*(0.200991122197811 + ...
                  pp.*(0.00650404577261471))));
            jk = [jk ; p + num./den];
        end
        phik = jk(1:npts)/(n+.5);
        tnew = phik + (phik.*cot(phik)-1)./(8*phik*(n+.5)^2);
        t(idx) = tnew(end:-1:1);
    end

    % locate the boundary node
    mint = t(end-nbdy+1);
    idx = max(find(t<mint,1)-1,1);

    dt = inf; j = 0;
    % Newton iteration
    while norm(dt,inf) > sqrt(eps)/1000
        [vals ders] = feval_asy1(n,t,mint,1); % Evaluate via asy formulae
        dt = vals./ders;                      % Newton update
        t = t - dt;                           % Next iterate
        j = j + 1;
        dt = dt(1:idx-1);
        if j > 10, dt = 0; end
    end
    [vals,ders] = feval_asy1(n,t,mint,1);  % once more for good ders.
    t = t - vals./ders;                    % Newton update

    x = cos(t);
    w = 2./ders.^2;
    v = sin(t)./ders;

    % Flip using symetry for negative nodes
    if s
        x = [-x(end:-1:2) x].';  w = [w(end:-1:2) w];  v = -[v(end:-1:2) v].'; t = [-t(end:-1:2) t].';
    else
        x = [-x(end:-1:1) x].';  w = [w(end:-1:1) w];  v = [-v(end:-1:1) v].';  t = [-t(end:-1:1) t].';
    end
    
end


function [vals ders] = feval_asy1(n,t,mint,flag)
    % Evaluate 1st asymptotic formula (interior)

    M = 20;         % Max number of expansion terms.
    % Asymptotic expansion.
    c = cumprod((1:2:2*M-1)./(2:2:2*M));
    d = cumprod((1:2:2*M-1)./(2*n+3:2:2*(n+M)+1));
    c = [1 c.*d];  	% Coefficients in expansion.
    % How many terms required in the expansion?
    R = (8/pi)*c./(2*sin(mint)).^(.5:M+1)/10;
    R = R(abs(R)>eps); M = length(R); c = c(1:M);

    % Constant out the front ( C = sqrt(4/pi)*gamma(n+1)/gamma(n+3/2) )
    ds = -1/8/n; s = ds; j = 1;
    while abs(ds/s) > eps/100
        j = j+1;
        ds = -.5*(j-1)/(j+1)/n*ds;
        s = s + ds;
    end
    p2 = exp(s)*sqrt(4/(n+.5)/pi);
    g = [1 1/12 1/288 -139/51840 -571/2488320 163879/209018880 ...
         5246819/75246796800 -534703531/902961561600 ...
         -4483131259/86684309913600 432261921612371/514904800886784000];
    f = @(z) sum(g.*[1 cumprod(ones(1,9)./z)]);
    C = p2*(f(n)/f(n+.5));

    % Some often used vectors/matrices
    onesT = ones(1,length(t));
    onesM = ones(M,1);
    M05 = transpose((0:M-1)+.5);
    onesMcotT = onesM*cot(t);
    M05onesT = M05*onesT;
    twoSinT = onesM*(2*sin(t));
    denom = cumprod(twoSinT)./sqrt(twoSinT);

%     alpha = onesM*(n*t) + M05onesT.*(onesM*(t-.5*pi));
%     cosAlpha = cos(alpha);
%     sinAlpha = sin(alpha);

    % Taylor expansion of cos(alpha0);
%     if flag
        k = numel(t):-1:1;
        rho = n+.5;
        ta = double(single(t));    tb = t - ta;
        hi = rho*ta;               lo = rho*tb;
        pia = double(single(pi));
        pib = -8.742278000372485e-08; %pib = pi - pia;
        dh = (hi-(k-.25)*pia)+lo-(k-.25)*pib;
        tmp = 0;
        sgn = 1; fact = 1; DH = dh; dh2 = dh.*dh;
        for j = 0:20
            dc = sgn*DH/fact;
            tmp = tmp + dc;
            sgn = -sgn;
            fact = fact*(2*j+3)*(2*j+2);
            DH = DH.*dh2;
            if norm(dc,inf) < eps/2000, break, end
        end
        tmp(2:2:end) = -tmp(2:2:end);
        tmp = sign(cos((n+.5)*t(2)-.25*pi)*tmp(2))*tmp;
        cosAlpha(1,:) = tmp;
        
        tmp = 0; sgn = 1; fact = 1; DH = 1; dh2 = dh.*dh;
        for j = 0:20
            dc = sgn*DH/fact;
            tmp = tmp + dc;
            sgn = -sgn;
            fact = fact*(2*j+2)*(2*j+1);
            DH = DH.*dh2;
            if norm(dc,inf) < eps/2000, break, end
        end
        tmp(2:2:end) = -tmp(2:2:end);
        tmp = sign(sin((n+.5)*t(2)-.25*pi)*tmp(2))*tmp;
        sinAlpha(1,:) = tmp;

        sint = sin(t); cost = cos(t);
        for k = 2:M
            cosAlpha(k,:) = cosAlpha(k-1,:).*sint+sinAlpha(k-1,:).*cost;
            sinAlpha(k,:) = sinAlpha(k-1,:).*sint-cosAlpha(k-1,:).*cost;
        end
%     end

    % Sum up all the terms.
    vals = C*(c*(cosAlpha./denom));             
    numer = M05onesT.*(cosAlpha.*onesMcotT + sinAlpha) + n*sinAlpha;
    ders = -C*(c*(numer./denom)); % (dP/dtheta)

end

function [x w v t] = asy2(n,npts)
    % Boundary method

    if npts > ceil((n+1)/2), error('CHEBFUN:legpts:asy2:N', ...
            'NPTS must be <= N/2'); end

    % Approximation for Legendre roots (See Olver 1974)
    % Roots pf the Bessel function J_0 (Precomputed in Mathematica)
    jk = [2.404825557695773     5.520078110286310    8.653727912911012 ...
         11.791534439014281    14.930917708487785   18.071063967910922 ...
         21.211636629879258    24.352471530749302   27.493479132040254 ... 
         30.634606468431975    33.775820213573568].';
    phik = jk(1:npts)/(n+.5);
    t = phik + (phik.*cot(phik)-1)./(8*phik*(n+.5)^2);

    [tB1 A2 tB2 A3] = asy2_higherterms(0,0,t,n);
    dt = inf; j = 0;
    % Newton iteration
    while norm(dt,inf) > sqrt(eps)/200
        [vals ders] = feval_asy2(n,t,0);   % Evaluate via asy formula
        dt = vals./ders;                   % Newton update
        t = t + dt;                        % Next iterate
        j = j + 1; if j > 10, dt = 0; end  % Bail
    end
    % Once more for good ders.
    [vals, ders] = feval_asy2(n,t,1);      %#ok<ASGLU> 

    % flip
    t = t(npts:-1:1); ders = ders(npts:-1:1);
    % Revert to x-space
    x = cos(t); w = (2./ders.^2).'; v = sin(t)./ders;

    function [vals, ders] = feval_asy2(n,t,flag)
        % Evaluate 2nd asymptotic formula (boundary)
        
        % Useful constants
        rho = n + .5; rho2 = n - .5;
        
        % Evaluate the Bessel functions
        Ja = besselj(0,rho*t,0);
        Jb = besselj(1,rho*t,0);
        Jbb = besselj(1,rho2*t,0);
        if ~flag
            Jab = besselj(0,rho2*t,0);
        else
            % In the final step, perform accurate evaluation
            Jab = besseltaylor(-t,rho*t);
        end
        
        % Evaluate functions for recurrsive definition of coefficients.
        gt = .5*(cot(t) - 1./t);
        gtdt = .5*(-csc(t).^2 + 1./t.^2);
        tB0 = .25*gt;
        A1 = gtdt/8 - 1/8*gt./t - gt.^2/32;
        tB1t = tB1(t); A2t = A2(t);    % Higher terms

        % VALS:
        vals = Ja + Jb.*tB0/rho + Ja.*A1/rho^2 + Jb.*tB1t/rho^3 + Ja.*A2t/rho^4;
        % DERS:
        vals2 = Jab + Jbb.*tB0/rho2 + Jab.*A1/rho2^2 + Jbb.*tB1t/rho2^3 + Jab.*A2t/rho2^4;
        
        % Higher terms (not needed for n > 1000).
        tB2t = tB2(t);
        A3t = A3(t);
        vals = vals + Jb.*tB2t/rho^5 + Ja.*A3t/rho^6;
        vals2 = vals2 + Jbb.*tB2t/rho2^5 + Jab.*A3t/rho2^6;
        
        % Relation for derivative
        ders = n*(-cos(t).*vals + vals2)./sin(t);
        
        % Common factors
        denom = sqrt(t./sin(t));
        ders = ders.*denom;
        vals = vals.*denom;  
       
    end

end

function Ja = besseltaylor(t,z)
    % Accurate evaluation of Bessel function for asy2
    npts = numel(t);
    kmax = min(ceil(abs(log(eps)/log(norm(t,inf)))),30);
    H = bsxfun(@power,t,0:kmax).';
    % Compute coeffs in Taylor expansions about z (See NIST 10.6.7)
    [nu, JK] = meshgrid(-kmax:kmax, z);
    Bjk = besselj(nu,JK,0);
    nck = abs(pascal(floor(1.25*kmax),1)); nck(1,:) = []; % nchoosek    
    AA = [Bjk(:,kmax+1) zeros(npts,kmax)];
    fact = 1;
    for k = 1:kmax
        sgn = 1;
        for l = 0:k
            AA(:,k+1) = AA(:,k+1) + sgn*nck(k,l+1)*Bjk(:,kmax+2*l-k+1);
            sgn = -sgn;
        end
        fact = k*fact;
        AA(:,k+1) = AA(:,k+1)/2^k/fact;
    end
    % Evaluate Taylor series
    Ja = zeros(npts,1);
    for k = 1:npts
        Ja(k,1) = AA(k,:)*H(:,k);
    end
end

function [tB1 A2 tB2 A3 tB3 A4] = asy2_higherterms(a,b,theta,n)
    % Compute the higher order terms in asy2 boundary formula 

    % The constants a = alpha and b = beta
    A = (.25-a^2); B = (.25-b^2); % These are more useful

    % For now, just work on half of the domain
    % c = pi/2; N = 30;
    c = max(max(theta),.5); 
    if n < 30, N = ceil(20-(n-20)); else N = 10; end
    if n > 30 && c > pi/2-.5, N = 15; end
    N1 = N-1;

    % 2nd-kind Chebyshev points and barycentric weights
    t = .5*c*(sin(pi*(-N1:2:N1)/(2*N1)).'+1);        
    v = [.5 ; ones(N1,1)]; v(2:2:end) = -1; v(end) = .5*v(end);

    % The g's
    g = A*(cot(t/2)-2./t)-B*tan(t/2);
    gp = A*(2./t.^2-.5*csc(t/2).^2)-.5*(.25-b^2)*sec(t/2).^2;
    gpp = A*(-4./t.^3+.25*sin(t).*csc(t/2).^4)-4*B*sin(t/2).^4.*csc(t).^3;
    g(1) = 0; gp(1) = -A/6-.5*B; gpp(1) = 0;

    % B0
    B0 = .25*g./t;
    B0p = .25*(gp./t-g./t.^2);
    B0(1) = .25*(-A/6-.5*B);
    B0p(1) = 0;

    % A1
    A10 = a*(A+3*B)/24;
    A1 = .125*gp - (1+2*a)/2*B0 - g.^2/32 - A10;
    A1p = .125*gpp - (1+2*a)/2*B0p - gp.*g/16;
    A1p_t = A1p./t;
    A1p_t(1) = -A/720-A^2/576-A*B/96-B^2/64-B/48+a*(A/720+B/48);

    % Make f accurately
    fcos = B./(2*cos(t/2)).^2;
    f = -A*(1/12+t.^2/240+t.^4/6048+t.^6/172800+t.^8/5322240 + ...
        691*t.^10/118879488000+t.^12/5748019200+3617*t.^14/711374856192000 + ...
        43867*t.^16/300534953951232000);
    idx = t>.5;
    ti = t(idx);
    f(idx) = A.*(1./ti.^2 - 1./(2*sin(ti/2)).^2);
    f = f - fcos;

    % Integrals for B1
    C = cumsummat(N)*(.5*c); 
    D = diffmat(N)*(2/c);
    I = (C*A1p_t);
    J = (C*(f.*A1));

    % B1
    tB1 = -.5*A1p - (.5+a)*I + .5*J;   
    
    tB1(1) = 0;
    B1 = tB1./t;
    B1(1) = A/720+A^2/576+A*B/96+B^2/64+B/48+a*(A^2/576+B^2/64+A*B/96)-a^2*(A/720+B/48);

    % A2
    K = C*(f.*tB1);
    A2 = .5*(D*tB1) - (.5+a)*B1 - .5*K;
    A2 = A2 - A2(1);

    if nargout < 3
        % Make function for output
        tB1 = @(theta) bary(theta,tB1,t,v);
        A2 = @(theta) bary(theta,A2,t,v);
    end

    % A2p
    A2p = D*A2;
    A2p = A2p - A2p(1);
    A2p_t = A2p./t;
    % Extrapolate point at t = 0
    w = pi/2-t(2:end);
    w(2:2:end) = -w(2:2:end);
    w(end) = .5*w(end);
    A2p_t(1) = sum(w.*A2p_t(2:end))/sum(w);

    % B2
    tB2 = -.5*A2p - (.5+a)*(C*A2p_t) + .5*C*(f.*A2);
    B2 = tB2./t;
    % Extrapolate point at t = 0
    B2(1) = sum(w.*B2(2:end))/sum(w);

    % A3
    K = C*(f.*tB2);
    A3 = .5*(D*tB2) - (.5+a)*B2 - .5*K;
    A3 = A3 - A3(1);

    if nargout < 6
        % Make function for output
        tB1 = @(theta) bary(theta,tB1,t,v);
        A2 = @(theta) bary(theta,A2,t,v);
        tB2 = @(theta) bary(theta,tB2,t,v);
        A3 = @(theta) bary(theta,A3,t,v);
        return
    end

    % A2p
    A3p = D*A3;
    A3p = A3p - A3p(1);
    A3p_t = A3p./t;
    % Extrapolate point at t = 0
    w = pi/2-t(2:end);
    w(2:2:end) = -w(2:2:end);
    A3p_t(1) = sum(w.*A3p_t(2:end))/sum(w);

    % B2
    tB3 = -.5*A3p - (.5+a)*(C*A3p_t) + .5*C*(f.*A3);
    B3 = tB3./t;
    % Extrapolate point at t = 0
    B3(1) = sum(w.*B3(2:end))/sum(w);

    % A3
    K = C*(f.*tB3);
    A4 = .5*(D*tB3) - (.5+a)*B3 - .5*K;
    A4 = A4 - A4(1);

    % Make function for output
    tB1 = @(theta) bary(theta,tB1,t,v);
    A2 = @(theta) bary(theta,A2,t,v);
    tB2 = @(theta) bary(theta,tB2,t,v);
    A3 = @(theta) bary(theta,A3,t,v);
    tB3 = @(theta) bary(theta,tB3,t,v);
    A4 = @(theta) bary(theta,A4,t,v);    

end

  
function Q = cumsummat(N)
% CUMSUMMAT  Chebyshev integration matrix.
% Q = CUMSUMMAT(N) is the matrix that maps function values at N Chebyshev
% points to values of the integral of the interpolating polynomial at
% those points, with the convention that the first value is zero. 

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

N = N-1;

persistent cache    % stores computed values for fast return
if isempty(cache), cache = {}; end    % first call

if length(cache) >= N && ~isempty(cache{N})
  Q = cache{N};
  return
else
  cache{N} = [];
end

% Matrix mapping coeffs -> values.
T = cp2cdm(N);

% Matrix mapping values -> coeffs.
Tinv = cd2cpm(N);

% Matrix mapping coeffs -> integral coeffs. Note that the highest order
% term is truncated. 
k = 1:N;
k2 = 2*(k-1);  k2(1) = 1;  % avoid divide by zero
B = diag(1./(2*k),-1) - diag(1./k2,1);
v = ones(N,1); v(2:2:end) = -1;
B(1,:) = sum( diag(v)*B(2:N+1,:), 1 );
B(:,1) = 2*B(:,1); 

Q = T*B*Tinv;               
cache{N} = Q;

end

function T = cp2cdm(N)
% Values of Cheb. polys at Cheb nodes, x(n)=-cos(pi*n/N).
theta = pi*(N:-1:0)'/N;
T = cos( theta*(0:N) );
end

function C = cd2cpm(N)
% Three steps: Double the data around the circle, apply the DFT matrix,
% and then take half the result with 0.5 factor at the ends.
theta = (pi/N)*(0:2*N-1)';
F = exp( -1i*theta*(0:2*N-1) );  % DFT matrix
rows = 1:N+1;  % output upper half only
% Impose symmetries on data and coeffs.
C = real( [ F(rows,N+1) F(rows,N:-1:2)+F(rows,N+2:2*N) F(rows,1) ] );
C = C/N;  C([1 N+1],:) = 0.5*C([1 N+1],:);
end
  
function [x w v] = chebpts(n,d,kind)
%CHEBPTS  Chebyshev points in [-1,1].
%   CHEBPTS(N) returns N Chebyshev points of the 2nd-kind in [-1,1].
%
%   CHEBPTS(N,D), where D is vector of length 2 and N is a scalar integer,
%   scales the nodes and weights for the interval [D(1) D(2)]. If the
%   interval is infinite, the map is chosen to be the default 'unbounded
%   map' with mappref('parinf') = [1 0] and mappref('adaptinf') = 0. If
%   length(D) > 2 and N a vector of length(D)-1, then CHEBPTS returns a
%   column vector of the stacked N(k) Chebyshev points on the subintervals
%   D(k:k+1). If length(N) is 1, then D is treated as [D(1) D(end)].
%
%   [X W] = CHEBPTS(N,D) returns also a row vector of the (scaled) weights
%   for Clenshaw-Curtis quadrature (computed using [1]). (For nodes and
%   weights of Gauss-Chebyshev quadrature, use [X W] = JACPTS(N,-.5,-.5,D))
%
%   [X W V] = CHEBPTS(N,D) returns, in addition to X and W, the barycentric
%   weights V corresponding to the Chebyshev points X.
%
%   [X W V] = CHEBPTS(F) returns the Chebyshev nodes and weights
%   corresponding to the domain and length of the chebfun F.
%
%   [X W V] = CHEBPTS(N,KIND) or CHEBPTS(N,D,KIND) returns Chebyshev points
%   and weights of the 1st-kind if KIND = 1 and 2nd-kind if KIND = 2
%   (default). (Note that if KIND is not supplied, chebpts will always
%   return 2nd-kind points, regardless of the value of 'chebkind' in
%   chebfunpref.).
%
%   See also legpts, jacpts, lagpts, and hermpts.

%   Copyright 2011 by The University of Oxford and The Chebfun Developers. 
%   See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

%   [1] Jörg Waldvogel, "Fast construction of the Fejér and Clenshaw-Curtis
%   quadrature rules", BIT Numerical Mathematics 46 (2006), pp 195-202.

% Intialise
x = []; w = []; v = []; scale = false;

% Shortcut for simple input (2nd-kind points on [-1 1]).
if nargin == 1 && isnumeric(n) && numel(n) == 1 && n > 0
    if n == 1, x = 0; w = 2; v = 1; return, end % Special case
    m = n-1;
    x = sin(pi*(-m:2:m)/(2*m)).';       % Chebyshev points 
    if nargout > 1                      % Quadrature weights            
        w = weights2(n);
    end
    if nargout > 2                      % Barycentric weights
        v = [.5 ; ones(n-1,1)]; v(2:2:end) = -1; v(end) = .5*v(end);
    end
    return
end

% Parse the inputs
if isa(n,'chebfun')
    if numel(n) > 1
        error('CHEBFUN:chebpts:quasi',...
            'chebpts does not support quasi-matrices.');
    end
    nn = zeros(n.nfuns,1);
    for k = 1:n.nfuns
        nn(k) = n.funs(k).n;
    end
    if nargin == 1, d = 2; end % Default to 2nd-kind
    if nargout == 1
        x = chebpts(nn,n.ends,d);
    elseif nargout == 2
        [x w] = chebpts(nn,n.ends,d);
    else
        [x w v] = chebpts(nn,n.ends,d);
    end        
    return
elseif nargin == 1
    d = [-1 1];
    kind = 2;
elseif nargin == 2
    if isa(d,'domain')
       scale = true;
       kind = 2;
    elseif length(d) == 1
       kind = d;
       d = [-1 1];
    else
       scale = true;
       kind = 2; 
    end
elseif nargin == 3
    scale = true; 
end
if isa(d,'domain')
    d = d.endsandbreaks;   
end
if isempty(d) || ~any(n)
    return % Return empty vector if domain is empty or n == 0
end
if numel(n) == 1
    d = d([1 end]);
end

% Deal with the piecewise case (where d has breakpoints and n is a vector).
numints = numel(d)-1; 
if numints > 1
    if numel(n) ~= numints
        error('CHEBFUN:chebpts:numints', ...
            'Vector N does not match domain D.'); 
    end
    csn = cumsum([0 ; n(:)]);
    x = zeros(csn(end),1);
    if nargout == 1
        for k = 1:numints
           idxk = csn(k)+1:csn(k+1);
           x(idxk) = chebpts(n(k),d(k:k+1),kind);
        end
    elseif nargout == 2
        w = zeros(1,csn(end));
        for k = 1:numints
           idxk = csn(k)+1:csn(k+1);
           [x(idxk) w(idxk)] = chebpts(n(k),d(k:k+1),kind);
        end
    else
        w = zeros(1,csn(end)); v = zeros(csn(end),1);
        for k = 1:numints
            idxk = csn(k)+1:csn(k+1);
            [x(idxk) w(idxk) v(idxk)] = chebpts(n(k),d(k:k+1),kind);
        end
    end
    return
end    

if numel(n) > 1, 
    error('CHEBFUN:chebpts:vecn','Vector N does not match domain D.');
end

% Avoid unnecessary scaling
if (d(1)==-1 && d(2)==1), scale = false; end

% Allow strings to determine which kind of points
if ischar(kind)
    if      strcmpi(kind,'1st'), kind = 1;
    elseif  strcmpi(kind,'2nd'), kind = 2; end
end

if n < 0, 
    error('CHEBFUN:chebpts:posinpt',...
        'Input should be a nonnegative number');
elseif n == 1,
    x = 0; w = 2; v = 1;
else
    m = n-1;
    if kind == 1
        x = sin(pi*(-m:2:m)/(2*m+2)).';      % 1st-kind Chebyshev points
        if nargout > 1  % Quadrature weights
            w = weights1(n);
        end
        if nargout > 2  % Barycentric weights
            v = sin((2*(0:n-1)+1)*pi/(2*n)).'; v(2:2:end) = -v(2:2:end);
            if ~mod(n,2), v = v./max(abs(v)); end
        end
    else
        x = sin(pi*(-m:2:m)/(2*m)).';        % 2nd-kind Chebyshev points
        if nargout > 1  % Quadrature weights            
            w = weights2(n);
        end
        if nargout > 2  % Barycentric weights
            v = [.5 ; ones(n-1,1)]; v(2:2:end) = -1; v(end) = .5*v(end);
        end
    end
end

% Rescale if d is provided:
if scale   
    if ~any(isinf(d))   % Finite interval
        dab05 = .5*diff(d);
        x = x*dab05 + (d(1) + dab05);
        if ( kind == 2 ) 
            x([1,end]) = d([1,end]);
        end
        w = dab05*w;
    else                % Infinite interval
        m = maps(fun,{'unbounded'},d); % Use default map
        if nargout > 1  % Quadrature weights
            w = w.*m.der(x.');
            if isinf(d(1)), w(1) = 0; end
            if isinf(d(end)), w(end) = 0; end
        end
        x = m.for(x);
        if kind == 2    % Force endpoints for 2nd-kind points
            x([1 end]) = d([1 end]);
        end
    end        
end

function w = weights1(n) % 1st-kind Chebyshev weights
% Jörg Waldvogel, "Fast construction of the Fejér and Clenshaw-Curtis
% quadrature rules", BIT Numerical Mathematics 43 (1), p. 001-018 (2004).
% http://www2.maths.ox.ac.uk/chebfun/and_beyond/programme/slides/wald.pdf
if n == 1
    w = 2;
else
    % new
    L = 0:n-1; r = 2./(1-4*min(L,n-L).^2); s1 = sign(n/2-L); % Aux vecs
    w = real(ifft(s1.*r.*exp(1i*pi/n*L)));  % Fejer weights
    
%     % old
%     l = floor(n/2)+1;
%     K = 0:n-l;   
%     v = [2*exp(1i*pi*K/n)./(1-4*K.^2)  zeros(1,l)];
%     w = real(ifft(v(1:n) + conj(v(n+1:-1:2))));
end
end

function w = weights2(n) % 2nd-kind Chebyshev wieghts
% Jörg Waldvogel, "Fast construction of the Fejér and Clenshaw-Curtis 
% quadrature rules", BIT Numerical Mathematics 43 (1), p. 001-018 (2004).
% http://www2.maths.ox.ac.uk/chebfun/and_beyond/programme/slides/wald.pdf
if n == 1
    w = 2;
else
    % new
    n = n-1;
    u0 = 1/(n^2-1+mod(n,2));                      % Boundary weights
    L = 0:n-1; r = 2./(1-4*min(L,n-L).^2);        % Auxiliary vectors
    w = [ifft(r-u0) u0];                          % C-C weights
    
%     % old
%     m = n-1;  
%     c = zeros(1,n);
%     c(1:2:n) = 2./[1 1-(2:2:m).^2 ]; 
%     f = real(ifft([c(1:n) c(m:-1:2)]));
%     w = [f(1) 2*f(2:m) f(n)];  
end
end
end
  
  
function fx = bary(x,gvals,xk,ek)
% BARY  Barycentric interpolation with arbitrary weights/nodes.
%  P = BARY(X,GVALS,XK,EK) interpolates the values GVALS at nodes
%  XK in the point X using the barycentric weights EK.
%
%  P = BARY(X,GVALS) assumes Chebyshev nodes and weights.
%
%  All inputs should be column vectors.

%  Copyright 2011 by The University of Oxford and The Chebfun Developers.
%  See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

warning('off', 'MATLAB:divideByZero'); % TODO: Delete this.

% Parse inputs
n = length(gvals);
bary1flag = 0;  % If true, possibly use type-1 barycentric formula
if n == 1,                % The function is a constant
    fx = gvals*ones(size(x));
    return
end
if any(isnan(gvals)),     % The function is NaN
    fx = NaN(size(x));
    return
end
if nargin < 3,            % Default to Chebyshev nodes
    bary1flag = 1;
    xk = chebpts(n);
end
if nargin < 4,            % Default to Chebyshev weights
    ek = [.5 ; ones(n-1,1)];
    ek(2:2:end) = -1;
    ek(end) = .5*ek(end);
end
if bary1flag,             % Call a barycentric formula of type 1 or 2
    ind1 = find(imag(x) | x < -1 | x > 1);
    if ~isempty(ind1),
        fx = NaN*x;
        fx(ind1) = kind1(x(ind1),gvals,xk,ek);
        ind2 = find(isnan(fx));
        fx(ind2) = kind2(x(ind2),gvals,xk,ek);
    else
        fx = kind2(x,gvals,xk,ek);
    end
else
    fx = kind2(x,gvals,xk,ek);
end
% Try to clean up NaNs
for i = find(isnan(fx(:)))',
    indx = find(x(i)==xk,1);
    if ~isempty(indx),
        fx(i) = gvals(indx);
    end
end
end



function fx = kind2(x,gvals,xk,ek)
% Evaluate the second-kind barycentric formula. Typically
% this is the standard for evaluating a barycentric interpolant
% on the interval.
if numel(x) < length(xk), % Loop over evaluation points
    fx = zeros(size(x));  % Initialise return value
    for i = 1:numel(x),
        xx = ek ./ (x(i)-xk);
        fx(i) = (xx.'*gvals) / sum(xx);
    end
else                      % Loop over barycentric nodes
    num = zeros(size(x)); denom = num; % initialise
    for i = 1:length(xk),
        y = ek(i) ./ (x-xk(i));
        num = num + (gvals(i)*y);
        denom = denom + y;
    end
    fx = num ./ denom;
end
end



function fx = kind1(x,gvals,xk,ek)
% Evaluate the first-kind barycentric formula. Typically we
% use this formula for evaluating a polynomial outside the interval.
% If the number of nodes is >=600, we compute the log of the
% nodal polynomial in order to avoid under-/ overflow.
% This method is only called with Chebyshev nodes xk on [-1,1]!
n = length(xk);
scale = 2; x = scale*x; xk = scale*xk;
fx = zeros(size(x)); % Initialise return value
if numel(x) < n,     % Loop over evaluation points
    for i = 1:numel(x),
        fx(i) = (ek./(x(i)-xk)).' * gvals;
    end
else                 % Loop over interpolation nodes
    for i = 1:n,
        y = ek(i) ./ (x-xk(i));
        fx = fx + gvals(i)*y;
    end
end
% Evaluate nodal polynomial ell
if n < 600,
    ell = ones(size(x));
    if numel(x) < n, % Loop over evaluation points
        for i = 1:numel(x),
            ell(i) = prod(x(i)-xk);
        end
    else             % Loop over interpolation nodes
        for i = 1:n,
            ell = ell .* (x-xk(i));
        end
    end
else
    ell = zeros(size(x));
    if numel(x) < n, % Loop over evaluation points
        for i = 1:numel(x),
            ell(i) = sum(log(x(i)-xk));
        end
    else             % Loop over interpolation nodes
        for i = 1:n,
            ell = ell + log(x-xk(i));
        end
    end
    ell = exp(ell);
    if isreal(x) && isreal(gvals) && isreal(xk) && isreal(ek),
        ell = real(ell);
    end
end
fx = fx .* ell * (1/(scale*(1-n))*(-2/scale)^(n-2));
end
  
  
function D = diffmat(N,k)  
% DIFFMAT  Chebyshev differentiation matrix
% D = DIFFMAT(N) is the matrix that maps function values at N Chebyshev
% points to values of the derivative of the interpolating polynomial at
% those points. 
%
% D = DIFFMAT(N,K) is the same, but for the Kth derivative.
%
% The matrices are computed using the 'hybrid' formula of Schneider & 
% Werner [1] and Welfert [2] proposed by Tee [3].

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% References:
%  [1] Schneider, C. and Werner, W., "Some new aspects of rational 
%   interpolation", Math. Comp. (47) 285--299, 1986.
%  [2] Welfert, B. D., "Generation of pseudospectral matrices I", SINUM,
%   (34) 1640--1657.
%  [3] Tee, T. W., "An adaptive rational spectral method for differential
%   equations with rapidly varying solutions", Oxford DPhil Thesis, 2006.

persistent cache    % stores computed values for fast return
if isempty(cache), cache = {}; end    % first call

if nargin < 2, k = 1; end

if N == 0, D = []; return, end
if N == 1, D = 0; return, end

if length(cache) >= N && length(cache{N}) >= k && ~isempty(cache{N}{k})
  D = cache{N}{k};
  return
else
  cache{N}{k} = [];
end

% construct Chebyshev grid and weights
x = chebpts(N);
w = [.5 ; ones(N-1,1)]; w(2:2:end) = -1; w(N) = .5*w(N);

ii = (1:N+1:N^2)';              % indices of diagonal
Dx = bsxfun(@minus,x,x');       % all pairwise differences
Dx(ii) = Dx(ii) + 1;            % add identity
Dxi = 1./Dx;                    % reciprocal 
Dw = bsxfun(@rdivide,w.',w);    % pairwise divisions
Dw(ii) = Dw(ii) - 1;            % subtract identity

% k = 1
if ~isempty(cache{N}{1})
    D = cache{N}{1};                            % recover from cache
else
    D = Dw .* Dxi;
    D(ii) = 0; D(ii) = - sum(D,2);              % negative sum trick
    cache{N}{1} = D;                            % store in cache    
end

if k == 1, return, end

% k = 2
if k > 1 && ~isempty(cache{N}{2})
    D = cache{N}{2};                            % recover from cache
elseif k > 1
    D = 2*D .* (repmat(D(ii),1,N) - Dxi);
    D(ii) = 0; D(ii) = - sum(D,2);              % negative sum trick
    cache{N}{2} = D;                            % store in cache    
end

% higher orders
for n = 3:k
    if ~isempty(cache{N}{n})
        D = cache{N}{n};
    else
        D = n*Dxi .* (Dw.*repmat(D(ii),1,N) - D);
        D(ii) = 0; D(ii) = - sum(D,2);          % negative sum trick
        cache{N}{n} = D;                        % store in cache    
    end
end

if N < 2^11+2
  siz = whos('cache'); 
  if siz.bytes > cheboppref('maxstorage')
    cache = {};
  end
  cache{N}{k} = D;
end  

end



