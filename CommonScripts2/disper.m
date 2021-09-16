function k = disper(nyu, sigma, rho, g, h)
% nyu = 3.93;

if ~exist('sigma', 'var') sigma = 72.3;end  % surface tension 
if ~exist('rho', 'var') rho = 1.0; end % fluid mass density
if ~exist('g', 'var') g = 981.9; end % gravity acceleration
kk = 0; 
if ~exist('h', 'var')
    
    for nyu = nyu
        kk = kk + 1;
        kt = roots([sigma/rho 0 g -(nyu*2*pi)^2]);
        k(kk) = kt(3);
    end
else
    syms xk
    k = double(vpasolve(nyu == (((sigma/rho*xk^3+g*xk)*tanh(xk*h)).^.5)/2/pi, xk, [0 100]));
    
end