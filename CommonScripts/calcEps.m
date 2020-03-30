function eps = calcEps(tau, nyu)
if size(tau)~=size(nyu) error('size of nyu mast be equil size of tau'); end
for ii = 1:numel(tau)
    omega = nyu(ii)*2*pi;
    nu = 0.01;
    if nyu(ii)<2
        k = disper(nyu(ii), 70, 1, 981, 10);
    else
        k = disper(nyu(ii));    
    end
    L = 70;
    h = 2;
    gamma = sqrt(nu*k^2/omega);
    syms ep;
    var = 2*gamma^2*(1+1/gamma/k*(1/L/sqrt(2)))+gamma/2/sqrt(2)*ep^2/(ep^2-ep*sqrt(2)+1);
%     var = 2*gamma^2*(1+1/gamma/k*((2*h+L)/2/h/L/sqrt(2)))+gamma/2/sqrt(2)*ep^2/(ep^2-ep*sqrt(2)+1);
    ep = vpasolve(var==1/tau(ii)/omega);

    eps(ii) = double(ep(2));

%     vrel = 1/sqrt(ep^2-ep*sqrt(2)+1)
end

