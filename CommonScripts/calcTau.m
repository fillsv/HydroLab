function tau = calcTau(eps, nyu)
% if size(tau)~=size(nyu) error('size of nyu mast be equil size of tau'); end
for ii = 1:numel(nyu)
    omega = nyu(ii)*2*pi;
    nu = 0.01;
    k = disper(nyu(ii));
    L = 70;
    gamma = sqrt(nu*k^2/omega);
%     eps = 1
    var = 2*gamma^2*(1+1/gamma/k/L/sqrt(2))+gamma/2/sqrt(2)*eps^2/(eps^2-eps*sqrt(2)+1);
    tau(ii) = 1/var/omega;

%     eps(ii) = double(eps(2));

%     vrel = 1/sqrt(ep^2-ep*sqrt(2)+1)
end

