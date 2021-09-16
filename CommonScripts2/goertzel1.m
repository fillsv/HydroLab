function S_dsplib = goertzel1(s, k)
%
%                W - z^-1
% H(z) = ---------------------------
%          1 - alpha * z^-1 + z^-2
%
%
%	W = W_N^(-k) = exp(2i*pi*k/N).
%	alpha = 2*cos(2*pi*k/N).
%
if size(s,1) == 1
%     s = s';
    N = size(s,2);
else
    N = size(s,1);
end

k = k - 1;
% s = s';
% size(s)
% N = size(s,1);
w = exp(2i*pi*k/N);

%alpha
alpha = 2*cos(-2*pi*k/N);

% Выход БИХ-фильтра с DSPLIB.ORG 
b_dsplib = [w, -1, 0];
a_dsplib = [1, -alpha, 1];
X = filter(b_dsplib, a_dsplib, s);
% size(X)
if size(s,1) == 1
%     s = s';
    S_dsplib = X(:,end)*exp(-1i*k*2*pi);
else
    S_dsplib = X(end,:,:,:)*exp(-1i*k*2*pi);
end

% S_dsplib = permute(X(end, :)*exp(-1i*k*2*pi), [2,1]);