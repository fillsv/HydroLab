clear x y u v
Lx = 66;
Ly = 66;
nx = 130;
ny = nx
dx = Lx/nx;
dy = Ly/ny;
freq = 12;
dt = 1/freq;
nyu = 3;
omega = nyu*2*pi;
k = disper(nyu);
t = dt:dt:10;
[px py] = meshgrid( -(Lx/2-dx/2):dx:(Lx/2-dx/2), -(Ly/2-dy/2):dy:(Ly/2-dy/2));
% size(px)
kk = 0;
k0 = k;
% k = k*ones(size(px))+3e-2*(rand(size(px))-.5);
if exist('Err', 'var') == 0
    Err = 10;
end
d = 9.8+1.49/1.38;
n = 1.38;
vt = zeros(1, numel(t));
Ax = 1;
Ay = 1;
Ao = 2;
for tt = t
    kk = kk + 1;
%     vx = Ax*sin(k.*px)*sin(omega*tt);
    vx = Ax/2*sin(k.*px-omega*tt);
    vx = vx+Ax/2*sin(-k.*px-omega*tt);
    vx = vx+Ao/2/k0*cos(k.*py).*sin(k.*px);
    vx = vx+1e-3*(rand(size(px))-.5);
    vy = Ay/2*sin(k.*py-omega*tt);
    vy = vy+Ay/2*sin(-k.*py-omega*tt);
    vy = vy-Ao/2/k0*sin(k.*py).*cos(k.*px);
    vy = vy+1e-3*(rand(size(py))-.5);
%     vy = 1*sin(k.*py)*sin(omega*tt+.9*pi/2)+1e-4*(rand(size(px))-.5);
%     vx = vx.^1.05;
%     vt(kk) = sin(omega*tt);
%      vx = 1/2/k0*cos(k.*py).*sin(k.*px) + 2*sin(k.*px)*sin(omega*tt)+1e1*(rand(size(px))-.5);
%      vy = -1/2/k0*sin(k.*py).*cos(k.*px) + sin(k.*py)*sin(omega*tt)+1e1*(rand(size(px))-.5);

    x{kk,1} = px;
    y{kk,1} = py;
%     u1{kk,1} = vx*dt*(k*d*(n-1)/n);
%     v1{kk,1} = vy*dt*(k*d*(n-1)/n);
    u{kk,1} = vx*dt;
    v{kk,1} = vy*dt;
    
end
pathout= 'testsignal/mat/01/';
if ~exist(pathout, 'dir')
    mkdir(pathout)
end
save([pathout 'MVI_01.mat'], 'x', 'y', 'u', 'v');

% padFactor = 4;
% [Vx Vy] = findWaveVelocity(x1, y1, u1, v1, 1:48, 1, nyu, 'baseForm.mat', 1);
% VxVy=Vx*Vy
% 
% [rezx rezy rez_kx rez_ky rezFreq] = showFFTWaves_pars(x1, y1, u1, v1,...
%           1:48, nyu, padFactor, 'baseForm.mat',...
%           'vortFFTProfile.mat',  1);
% %         rez = rez./(k*d*(n-1)/n);
% Vx=onePeak(rezx, rez_kx, rez_ky, k, 0, padFactor) + onePeak(rezx, rez_kx, rez_ky, -k, 0, padFactor);
% Vy=onePeak(rezy, rez_kx, rez_ky, 0, k, padFactor) + onePeak(rezy, rez_kx, rez_ky, 0, -k, padFactor);
%         VxVy = Vx*Vy
% % return
% [rez rez_kx rez_ky rezFreq] = showFFTVxVy_pars(x1, y1, u1, v1,...
%           1:48, nyu, padFactor, 'baseForm.mat',...
%           'vortFFTProfile.mat',  1);
% %         rez = rez./(k*d*(n-1)/n);
%         Vxy=sumPeak(rez, rez_kx, rez_ky, k, padFactor);
%         
%         Vxy
%         
%         
% %  [omega kx ky fomegaAll] = showFFTVorticity_pars(x1, y1, u1, v1,...
% %           1:48, padFactor, 'baseForm.mat',...
% %           'vortFFTProfile.mat');
% % Omega=sumPeak(omega, kx, ky, k, padFactor);
% % Omega
 
% function [A xmax1 ymax1] = onePeak(fo, kx, ky, kfindx, kfindy, padFactor)
%     fo = abs(fo);
%     [fmax1 xmax1 ymax1] = findPeak(fo, kx, ky, kfindx,  kfindy, .06, .06);
%     ind1 = find(abs((kx-xmax1)+i*(ky-ymax1))<.2);
%     numelo = numel(fo)/padFactor^2;
%     A = (sum(fo([ind1;]).^2)/numel(fo)/numelo).^.5;
% end        
% function [A xmax1 ymax1] = sumPeak(fo, kx, ky, kfind, padFactor)
%     fo = abs(fo);
%     [fmax1 xmax1 ymax1] = findPeak(fo, kx, ky,  kfind,  kfind, .06, .06);
%     [fmax2 xmax2 ymax2] = findPeak(fo, kx, ky, -kfind,  kfind, .06, .06);
%     [fmax3 xmax3 ymax3] = findPeak(fo, kx, ky,  kfind, -kfind, .06, .06);
%     [fmax4 xmax4 ymax4] = findPeak(fo, kx, ky, -kfind, -kfind, .06, .06);
%     ind1 = find(abs((kx-xmax1)+i*(ky-ymax1))<.2);
%     ind2 = find(abs((kx-xmax2)+i*(ky-ymax2))<.2);
%     ind3 = find(abs((kx-xmax3)+i*(ky-ymax3))<.2);
%     ind4 = find(abs((kx-xmax4)+i*(ky-ymax4))<.2);
%     numelo = numel(fo)/padFactor^2;
%     A = 2*(sum(fo([ind1; ind2; ind3; ind4]).^2)/numel(fo)/numelo).^.5;
% end