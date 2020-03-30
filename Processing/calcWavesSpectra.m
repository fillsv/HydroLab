function frameWavesSpectra = calcWavesSpectra(frame, padFactor, num)
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end    

    if(exist('padFactor', 'var')==0) 
        padFactor=4; 
    end

    px = frame.px{num(1)};
    py = frame.py{num(1)};    
    freq = frame.freq;
    [npy, npx] = size(px);

    sigma = 73;
    rho = 1;
    g = 981.9;
    h = 10;
    
    [kx ky dkx dky] = genkxky(px, py, padFactor);

    frameWavesSpectra.px = px;
    frameWavesSpectra.py = py;    
    frameWavesSpectra.kx = kx;
    frameWavesSpectra.ky = ky;
    frameWavesSpectra.Lx = frame.Lx;    
    frameWavesSpectra.Ly = frame.Ly;    
    frameWavesSpectra.freq = frame.freq;    
    frameWavesSpectra.padFactor = padFactor;
    frameWavesSpectra.tt(1:numel(num)) = frame.tt(num);
    
    vx = cat(3, frame.vx{num});
    vy = cat(3, frame.vy{num});

    vx = permute(vx, [3 1 2]);
    vy = permute(vy, [3 1 2]);
    
%     w = repmat(hann(size(vx,1)), [1 size(vx,2) size(vx,3)]);
    
    vx(find(isnan(vx))) = 0;
    vy(find(isnan(vy))) = 0;
    
%     pp = 1 + nyuCrop./frame.freq*numel(num);

%     k = disper(nyuCrop, sigma, rho, g, h);

%     w23 = hann(size(vx, 2))*hann(size(vx,3))'; 
%     w1 = hann(size(vx,1));
%     w = repmat(w1, [1, size(vx,2), size(vx,3)]);
    d = 9.4+1.48/1.33;
    n = 1.33;
%     padFac 
    padFactorTime = 4;
    dnu = 1/numel(num)*freq/padFactorTime;
    nyu = (0:numel(num)*padFactorTime/2-1)*dnu;

    k = [0 disper(nyu(2:end))];
    k = k';

%     max(nyu)
    w23 = hann(size(vx, 2))*hann(size(vx,3))'; 
    w1 = hann(size(vx,1));
    w = permute(w23.*permute(w1,[3,2,1]),[3 1 2]);
%     size(w)
%     size(nu)
    fvx = fftn(vx.*w, [padFactorTime*size(vx,1), padFactor*size(vx,2), padFactor*size(vx,3)]);
    fvx = fvx/sqrt(mean(w23(:).^2))/mean(w1);
    fvx = abs(fftshift(fftshift(fvx, 2), 3));
    fvy = fftn(vy.*w, [padFactorTime*size(vy,1), padFactor*size(vy,2), padFactor*size(vy,3)]);
    fvy = fvy/sqrt(mean(w23(:).^2))/mean(w1);
    fvy = abs(fftshift(fftshift(fvy, 2), 3));

    
    fvx = 2*fvx; % fft(t)
    fvy = 2*fvy;
    fv = abs(fvx+i*fvy);
    dk = dkx*padFactor*2;
    numelo= numel(fvx)/padFactor^2/padFactorTime;
    stept = round(padFactorTime*0);
    for ii = 1+stept:numel(nyu)-stept
        ind = find((abs(abs(kx+i*ky)-k(ii)))<dk);
%         indx = find((abs(abs(kx+i*ky-k(ii)))<dk)|(abs(abs(kx+i*ky+k(ii)))<dk));
%         indy = find((abs(abs(kx+i*ky-i*k(ii)))<dk)|(abs(abs(kx+i*ky+i*k(ii)))<dk));
%         tmp(ii) = numel(indx);
%         if(nyu(ii)==3)
% %             fu = fv;
%             fu1 = fu;
%             fu1(ii, indx) = fu1(ii,indx)*0;
%             surf(kx, ky,squeeze(fu(ii,:,:)+fu1(ii,:,:)), 'linestyle', 'none');
%             view(2)
%             caxis([-max(max(squeeze(fu(ii,:,:))))*0 max(max(squeeze(fu(ii,:,:))))])
%             [fup, fukx, fuky] = findPeak(squeeze(fu(ii, :, :)), kx, ky, 0, 0, 5,5);
%             disp([fukx, fuky])
%             [fvp, fvkx, fvky] = findPeak(squeeze(fv(ii, :, :)), kx, ky, 0, 0, 5,5);            
%             disp([fvkx, fvky])
%             clear all
%             error('1')
%         end
%         Au(ii) = sqrt(sum(sum(fvx(ii-stept:ii+stept,indx).^2,2)/numel(fvx)/numelo));
        Av(ii) = sqrt(sum(sum(fv(ii-stept:ii+stept,ind).^2,2)/numel(fv)/numelo));        
        
    end
  
%     plot(tmp) 
%     return
%     Au(1:stept) = 0;
%     Au(numel(nyu)-stept+1:numel(nyu)) = 0;
    Av(1:stept) = 0;
    Av(numel(nyu)-stept+1:numel(nyu)) = 0;

%     spx = Au';
%     spy = Av';
    Av = Av*sqrt(2); %becouse of normal fft of time
    frameWavesSpectra.vsp = Av';%(spx.^2+spy.^2).^.5;    
    
    
    frameWavesSpectra.nyu =nyu;
    
end

function [kx ky dkx dky] = genkxky(px, py, padFactor)

    nky = size(px,1)*padFactor;
    nkx = size(px,2)*padFactor;

    dkx = 2*pi/(px(1,2)-px(1,1))/nkx;
    dky = 2*pi/(py(2,1)-py(1,1))/nky;

    nyb = floor(nky/2);
    nxb = floor(nkx/2);
    nye = ceil(nky/2)-1;
    nxe = ceil(nkx/2)-1;
    [kx, ky] = meshgrid((-nxb:nxe)*dkx, (-nyb:nye)*dky);

end



