function frameDisper = calcDisper(frame, padFactor, nyuCrop, kstep, num)
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end    

    if(exist('padFactor', 'var')==0) 
        padFactor=4; 
    end

    px = frame.px{num(1)};
    py = frame.py{num(1)};    
    
    [npy, npx] = size(px);

    sigma = 73;
    rho = 1;
    g = 981.9;
    h = 10;
    
    [kx ky] = genkxky(px, py, padFactor);
%     frameDisper.px = px;
%     frameDisper.py = py;    
%     frameDisper.kx = kx;
%     frameDisper.ky = ky;
    frameDisper.Lx = frame.Lx;    
    frameDisper.Ly = frame.Ly;    
    frameDisper.freq = frame.freq;    
    frameDisper.padFactor = padFactor;
    frameDisper.tt(1:numel(num)) = frame.tt(num);
    vx = cat(3, frame.vx{num});
    vy = cat(3, frame.vy{num});
    ll = 0;
    vx = permute(vx, [3 1 2]);
    vy = permute(vy, [3 1 2]);
    
    w = repmat(hann(size(vx,1)), [1 size(vx,2) size(vx,3)]);
    w23 = hann(size(vx, 2))*hann(size(vx,3))'; 
    
    vx(find(isnan(vx))) = 0;
    vy(find(isnan(vy))) = 0;
    kmax = max(kx(:));
    numk = floor(kmax/kstep);
    frameDisper.vok = zeros(numel(nyuCrop), numk);
    frameDisper.nyuCrop = nyuCrop;
    frameDisper.nyu = repmat(nyuCrop', [1 numk]);

    for nyuCrop = nyuCrop
        ll = ll +1;
        disp(nyuCrop)
        pp = 1 + nyuCrop./frame.freq*numel(num);

        vox = goertzel1(vx.*w, pp)/mean(w(:));
        voy = goertzel1(vy.*w, pp)/mean(w(:));

        vox = 2*squeeze(vox);
        voy = 2*squeeze(voy);    

        fft2vox = fft2(vox.*w23, padFactor*size(vx,2), padFactor*size(vx,3));
        fft2vox = fft2vox/mean(w23(:));
        fft2vox = fftshift(fft2vox);
        fft2vox = fft2vox/numel(vx);
        fft2voy = fft2(voy.*w23, padFactor*size(vx,2), padFactor*size(vx,3));
        fft2voy = fft2voy/mean(w23(:));
        fft2voy = fftshift(fft2voy);
        fft2voy = fft2voy/numel(vy);    
        fft2vo = abs(abs(fft2vox)+i*abs(fft2voy));
        [mrfft2vo k] = sumRadius(fft2vo, kx, ky, kstep);
        frameDisper.vok(ll, :) = mrfft2vo;
%         frameDisper.k = k;
%     fft2ox = fft2ox;%/(k(pp)*d*(n-1)/n);%squeeze(fu(pp,:,:));
%     fft2oy = fft2oy;%/(k(pp)*d*(n-1)/n);%squeeze(fv(pp,:,:));    
    end
    frameDisper.k = repmat(k, [size(frameDisper.vok, 1) 1]);
    
%     frameWaves.fft2voxy = fft2voxy;
%     frameWaves.fft2vox = fft2vox;
%     frameWaves.fft2voy = fft2voy;
%     frameWaves.vox = vox;
%     frameWaves.voy = voy;
%     frameWaves.voxy = voxy;
%     
end

function [vok, k] = sumRadius(vo, kx, ky, kstep)
    %
    
    kmax = max(kx(:));
    kr = (kx.^2+ky.^2).^.5;
    [k_temp, k_ind] = sort(kr(:));
    numk = floor(kmax/kstep);
    k = ((1:numk)-.5)*kstep;
    vok = zeros(numk,1);
    num_step = 2*pi*numel(kr)/(numk^2+numk)/4;
    for ii = 1:numk
        in = k_ind(round(num_step*ii*(ii-1)/2+1):round(num_step*ii*(ii+1)/2));
        vok(ii) = max(vo(in));
    end

end

function [kx ky] = genkxky(px, py, padFactor)

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



