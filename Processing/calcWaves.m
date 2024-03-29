function frameWaves = calcWaves(frame, padFactor, nyuCrop, num)
    if ~exist('frame', 'var')
        disp('frameWaves = calcWaves(frame, padFactor, nyuCrop, num)');
        return
    end
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

    frameWaves.nyuCrop = nyuCrop;
    frameWaves.px = px;
    frameWaves.py = py;    
    frameWaves.kx = kx;
    frameWaves.ky = ky;
    frameWaves.Lx = frame.Lx;    
    frameWaves.Ly = frame.Ly;    
    frameWaves.freq = frame.freq;    
    frameWaves.padFactor = padFactor;
    frameWaves.tt(1:numel(num)) = frame.tt(num);
    
    vx = cat(3, frame.vx{num});
    vy = cat(3, frame.vy{num});

    vx = permute(vx, [3 1 2]);
    vy = permute(vy, [3 1 2]);
    
    w = repmat(hann(size(vx,1)), [1 size(vx,2) size(vx,3)]);
    
    vx(find(isnan(vx))) = 0;
    vy(find(isnan(vy))) = 0;
    
    pp = 1 + nyuCrop./frame.freq*numel(num);

%     k = disper(nyuCrop, sigma, rho, g, h);
    w23 = hann(size(vx, 2))*hann(size(vx,3))'; 
    w1 = hann(size(vx,1));
    w = repmat(w1, [1, size(vx,2), size(vx,3)]);
    vox = goertzel1(vx.*w, pp)/mean(w1)/size(vx,1);
    voy = goertzel1(vy.*w, pp)/mean(w1)/size(vx,1);

    vox = 2*squeeze(vox);
    voy = 2*squeeze(voy);    

    voxy = abs((abs(vox)+i*abs(voy)).^2);
    fft2voxy = fft2(voxy.*w23, padFactor*size(vx,2), padFactor*size(vx,3));
    fft2voxy = fft2voxy/mean(w23(:));
    fft2voxy = fftshift(fft2voxy);
    fft2voxy = fft2voxy/size(vx,2)/size(vx,3);
    
    fft2vox = fft2(vox.*w23, padFactor*size(vx,2), padFactor*size(vx,3));
    fft2vox = fft2vox/mean(w23(:));
    fft2vox = fftshift(fft2vox);
    fft2vox = fft2vox/numel(vx);
    fft2voy = fft2(voy.*w23, padFactor*size(vx,2), padFactor*size(vx,3));
    fft2voy = fft2voy/mean(w23(:));
    fft2voy = fftshift(fft2voy);
    fft2voy = fft2voy/numel(vy);    

    frameWaves.fft2voxy = fft2voxy;
    frameWaves.fft2vox = fft2vox;
    frameWaves.fft2voy = fft2voy;
    frameWaves.vox = vox;
    frameWaves.voy = voy;
    frameWaves.voxy = voxy;
    
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



