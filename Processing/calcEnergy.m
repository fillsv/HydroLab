function frameEnergy = calcEnergy(frame, padFactor, num)
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end    

    if(exist('padFactor', 'var')==0) 
        padFactor=4; 
    end

    px = frame.px{1};
    py = frame.py{1};    
    [npy, npx] = size(px);

    dx = px(1,2)-px(1,1);
    dy = py(2,1)-py(1,1);
    
    nky = size(px,1)*padFactor;
    nkx = size(px,2)*padFactor;
    
    [kx ky] = genkxky(px, py, padFactor);

    frameEnergy.kx = kx;
    frameEnergy.ky = ky;
    frameEnergy.Lx = frame.Lx;    
    frameEnergy.Ly = frame.Ly;    
    frameEnergy.freq = frame.freq;    
    frameEnergy.padFactor = padFactor;
    
    w = window(@hann, npy)*window(@hann, npx)';

    vx = cat(3, frame.vx{num});
    vy = cat(3, frame.vy{num});

    vx(find(isnan(vx))) = 0;
    vy(find(isnan(vy))) = 0;
%     E = zeros([nky, nkx, numel(num)]);
    for ii = 1:numel(num)
        
        fft2vx = fftshift(abs(fft2(vx(:,:,ii).*w, nky, nkx)))*dx*dy/2/pi/sqrt(mean(w(:).^2));
        fft2vy = fftshift(abs(fft2(vy(:,:,ii).*w, nky, nkx)))*dx*dy/2/pi/sqrt(mean(w(:).^2));

        frameEnergy.E{ii} = .5*(fft2vx.^2 + fft2vy.^2);

    end
    frameEnergy.tt(1:numel(num)) = frame.tt(num);
    clear frame fft2vx fft2vy vx vy w kx ky px py
    return
%     frameEnergy.E{1} = E;
    
%     fft2vx = fftshift(abs(fft2(vx.*w, nky, nkx)))*sqrt(dpx*dpy)/sqrt(mean(w(:).^2));
%     fft2vy = fftshift(abs(fft2(vy.*w, nky, nkx)))*sqrt(dpx*dpy)/sqrt(mean(w(:).^2));
%     E = .5*(fft2vx.^2 + fft2vy.^2);
%     for ii = 1:numel(num)
%         frameEnergy.E{ii} = E(:,:,ii);
%     end
    
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