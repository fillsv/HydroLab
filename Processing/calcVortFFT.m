function frameFFTVort = calcVortFFT(frame, padFactor, num)

    if ~exist('frame', 'var')
        disp('function frameVortFFT = calcFFTVort(frame, padFactor, num)');
        return
    end
    
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end    
    if ~isfield(frame, 'omega')
        frameOmega = calcVort(frame, num);
    else
        frameOmega = frame;
    end

    if(exist('padFactor', 'var')==0) 
        padFactor=4; 
    end
    dt = 1/frame.freq;

    px = frameOmega.px{num(1)};
    py = frameOmega.py{num(1)};    
    [npy, npx] = size(px);


    [kx ky] = genkxky(px, py, padFactor);

    frameFFTVort.freq = frame.freq;    
    frameFFTVort.Lx = frame.Lx;    
    frameFFTVort.Ly = frame.Ly;    
    frameFFTVort.kx = kx;
    frameFFTVort.ky = ky;
    frameFFTVort.padFactor = padFactor;
    
    w = window(@hann, npy)*window(@hann, npx)';
    omega = cat(3, frameOmega.omega{num});
    omega(find(isnan(omega))) = 0;

    [ox oy ot] = size(omega);
    fomega = fftshift(fft2(omega.*w, padFactor*ox, padFactor*oy)/(ox*oy)/mean(w(:)));
    
    for ii = 1:numel(num)
        frameFFTVort.fomega{ii} = fomega(:,:,ii);
    end
    frameFFTVort.tt(1:numel(num)) = frame.tt(num);
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