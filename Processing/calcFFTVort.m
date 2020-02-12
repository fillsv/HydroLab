function frameFFTOmega = calcFFTVort(frame, padFactor, num)
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end    

    frameOmega = calcVort(frame, num)

    if(exist('padFactor', 'var')==0) 
        padFactor=4; 
    end
%     Lx = frame.Lx; %cm.
%     Ly = frame.Ly; %cm.
    dt = 1/frame.freq;

    px = frameOmega.px{1};
    py = frameOmega.py{1};    
    [npy, npx] = size(px);

    o = zeros(npy*padFactor, npx*padFactor);

    dkx = 2*pi/(px(1,2)-px(1,1))/size(o, 2);
    dky = 2*pi/(py(2,1)-py(1,1))/size(o, 1);

    nyb = floor(size(o, 1)/2);
    nxb = floor(size(o, 2)/2);
    nye = ceil(size(o, 1)/2)-1;
    nxe = ceil(size(o, 2)/2)-1;
    [kx, ky] = meshgrid((-nxb:nxe)*dkx, (-nyb:nye)*dky);
    
    frameFFTOmega.kx = kx;
    frameFFTOmega.ky = ky;

    w = window(@hann, npy)*window(@hann, npx)';
    jj = 0;
%     size(w)
    omega = cat(3, frameOmega.omega{:});
    omega(find(isnan(omega))) = 0;

    omega = omega(2:end, 2:end,:);
    [ox oy ot] = size(omega);
    fomega = (fftshift(fft2(omega.*w, padFactor*ox, padFactor*oy)/(ox*oy)/mean(w(:))));
    
    for ii = 1:numel(num)
        frameFFTOmega.fomega{ii} = fomega(:,:,ii);
    end