function amplFFT = calcWavesFFTMat(frame, padFactor, nyuCrop, abskx, absky, num)
    if ~exist('frame', 'var')
        disp('function calcWavesFFT(frame, padFactor, nyuCrop, abskx, absky, num)');
        return
    end
    if ~isfield(frame, 'fft2vox')
        if isfield(frame, 'vx')
            if ~exist('num', 'var')
                num = 1:numel(frame.vx);
            end 
        else
            warning('frame is incorrect!')
            return
        end
    end
    
%     frame = frame';
    [ny nx] = size(frame);
    cmax = 0;
    for ikx = [1 2]
        for iky = [1 2]
            for ii = 1:nx;
                for jj = 1:ny;
                    fr = calcWaves(frame(jj, ii), padFactor, nyuCrop, num);
                    value_kx = abskx*(ikx-1.5)*2;
                    value_ky = absky*(iky-1.5)*2;                    
                    ind = findnear2(fr.kx, fr.ky, value_kx, value_ky);
                    f = abs(abs(fr.fft2vox(ind)) + i*abs(fr.fft2voy(ind)));
                    amplFFT(jj,ii,iky,ikx) = f;
%                     disp([ikx iky ind value_kx value_ky])
                end
            end
        end
    end
    
    
    
function ind = findnear2(kx, ky, value_kx, value_ky)
    ind = find(min(abs(kx(:)-value_kx + i*(ky(:)-value_ky)))==abs(kx(:)-value_kx + i*(ky(:)-value_ky)));
    ind = ind(1);
    