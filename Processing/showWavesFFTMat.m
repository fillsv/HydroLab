function showWavesFFTMat(frame, maxk, padFactor, nyuCrop, num)
    if ~exist('frame', 'var')
        disp('function showWavesFFT(frame, maxk, padFactor, nyuCrop, num)');
        return
    end
    cutoff_percent = 0.002;
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
%     cmax = 0;
    for ii = 1:nx;
        fprintf('%02d', ii)
        for jj = 1:ny;
            hs(ii,jj) = subplot(nx, ny, nx+1-ii+(ny-jj)*ny);   
            showWavesFFT(frame(jj, ii), maxk, padFactor, nyuCrop, num);
            xlabel('')
            ylabel('')
            set(gca, 'fontSize', 1)
            h = colorbar;
            cmax(ii, jj) = h.Limits(2);
            colorbar('off')
            fprintf('.')
%             if ny+1-ii+(nx-jj)*ny==1 
%                 disp([ii, jj]); 
%             end
        end
        fprintf('\n')
    end
    for ii = 1:nx*ny;
        caxis(hs(ii),[0 max(cmax(:))]);
    end
    drawnow;
    
    
