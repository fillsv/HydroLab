function showamplFFT(amplFFT)
    if ~exist('amplFFT', 'var')
        disp('function showamplFFT(amplFFT)');
        return
    end
    amplFFT(2:end+1, 2:end+1,:,:) = amplFFT;
%     amplFFT(1, end+1,:,:) = 0;
    cmax = 0;
    for jj = 1:2
        for ii = 1:2
            hs(ii,jj) = subplot(2,2,ii+(2-jj)*2);
            surf(amplFFT(end:-1:1,end:-1:1,jj,ii), 'linestyle', 'none'); view(2);
            h = colorbar;
            cmax = max([h.Limits(2) cmax]);
            colorbar('off')
            xlim([1 size(amplFFT,2)])
            ylim([1 size(amplFFT,2)])
        end
    end
    
    for ii = 1:4;
        caxis(hs(ii),[0 cmax]);
    end

    drawnow;