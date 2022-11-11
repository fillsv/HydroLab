function showWavesFFT(frame, maxk, padFactor, nyuCrop, num)
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
            frame = calcWaves(frame, padFactor, nyuCrop, num);
        else
            warning('frame is incorrect!')
            return
        end
    end


%     if exist('num', 'var')==0 
%         num = 1:numel(frame.fft2vox);
%     end 

    if exist('maxk', 'var')==0 
        maxk = 2;
    end 

    kx = frame.kx;
    ky = frame.ky;
    

    rez = abs(abs(frame.fft2vox)+i*abs(frame.fft2voy));

    ixy = find(abs(kx)<maxk & abs(ky)<maxk);
    pr_right = 14;
    pr_top = 12;  
    newxsize = length(find(abs(kx(1,:))<maxk));
    newysize = length(find(abs(ky(:,1))<maxk));
    newkx = reshape(kx(ixy),[newysize, newxsize]);
    newky = reshape(ky(ixy),[newysize, newxsize]);
    newrez = reshape(rez(ixy),[newysize, newxsize]);
    surf(newkx, newky, newrez, 'LineStyle', 'None'); view(2);
    set(gca, 'fontsize', 20);
    set(gcf, 'PaperPosition', [0, 0, pr_right, pr_top]);
    colorbar( 'fontsize', 20);
    xlim([-maxk maxk]);
    ylim([-maxk maxk]);   
    xlabel('k_x, cm^{-1}', 'fontsize', 20);
    ylabel('k_y, cm^{-1}', 'fontsize', 20);
        