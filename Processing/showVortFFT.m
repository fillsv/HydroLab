function showFFTVort(frame, maxk, padFactor, num)
    if ~exist('frame', 'var')
        disp('function showFFTVort(frame, maxk, padFactor, num)');
        return
    end
    if exist('padFactor', 'var')==0 
        padFactor = 4;
    end 

    if ~isfield(frame, 'vox')
        if isfield(frame, 'vx')
            if exist('num', 'var')==0 
                num = 1:numel(frame.vx);
            end 
            frame = calcVortFFT(frame, padFactor, num);
        else
            warning('frame is incorrect!')
            return
        end
    end

    if exist('num', 'var')==0 
        num = 1:numel(frame.fomega);
    end 

    if exist('maxk', 'var')==0 
        maxk = 2;
    end 

    kx = frame.kx;
    ky = frame.ky;
    
    rez = cat(3, frame.fomega{:});
    rez = abs(mean(rez, 3));

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
        