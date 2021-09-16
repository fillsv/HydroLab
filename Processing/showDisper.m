function showDisper(frame, maxk, padFactor, nyuCrop, kstep, num)

    cutoff_percent = 0.002;
    if ~isfield(frame, 'vok')
        if isfield(frame, 'vx')
            if exist('num', 'var')==0 
                num = 1:numel(frame.vx);
            end 
            frame = calcDisper(frame, padFactor, nyuCrop, kstep, num);
        else
            warning('frame is incorrect!')
            return
        end
    end

    if exist('maxk', 'var')==0 
        maxk = 20;
    end 

    k = frame.k;
    nyu = frame.nyu;
    

    rez = frame.vok;

    ixy = find(abs(k)<maxk);
    pr_right = 14;
    pr_top = 12;  
    newnyusize = size(nyu, 1);
    newksize = length(find(abs(k(1,:))<maxk));
    newnyu = reshape(nyu(ixy),[newnyusize, newksize]);
    newk = reshape(k(ixy),[newnyusize, newksize]);
    newrez = reshape(rez(ixy),[newnyusize, newksize]);
    surf(newnyu, newk, log10(newrez), 'LineStyle', 'None'); view(2);
    set(gca, 'fontsize', 20);
    set(gcf, 'PaperPosition', [0, 0, pr_right, pr_top]);
    colorbar( 'fontsize', 20);
    xlim([min(newnyu(:)) max(newnyu(:))]);
    ylim([min(newk(:)) max(newk(:))]);   
    xlabel('Freq, Hz', 'fontsize', 20);
    ylabel('k, cm^{-1}', 'fontsize', 20);
        