function showWaves(frame, nyuCrop, num)
    if ~exist('frame', 'var')
        disp('function showWaves(frame, nyuCrop, num)');
        return
    end
    cutoff_percent = 0.002;
    if ~isfield(frame, 'vox')
        if isfield(frame, 'vx')
            if exist('num', 'var')==0 
                num = 1:numel(frame.vx);
            end 
            frame = calcWaves(frame, 1, nyuCrop, num);
        else
            warning('frame is incorrect!')
            return
        end
    end

    px = frame.px;
    py = frame.py;
    Lx = frame.Lx;
    Ly = frame.Ly;
    
    rez = abs(abs(frame.vox)+i*abs(frame.voy));
    pr_right = 14;
    pr_top = 12;  
    surf(px-Lx/2, py-Ly/2, rez, 'LineStyle', 'None');
    view(2);
    xlim([px(1) px(end)]-Lx/2);
    ylim([py(1) py(end)]-Ly/2);
    set(gca, 'fontsize', 15);
    xlabel('x, cm', 'fontsize', 15);
    ylabel('y, cm', 'fontsize', 15);  
    colormap(parula)
    hstep = (max(max(rez) - min(min(rez))))/1e4;
    h = hist(rez(1:end), min(min(rez)):hstep:max(max(rez)));
    c = cumsum(h);
    xs = min(min(rez)):hstep:max(max(rez));
    xmin = xs(max(find(c < cutoff_percent*c(end))));
    xmax = xs(min(find(c > (1-cutoff_percent)*c(end))));
    xboundary = max(abs([xmin xmax]));
    caxis([0 xboundary(1)] );
    colorbar( 'fontsize', 15);
    title(nyuCrop)
    set(gcf, 'PaperPosition', [0, 0, pr_right, pr_top]);
