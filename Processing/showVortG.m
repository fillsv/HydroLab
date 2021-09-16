function frame = showVortG(frame, num)
    cutoff_percent = 0.01;
    if exist('frame', 'var')==0 
        disp('function frame = showVortG(frame, num)');
        return
    end
    if ~isfield(frame, 'omega')
        if isfield(frame, 'vx')
            if exist('num', 'var')==0 
                num = 1:numel(frame.vx);
            end 
            frame = calcVort(frame, num);
            num = 1:numel(frame.px);
        else
            warning('frame is incorrect!')
            return
        end
    end    
    
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end
    px = frame.px{1};%(1:end-1, 1:end-1);
    py = frame.py{1};%(1:end-1, 1:end-1);
    px = px - mean(px(:));
    py = py - mean(py(:));    
    omega = cat(3, frame.omega{num});
    omega = nanmean(omega, 3);
    sigma = 2;
    omega = imgaussfilt(omega, sigma, 'padding', 'symmetric');    
%     omega = abs(omega);
    pr_right = 14;
    pr_top = 12;  
    surf(px, py , omega, 'LineStyle', 'None');
    view(2);
    xlim([px(1) px(end)]);
    ylim([py(1) py(end)]);
    set(gca, 'fontsize', 15);
    xlabel('x, cm', 'fontsize', 15);
    ylabel('y, cm', 'fontsize', 15);  
    colormap(parula)
    hstep = (max(max(omega) - min(min(omega))))/1e4;
    h = hist(omega(1:end), min(min(omega)):hstep:max(max(omega)));
    c = cumsum(h);
    xs = min(min(omega)):hstep:max(max(omega));
    xmin = xs(max(find(c < cutoff_percent*c(end))));
    xmax = xs(min(find(c > (1-cutoff_percent)*c(end))));
    xboundary = max(abs([xmin xmax]));
    caxis([-xboundary(1) xboundary(1)] );
%     caxis([0 xboundary(1)] );
    colorbar( 'fontsize', 15);
    set(gcf, 'PaperPosition', [0, 0, pr_right, pr_top]);


