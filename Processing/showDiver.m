function frame = showDiver(frame, num)
    cutoff_percent = 0.01;
    if exist('frame', 'var')==0 
        disp('function frame = showVort(frame, num)');
        return
    end
    if ~isfield(frame, 'diver')
        if isfield(frame, 'vx')
            if exist('num', 'var')==0 
                num = 1:numel(frame.vx);
            end 
            frame = calcDiver(frame, num);
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
    diver = cat(3, frame.diver{num});
    diver = nanmean(diver, 3);    
%     omega = abs(omega);
    pr_right = 14;
    pr_top = 12;  
    surf(px, py , diver, 'LineStyle', 'None');
    view(2);
    xlim([px(1) px(end)]);
    ylim([py(1) py(end)]);
    set(gca, 'fontsize', 15);
    xlabel('x, cm', 'fontsize', 15);
    ylabel('y, cm', 'fontsize', 15);  
    colormap(parula)
    hstep = (max(max(diver) - min(min(diver))))/1e4;
    h = hist(diver(1:end), min(min(diver)):hstep:max(max(diver)));
    c = cumsum(h);
    xs = min(min(diver)):hstep:max(max(diver));
    xmin = xs(max(find(c < cutoff_percent*c(end))));
    xmax = xs(min(find(c > (1-cutoff_percent)*c(end))));
    xboundary = max(abs([xmin xmax]));
    caxis([-xboundary(1) xboundary(1)] );
%     caxis([0 xboundary(1)] );
    colorbar( 'fontsize', 15);
    set(gcf, 'PaperPosition', [0, 0, pr_right, pr_top]);

