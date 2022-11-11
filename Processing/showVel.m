function showVel(frame, num)
    cutoff_percent = 0.01;

    if ~exist('frame', 'var')
        disp('function showVel(frame, num)');
        return
    end
    
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end    
      
    oVx = nanmean(cat(3, frame.vx{num}),3);
    oVy = nanmean(cat(3, frame.vy{num}),3);
    px = frame.px{num(1)};
    py = frame.py{num(1)};
    px = px - mean(px(:));
    py = py - mean(py(:));

    oVx = mean(oVx, 3);
    oVy = mean(oVy, 3);
    V = abs(oVx+i*oVy);
%     V = oVy;
    pr_right = 14;
    pr_top = 12;  
    surf(px, py , V, 'LineStyle', 'None');
    view(2);
    xlim([px(1) px(end)]);
    ylim([py(1) py(end)]);
    set(gca, 'fontsize', 15);
    xlabel('x, cm', 'fontsize', 15);
    ylabel('y, cm', 'fontsize', 15);  
    colormap(parula)
    hstep = (max(max(V) - min(min(V))))/1e4;
    h = hist(V(1:end), min(min(V)):hstep:max(max(V)));
    c = cumsum(h);
    xs = min(min(V)):hstep:max(max(V));
    xmin = xs(max(find(c < cutoff_percent*c(end))));
    xmax = xs(min(find(c > (1-cutoff_percent)*c(end))));
    xboundary = max(abs([xmin xmax]));
%     caxis([-xboundary(1) xboundary(1)] );
    caxis([0 xboundary(1)] );
    colorbar( 'fontsize', 15);
    set(gcf, 'PaperPosition', [0, 0, pr_right, pr_top]);