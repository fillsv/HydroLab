function showLyap8(frame, num, direction)
    if exist('frame', 'var')==0 
        disp('function frame = showLyap8(frame, num, direction)');
        return
    end
    if ~isfield(frame, 'lyap_b')
        warning('frame is incorrect!')
        return
    end    
    if ~exist('direction', 'var')
        direction = 'b';
    end

    px = frame.pxL{1};
    py = frame.pyL{1};
    px = px - mean(px(:));
    py = py - mean(py(:));
    lb = double(frame.lyap_b{num});
    lf = double(frame.lyap_f{num});
    lb(find(lb == 255)) = NaN;
    lf(find(lf == 255)) = NaN;
    lb = lb/254;
    lf = lf/254;
    lb = lb*(frame.lyap_lb_max(num)-frame.lyap_lb_min(num));
    lb = lb +frame.lyap_lb_min(num);
    lf = lf*(frame.lyap_lf_max(num)-frame.lyap_lf_min(num));
    lf = lf +frame.lyap_lf_min(num);
    l = zeros(size(lf));
    if ~isempty(strfind(direction, 'f'))
        l = l + lf;
    end
    if ~isempty(strfind(direction, 'b'))
        l = l + lb;
    end    
    surf(px, py, l, 'linestyle', 'none');view(2)    
    colormap('gray')
    set(gca, 'fontsize', 15);
    xlabel('x, cm', 'fontsize', 15);
    ylabel('y, cm', 'fontsize', 15);  
 
    xlim([-1 1]*frame.Lx/2)
    ylim([-1 1]*frame.Ly/2)
    colorbar
    