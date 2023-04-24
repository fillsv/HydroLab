function frame = showLyap(frame, num, direction)
    if exist('frame', 'var')==0 
        disp('function frame = showVort(frame, num, direction))');
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
%     py = -py;
    py = py - mean(py(:));
    lb = squeeze(frame.lyap_b{num});
    lf = squeeze(frame.lyap_f{num});
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
   
%     hold on
%     showVelQ(frame, 250, 5);
%     hold off

    % colormap('palura')
 
 