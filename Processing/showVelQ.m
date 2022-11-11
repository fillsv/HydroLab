function showVelQ(frame, num, mult, center)

    if ~exist('frame', 'var')
        disp('function showVelQ(frame, num, mult,center)');
        return
    end
    
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end    
    if exist('center', 'var')==0 
        center = 1;
    end    
    
    if exist('mult', 'var')==0 
        mult = 4;
    end        
    oVx = nanmean(cat(3, frame.vx{num}),3);
    oVy = nanmean(cat(3, frame.vy{num}),3);
    px = frame.px{num(1)};
    py = frame.py{num(1)};
    if center
        px = px - mean(px(:));
        py = py - mean(py(:));
    end
    sizex = floor(size(oVx,2)/mult)*mult;
    sizey = floor(size(oVx,1)/mult)*mult;
    oVx = avermult(oVx(1:sizey,1:sizex), mult);
    oVy = avermult(oVy(1:sizey,1:sizex), mult);
    px = avermult(px(1:sizey,1:sizex), mult);
    py = avermult(py(1:sizey,1:sizex), mult);
    
    scale = 1;
    q = quiver3(px, py, 100*ones(size(px)), oVx*scale, oVy*scale,zeros(size(px)), 'r'); view(2)
%     title(num2str(max(max(abs(oVx+i*oVy)))))
    set(q, 'AutoScale', 'on')
    set(q, 'AutoScaleFactor', .8)  
    if center
        xlim([-1 1]*frame.Lx/2);
        ylim([-1 1]*frame.Ly/2);
    else
        xlim([min(px(:)) max(px(:))])
        ylim([min(py(:)) max(py(:))])
    end
    set(gca, 'fontsize', 15);
    pr_right = 14;
    pr_top = 12;  
    set(gcf, 'PaperPosition', [0, 0, pr_right, pr_top]);

    xlabel('x, cm', 'fontsize', 15);
    ylabel('y, cm', 'fontsize', 15);  
%     drawnow;

