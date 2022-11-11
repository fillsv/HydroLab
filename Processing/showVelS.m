function showVelS(frame, num, zoom, mult)

    if ~exist('frame', 'var')
        disp('function showVelL(frame, num, zoom, mult)');
        return
    end
    
    if ~exist('zoom', 'var')
        zoom = 3;
        
    end
    
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end    
    
    if exist('mult', 'var')==0 
        mult = 1;
    end        
    oVx = nanmean(cat(3, frame.vx{num}),3);
    oVy = nanmean(cat(3, frame.vy{num}),3);
    px = frame.px{num(1)};
    py = frame.py{num(1)};
    px = px - mean(px(:));
    py = py - mean(py(:));

    sizex = floor(size(oVx,2)/mult)*mult;
    sizey = floor(size(oVx,1)/mult)*mult;
    oVx = avermult(oVx(1:sizey,1:sizex), mult);
    oVy = avermult(oVy(1:sizey,1:sizex), mult);
    px = avermult(px(1:sizey,1:sizex), mult);
    py = avermult(py(1:sizey,1:sizex), mult);
    
    scale = 1;
    %q = streamslice(px, py, 100*ones(size(px)), oVx*scale, oVy*scale,zeros(size(px)), 'r'); view(2)
    xmin = min(min(px));%-51; %
    xmax = max(max(px));%51; %
    deltax= (xmax-xmin)/(length(px)-1);
    x=xmin:deltax:xmax;
    
    ymin = min(min(py));%-51; %
    ymax = max(max(py));%51; %
    deltay= (ymax-ymin)/(length(py)-1);
    y=ymin:deltay:ymax;
    
    [xx, yy] = meshgrid(x,y);
    xx=double(xx);
    yy=double(yy);
%     xx=round(xx,4);
%     yy=round(yy,4);
    
    q = streamslice(xx, yy, oVx*scale, oVy*scale ,zoom); view(2)
    
%     frame = calcVort(frame, num);           
%     omega = cat(3, frame.omega{num});
    
    for ii=1:length(q)
    %zi = interp2(omega,q(ii).XData, q(ii).YData);
%     q(ii).ZData = zi;

    q(ii).ZData = 100*ones(size(q(ii).XData));
    end
    set(q,'LineWidth',1)
    set(q,'Color','r');
    %%q = quiver3(px, py, 100*ones(size(px)), oVx*scale, oVy*scale,zeros(size(px))); view(2)

    title(num2str(max(max(abs(oVx+i*oVy)))))
%     set(q, 'AutoScale', 'on')
%     set(q, 'AutoScaleFactor', .8)        
    xlim([-1 1]*frame.Lx/2);
    ylim([-1 1]*frame.Ly/2);
    set(gca, 'fontsize', 15);
    pr_right = 14;
    pr_top = 12;  
    set(gcf, 'PaperPosition', [0, 0, pr_right, pr_top]);

    xlabel('x, cm', 'fontsize', 15);
    ylabel('y, cm', 'fontsize', 15);  
%     drawnow;
%@D
