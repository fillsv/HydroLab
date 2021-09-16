function h = showStream(frame, num, step, distance)

    if ~exist('frame', 'var')
        disp('function showStream(frame, num, step, distance)');
        return
    end
    
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end    
    
    if exist('step', 'var')==0 
        step = 6;
    end   
    
    if exist('distance', 'var')==0 
        distance = [.3 40];
    end   
    
    vx = double(nanmean(cat(3, frame.vx{num}),3));
    vy = double(nanmean(cat(3, frame.vy{num}),3));
    vz = zeros(size(vx));
    
    px = frame.px{num(1)};
    py = frame.py{num(1)};
    px = double(px - mean(px(:)));
    py = double(py - mean(py(:)));
    pz = 100*ones(size(px));
    
    px(:,:,2) = px;
    py(:,:,2) = py;
    pz(:,:,2) = pz+1;    
    vx(:,:,2) = vx;
    vy(:,:,2) = vy;
    vz(:,:,2) = vz;
    
    pxs = px(1:step:end,1:step:end);
    pys = py(1:step:end,1:step:end);
    pzs = 100*ones(size(pxs));
    
%     clf
    h = streamline(px, py, pz, vx, vy, vz, pxs, pys, pzs, distance);
    xlim([-1 1]*frame.Lx/2);
    ylim([-1 1]*frame.Ly/2);
    set(gca, 'fontsize', 15);
    pr_right = 14;
    pr_top = 12;  
    set(gcf, 'PaperPosition', [0, 0, pr_right, pr_top]);

    xlabel('x, cm', 'fontsize', 15);
    ylabel('y, cm', 'fontsize', 15);  

