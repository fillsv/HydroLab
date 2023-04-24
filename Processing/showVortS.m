function frame = showVortS(frame, num, zoom, mult)
    if exist('frame', 'var')==0 
        disp('function frame = showVortQ(frame, num, zoom, mult)');
        return
    end
    if ~exist('zoom', 'var')
        zoom = 3;
        
    end
    if exist('num', 'var')==0 
        num = 1:numel(frame.vx);
    end 
    if exist('mult', 'var')==0 
        mult = 4;
    end
    clf;
    showVort(frame, num);
    hold on;
    showVelS(frame, num, zoom, mult);
    hold off;
    %@D