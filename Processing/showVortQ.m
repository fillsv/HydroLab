function frame = showVortQ(frame, num, mult, arrowParam)
    if exist('frame', 'var')==0 
        disp('function frame = showVortQ(frame, num, mult)');
        return
    end
    if exist('num', 'var')==0 
        num = 1:numel(frame.vx);
    end 
    if isempty(num)
        num = 1:numel(frame.vx);
    end
    if exist('mult', 'var')==0 
        mult = 4;
    end
    if exist('arrowParam', 'var')==0 
        arrowParam = [];
    end
    showVort(frame, num);
    hold on;
    showVelQ(frame, num, mult, 1, arrowParam);
    hold off;