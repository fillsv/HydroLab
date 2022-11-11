function frame = showVortCoh(frame, num, keep, deep, wavelette)
    if exist('frame', 'var')==0 
        disp('function frame = showVortCoh(frame, num, keep, deep, wavelette)');
        return
    end
    if ~exist('keep', 'var')
        keep = 0;
        
    end
    if exist('num', 'var')==0 
        num = 1:numel(frame.vx);
    end 
    if exist('deep', 'var')==0 
        deep = 3;
    end
    if exist('wavelette', 'var')==0 
        wavelette = 'coif2';
    end
    
    clf;
    frameOmega = calcVort(frame,num);
    num = 1:numel(frameOmega.px);
    [frameCoh] = calcCoh_Incoh(frameOmega, num, keep, deep, wavelette );
    showVort(frameCoh, num);
%     hold on;
%     showVelS(frame, num, zoom, mult);

%@D