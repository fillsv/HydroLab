function frame = correctFrame(frame)
    px = frame.px{1};
    py = frame.py{1};
%     py = py - mean(py(:,1));
%     px = px - mean(px(1,:));
    dx = px(1,2) - px(1,1);
    dy = py(2,1) - py(1,1);
    dt = 1/75;%frame.freq;
    num = 1:numel(frame.vx);
    for nn = num
        vx = frame.vx{nn};
        vy = frame.vy{nn};
        dxvx = (vx(2:end-1, 3:end) - vx(2:end-1, 1:end-2))/2/dx; 
        dxvy = (vy(2:end-1, 3:end) - vy(2:end-1, 1:end-2))/2/dx;
        dyvx = (vx(3:end, 2:end-1) - vx(1:end-2, 2:end-1))/2/dy;
        dyvy = (vy(3:end, 2:end-1) - vy(1:end-2, 2:end-1))/2/dy;
        vx(2:end-1, 2:end-1) = vx(2:end-1, 2:end-1) - dt/2*...
            (vx(2:end-1, 2:end-1).*dxvx+vy(2:end-1, 2:end-1).*dyvx);
        vy(2:end-1, 2:end-1) = vy(2:end-1, 2:end-1) - dt/2*...
            (vx(2:end-1, 2:end-1).*dxvy+vy(2:end-1, 2:end-1).*dyvy);
        frame.vx{nn} = vx;
        frame.vy{nn} = vy;        
    end