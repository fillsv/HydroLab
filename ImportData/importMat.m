function frame = importMat(x,y,u,v,matInfo)
    freq = matInfo.freq;
    Lx = matInfo.Lx;
    Ly = matInfo.Ly;

    dt = 1/freq;
    
    px = x{1};
    py = y{1};
    px = px - px(1, 1) + (px(1, 2) - px(1, 1));
    py = py - py(1, 1) + (py(2, 1) - py(1, 1));
    maxpx = max(px(:));
    maxpy = max(py(:));
    filterLimit = 3.5;
    px = px*Lx/maxpx;
    py = py*Ly/maxpy;
 
    for ii = 1:numel(x)
        vx = u{ii};
        vy = v{ii};
       
        vx(find(isinf(vx)))=NaN;
        vy(find(isinf(vy)))=NaN;
        if exist('filterLimit', 'var') ~= 0 
            if filterLimit > 0
            vabs = abs(vx+i*vy);
            f5 = find(vabs>nanmean(vabs(:))*filterLimit);
            vx(f5) = NaN;
            vy(f5) = NaN;
            end
        end

        vx = vx*Lx/maxpx/dt;
        vy = vy*Ly/maxpy/dt;
        frame.px{ii,1} = px;
        frame.py{ii,1} = py;
        frame.vx{ii,1} = vx;
        frame.vy{ii,1} = vy;
    end
    frame.tt = matInfo.tt;
    frame.freq = freq;
    
