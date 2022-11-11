function aframe = averFrame(frame, num, mult)

    if ~exist('frame', 'var')
        disp('function averFrame(frame, num, mult)');
        return
    end
    
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end
    if isempty(num)
        num = 1:numel(frame.px);
    end
        
    aframe = sumVel(frame,1, num);
    kk = 0;
    for ii = 1:numel(num)
        kk = kk + 1;
        oVx = frame.vx{ii};
        oVy = frame.vy{ii};
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
        aframe.px{kk} = px;
        aframe.py{kk} = py;
        aframe.vx{kk} = oVx;
        aframe.vy{kk} = oVy;
    end
