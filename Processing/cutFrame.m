function cutframe = cutFrame(frame, xb, yb, num, center)
    if ~exist('frame', 'var')
        disp('function cutframe = cutFrame(frame, xb, yb, num)');
        return
    end
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end
    if exist('yb', 'var') == 0
        yb = xb;
    end
    if isempty(yb)
        yb = xb;
    end
    if isempty(num)
        num = 1:numel(frame.vx);
    end
        
    px = frame.px{1};
    py = frame.py{1};

    if exist('center', 'var')
        if center == 1
            px = px-mean(px, 'all');
            py = py-mean(py, 'all');
        end
    end

    
    cutframe = frame;
    cutframe.vx = [];
    cutframe.vy = [];
    cutframe.px = [];
    cutframe.py = [];
    cutframe.Lx = xb(2)-xb(1);
    cutframe.Ly = yb(2)-yb(1);

    ixy = find(px>=xb(1)&px<=xb(2) & py>=yb(1)&py<=yb(2));
    newxsize = length(find(px(1,:)>=xb(1)&px(1,:)<=xb(2)));
    newysize = length(find(py(:,1)>=yb(1)&py(:,1)<=yb(2)));
    px = reshape(px(ixy),[newysize, newxsize]);
    py = reshape(py(ixy),[newysize, newxsize]);
    kk =0;
    for ii = num
        kk = kk + 1;
        cutframe.vx{kk,1} = reshape(frame.vx{ii}(ixy),[newysize, newxsize]);
        cutframe.vy{kk,1} = reshape(frame.vy{ii}(ixy),[newysize, newxsize]);
        cutframe.px{kk,1} = px;
        cutframe.py{kk,1} = py;
    end
%     disp(kk)
end