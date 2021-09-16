function cutframe = cutFrame(frame, xb, yb, num)
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
    px = frame.px{1};
    py = frame.py{1};
    
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
        cutframe.vx{kk} = reshape(frame.vx{ii}(ixy),[newysize, newxsize]);
        cutframe.vy{kk} = reshape(frame.vy{ii}(ixy),[newysize, newxsize]);
        cutframe.px{kk} = px;
        cutframe.py{kk} = py;
    end
%     disp(kk)
end