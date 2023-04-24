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
    dx = px(1,2)-px(1,1);
    dy = py(2,1)-py(1,1);
    Lx = xb(2)-xb(1);
    Lx = ceil(Lx/dx)*dx; %Округяем до челого числа шагов
    cenX = ((xb(2)-xb(1))/2) + xb(1);
    xb(2) = cenX + Lx/2;
    xb(1) = cenX - Lx/2;
    
    Ly = yb(2)-yb(1);
    Ly = ceil(Ly/dy)*dy; %Округяем до челого числа шагов
    cenY = ((yb(2)-yb(1))/2) + yb(1);
    yb(2) = cenY + Ly/2;
    yb(1) = cenY - Ly/2;
    
    cutframe.Lx = xb(2)-xb(1);
    cutframe.Ly = yb(2)-yb(1);

%     ixy = find(px>=xb(1)&px<=xb(2) & py>=yb(1)&py<=yb(2));
    ixy = find(px>xb(1)&px<xb(2) & py>yb(1)&py<yb(2));
%     newxsize = length(find(px(1,:)>=xb(1)&px(1,:)<=xb(2)));
%     newysize = length(find(py(:,1)>=yb(1)&py(:,1)<=yb(2)));
    newxsize = length(find(px(1,:)>xb(1)&px(1,:)<xb(2)));
    newysize = length(find(py(:,1)>yb(1)&py(:,1)<yb(2)));
    px = reshape(px(ixy),[newysize, newxsize]);
    py = reshape(py(ixy),[newysize, newxsize]);
    kk =0;
    
    dx = px(1,2)-px(1,1);
    dy = py(2,1)-py(1,1);
    lineRX = [-300:300]*dx;
    lineRY = [-300:300]*dy;
%     lineRX = fix(lineRX*100);
%     lineRY = fix(lineRY*100);
%     lineRX = round(lineRX,2);
%     lineRY = round(lineRY,2);
%     ny = find(lineRX>round(xb(1),2) & lineRX<=round(round(xb(2),3),2));
%     nx = find(lineRY>round(yb(1),2) & lineRY<=round(round(yb(2),3),2));
%     ny = find(lineRX>fix(xb(1)*100) & lineRX<=fix(xb(2)*100));
%     nx = find(lineRY>fix(yb(1)*100) & lineRY<=fix(yb(2)*100));
%     ny = find(lineRX>xb(1) & lineRX<=xb(2));
%     nx = find(lineRY>yb(1) & lineRY<=yb(2));
%     nx = length(nx);
%     ny = length(ny);
    nx = round(Ly/dy);
    ny = round(Lx/dx);
%     disp(nx +":"+ ny);
    
%     nx = round(cutframe.Lx/dx);
%     ny = round(cutframe.Ly/dy);
    px_tmp = zeros(nx,ny)/0;
    py_tmp = zeros(nx,ny)/0;
    
    if(xb(2)-max(px(1,:)) < dx && yb(2)-max(py(:,1)) < dy)
        %frame in up right
        px_tmp(end-size(px,1)+1:end, end-size(px,2)+1:end) = px;
        ax = [size(px_tmp,1)-size(px,1)+1:size(px_tmp,1)];
        bx = [size(px_tmp,2)-size(px,2)+1:size(px_tmp,2)];
        py_tmp(end-size(py,1)+1:end, end-size(py,2)+1:end) = py;
        ay = [size(py_tmp,1)-size(py,1)+1:size(py_tmp,1)];
        by = [size(py_tmp,2)-size(py,2)+1:size(py_tmp,2)];
    elseif(xb(2)-max(px(1,:)) < dx && min(py(:,1))-yb(1) < dy)
        %frame in down right
        px_tmp(1:size(px,1), end-size(px,2)+1:end) = px;
        ax = [1:size(px,1)];
        bx = [size(px_tmp,2)-size(px,2)+1:size(px_tmp,2)];
        py_tmp(1:size(py,1), end-size(py,2)+1:end) = py;
        ay = [1:size(py,1)];
        by = [size(py_tmp,2)-size(py,2)+1:size(py_tmp,2)];
    elseif(min(px(1,:))-xb(1) < dx && min(py(:,1))-yb(1) < dy)
        %frame in down left        
        px_tmp(1:size(px,1), 1:size(px,2)) = px;
        ax = [1:size(px,1)];
        bx = [1:size(px,2)];
        py_tmp(1:size(py,1), 1:size(py,2)) = py;
        ay = [1:size(py,1)];
        by = [1:size(py,2)];
    elseif(min(px(1,:))-xb(1) < dx && yb(2)-max(py(:,1)) < dy)
        %frame in up left
        px_tmp(end-size(px,1)+1:end, 1:size(px,2)) = px;
        ax = [size(px_tmp,1)-size(px,1)+1:size(px_tmp,1)];
        bx = [1:size(px,2)];
        py_tmp(end-size(py,1)+1:end, 1:size(py,2)) = py;
        ay = [size(py_tmp,1)-size(py,1)+1:size(py_tmp,1)];
        by = [1:size(py,2)];
    end    
        b_px = nanmean(px_tmp,1);
        if(isnan(b_px(1)))
            for(ttt = length(b_px):-1:1)
                if(isnan(b_px(ttt)))
                   px_tmp(:,ttt) = b_px(ttt+1)-dx;  
                   b_px(ttt) = b_px(ttt+1)-dx;
                else
                   px_tmp(:,ttt) = b_px(ttt);
                end
            end
        else
            for(ttt = 1:length(b_px))
                if(isnan(b_px(ttt)))
                   px_tmp(:,ttt) = b_px(ttt-1)+dx;
                   b_px(ttt) = b_px(ttt-1)+dx;
                else
                   px_tmp(:,ttt) = b_px(ttt);
                end
            end
        end
        
        b_py = nanmean(py_tmp,2);
        if(isnan(b_py(1)))
            for(ttt = length(b_py):-1:1)
                if(isnan(b_py(ttt)))
                   py_tmp(ttt,:) = b_py(ttt+1)-dx;  
                   b_py(ttt) = b_py(ttt+1)-dx;
                else
                   py_tmp(ttt,:) = b_py(ttt);
                end
            end
        else
            for(ttt = 1:length(b_py))
                if(isnan(b_py(ttt)))
                   py_tmp(ttt,:) = b_py(ttt-1)+dx;
                   b_py(ttt) = b_py(ttt-1)+dx;
                else
                   py_tmp(ttt,:) = b_py(ttt);
                end
            end
        end
        
        px = px_tmp;
        py = py_tmp;
    
    
    for ii = num
        kk = kk + 1;
        cutframe.vx{kk,1} = zeros(nx,ny)/0;
        cutframe.vx{kk,1}(ax,bx) = reshape(frame.vx{ii}(ixy),[newysize, newxsize]);
%         cutframe.vx{kk,1} = reshape(frame.vx{ii}(ixy),[newysize, newxsize]);
        cutframe.vy{kk,1} = zeros(nx,ny)/0;
        cutframe.vy{kk,1}(ay,by) = reshape(frame.vy{ii}(ixy),[newysize, newxsize]);
%         cutframe.vy{kk,1} = reshape(frame.vy{ii}(ixy),[newysize, newxsize]);
        cutframe.px{kk,1} = px;
        cutframe.py{kk,1} = py;
    end
%     disp(kk)
end