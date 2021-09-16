function divideFrame = divideFrame(frame, xynum, mult,  num)
    if ~exist('frame', 'var')
        disp('function divideFrame = divideFrame(frame, xnum, ynum, num)');
        return
    end
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end
%     if exist('yb', 'var') == 0
%         yb = xb;
%     end
    px = frame.px{num(1)};
    py = frame.py{num(1)};
    
    [ny nx] = size(px);
    xnum = xynum;
    ynum = xynum;
    multx = mult;
    multy = mult;
    [xb step_nx] = divxy(nx, xnum, multx);
    [yb step_ny] = divxy(ny, ynum, multy);
%     step_ny = floor(ny/ynum);
    
    for jj = 1:numel(yb)
        for ii = 1:numel(xb)
            xbound = [px(1,1+xb(ii)) px(1,xb(ii)+step_nx-1)];
            ybound = [py(1+yb(jj), 1) py(step_ny+yb(jj)-1, 1)];
            divideFrame(jj, ii) = cutFrame(frame,xbound, ybound, num);
        end
    end
end


function [xb xstep] = divxy(nx, xnum, mult)
    xstep = round(nx/xnum);
    xnum_new = ((xnum-1)*mult+1);
    nx_new = (nx+1-xstep);
    xb = round((0:xnum_new-1)/xnum_new*nx_new)+1;
end
