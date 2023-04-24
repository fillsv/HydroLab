function frameL = calcLyap8(frame, num, dn, padFactorSpace, padFactorTime)
    if ~exist('dn', 'var') dn = 25; end
    if ~exist('num', 'var') num = []; end
    if isempty(num) 
        num = dn:numel(frame.px)-dn+1
    end
    px = (frame.px{num(1)});
    py = (frame.py{num(1)});
    minPx = min(px(1,:));
    maxPx = max(px(1,:));
    minPy = min(py(:,1));
    maxPy = max(py(:,1));
    nL = size(px)*padFactorSpace;
    dxL = (maxPx-minPx)/nL(2);
    dyL = (maxPy-minPy)/nL(1);
    [pxL pyL] = meshgrid(minPx:dxL:maxPx, minPy:dyL:maxPy);
    frame.pxL{1} = pxL;
    frame.pyL{1} = pyL;
    frameL = sumVel(frame,1,num);
    dt = 1/frame.freq/padFactorTime;
    frameL.dn = dn;
    kk = 0;
    tic
    for num = num
        if(toc>1) toc; end
        tic
        disp(num)
        kk = kk + 1;
        dr = 1e-4;
        tt = dn/frame.freq;
        nPhi = 10;

        [drf plf] = calcEval(frame, nPhi, dr, num+(0:dn-1), dt);
        [drb plb] = calcEval(frame, nPhi, dr, num+(0:-1:-dn+1), -dt);
        lf = -log(drf/abs(dr))./tt;
        lb = -log(drb/abs(dr))./tt;
%         size(plf)
%         size(lb)
        [lb lb_min lb_max] = double2uint8(lb);
        [lf lf_min lf_max] = double2uint8(lf);
        [plb plb_min plb_max] = double2uint8(plb);
        [plf plf_min plf_max] = double2uint8(plf);
       
        frameL.lyap_f{kk} = lf;
        frameL.lyap_b{kk} = lb;
        frameL.lyap_lf_min(kk) = lf_min;
        frameL.lyap_lf_max(kk) = lf_max;
        frameL.lyap_lb_min(kk) = lb_min;
        frameL.lyap_lb_max(kk) = lb_max;
        frameL.path_length_f{kk} = plf;
        frameL.path_length_f_min(kk) = plf_min;
        frameL.path_length_f_max(kk) = plf_max;
        frameL.path_length_b{kk} = plb;
        frameL.path_length_b_min(kk) = plb_min;
        frameL.path_length_b_max(kk) = plb_max;
%         end
    end
    toc
end

function [lb lb_min lb_max] = double2uint8(lb)
    lb(find(isinf(lb))) = NaN;
    lb_min = min(lb(:));
    lb_max = max(lb(:));
    lb = lb - min(lb(:));
    lb = round(lb/max(lb(:))*254);
    lb(find(isnan(lb))) = 255;
    lb = uint8(lb);
end


function [dr1 lr1 xx0 yy0 xx1 yy1] = calcEval(frame, nPhi, dr, num, dt)
    vx = cat(3, frame.vx{num});
    vy = cat(3, frame.vy{num});
    px = (frame.px{num(1)});
    py = (frame.py{num(1)});
    pxL = (frame.pxL{1});
    pyL = (frame.pyL{1});
    px1 = zeros([size(pxL) nPhi+1]);
    py1 = zeros([size(pxL) nPhi+1]);
    dphi = 2*pi/nPhi;
    kk = 1;
    px1(:,:,1) = pxL;
    py1(:,:,1) = pyL;
    for phi = [dphi:dphi:2*pi]
        kk = kk + 1;
        dx = real(dr*exp(i*phi));
        dy = imag(dr*exp(i*phi));
        px1(:,:,kk) = pxL + dx;
        py1(:,:,kk) = pyL + dy;
    end
    kk = 0;
    ttOrig = single(frame.tt(num));
    ttFull = ttOrig(1):dt:ttOrig(end);
    tt = repmat(ttOrig, [1, size(px)]);  
    tt = permute(tt, [2,3,1]);

    vx(find(isnan(vx))) = 0;
    vy(find(isnan(vy))) = 0; 
    lr1 = zeros(size(pxL));
    for ii = 1:numel(ttFull)
        kk = kk + 1; 
        ttCurr = ttFull(ii)*ones(size(pxL));
        num1 = findnear(ttOrig, ttCurr(1));
        vx1 = interp2(px, py, vx(:,:,num1), px1, py1, 'linear');
        vy1 = interp2(px, py, vy(:,:,num1), px1, py1, 'linear');
        lr1 = lr1 + dt*abs(vx1(:,:,1)+1i*vy1(:,:,1));
        px1 = px1 + vx1*dt;
        py1 = py1 + vy1*dt;
    end
    dr1 = zeros([size(pxL) nPhi]);
    for ii = 1:nPhi
        dx1 = px1(:,:,ii+1) - px1(:,:,1);
        dy1 = py1(:,:,ii+1) - py1(:,:,1);
        dr1(:,:,ii) = abs(dx1 + 1i*dy1);
    end
    dr1 = mean(dr1,3);
    
end

function [num alpha] = findTtNum(tt, ttCurr)
    indMin = max(find(tt<=ttCurr));
    indMax = min(find(tt>=ttCurr));
    num = min([indMin,indMax]:max([indMin,indMax]));
    if numel(num) == 1
        alpha = 1;
    else
        deltaAll = abs(diff(tt(num)));
        delta1 = abs(tt(num(1))-ttCurr);
        delta2 = abs(tt(num(2))-ttCurr);
        alpha(1,1,:) = [delta2/deltaAll delta1/deltaAll];
    end
end

