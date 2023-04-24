function frameL = calcLyap(frame, num, dn, padFactorSpace, padFactorTime)
    if ~exist('dn', 'var') dn = 25; end
    if ~exist('num', 'var') num = []; end
    if isempty(num) 
        num = dn:numel(frame.px)-dn+1
    end
    px = (frame.px{num(1)});
    py = (frame.py{num(1)});
    px = px - mean(px(1,:));
    py = py - mean(py(:,1));
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
    frameL.px = [];
    frameL.py = [];
    frameL.vx = [];
    frameL.vy = [];
    kk = 0;
    tic
    for num = num
        if kk > 0; toc; end
        tic
        disp(num)
        kk = kk + 1;
        dr = 1e-4;
        tt = dn/frame.freq;
        nPhi = 10;
        drf = calcEval(frame, nPhi, dr, num+(0:dn-1), dt);
        drb = calcEval(frame, nPhi, dr, num+(0:-1:-dn+1), -dt);
        lf = -log(drf/abs(dr))./tt;
        lb = -log(drb/abs(dr))./tt;
        frameL.lyap_f{kk} = single(lf);
        frameL.lyap_b{kk} = single(lb);
    end
    toc
end


function [dr1 xx0 yy0 xx1 yy1] = calcEval(frame, nPhi, dr, num, dt)
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
%     dt = 1/frame.freq;
    kk = 0;
    ttOrig = single(frame.tt(num));
    ttFull = ttOrig(1):dt:ttOrig(end);
    tt = repmat(ttOrig, [1, size(px)]);  
    tt = permute(tt, [2,3,1]);
%     px = cat(3, frame.px{num});
%     py = cat(3, frame.py{num});
%     vx = (frame.vx{ii});
%     vy = (frame.vy{ii});
    vx(find(isnan(vx))) = 0;
    vy(find(isnan(vy))) = 0;   
    [px1, py1] = moveVirtualParticles(frame, px1, py1, num, dt);
%     for ii = 1:numel(ttFull)
% %         ii
%         kk = kk + 1; 
%         ttCurr = ttFull(ii)*ones(size(pxL));
% %         [num1 alpha] = findTtNum(ttOrig, ttCurr(1));
% %         vx1 = interp2(px, py, sum(vx(:,:,num1).*alpha,3), px1, py1, 'linear');
% %         vy1 = interp2(px, py, sum(vy(:,:,num1).*alpha,3), px1, py1, 'linear');
%         num1 = findnear(ttOrig, ttCurr(1));
%         vx1 = interp2(px, py, vx(:,:,num1), px1, py1, 'linear');
%         vy1 = interp2(px, py, vy(:,:,num1), px1, py1, 'linear');
%         
% %         if numel(num1) == 1
% %             vx1 = interp2(px(:,:,num1), py(:,:,num1), vx(:,:,num1), px1, py1, 'linear');
% %             vy1 = interp2(px(:,:,num1), py(:,:,num1), vy(:,:,num1), px1, py1, 'linear');
% % %             vx0 = interp2(px(:,:,num1), py(:,:,num1), vx(:,:,num1), px0, py0, 'linear');
% % %             vy0 = interp2(px(:,:,num1), py(:,:,num1), vy(:,:,num1), px0, py0, 'linear');
% %         else
% %             vx1 = interp3(px(:,:,num1), py(:,:,num1), tt(:,:,num1), vx(:,:,num1), px1, py1, ttCurr, 'linear');
% %             vy1 = interp3(px(:,:,num1), py(:,:,num1), tt(:,:,num1), vy(:,:,num1), px1, py1, ttCurr, 'linear');
% %             vx0 = interp3(px(:,:,num1), py(:,:,num1), tt(:,:,num1), vx(:,:,num1), px0, py0, ttCurr, 'linear');
% %             vy0 = interp3(px(:,:,num1), py(:,:,num1), tt(:,:,num1), vy(:,:,num1), px0, py0, ttCurr, 'linear');
% %         end
%         px1 = px1 + vx1*dt;
%         py1 = py1 + vy1*dt;
% %         px0 = px0 + vx0*dt;
% %         py0 = py0 + vy0*dt;
% %         dx1 = px1-px0;
% %         dy1 = py1-py0;
% %         dr1(kk,:,:) = abs(dx1 + 1i*dy1);
% %         xx0(kk,:,:) = px0;
% %         yy0(kk,:,:) = py0;
% %         xx1(kk,:,:) = px1;
% %         yy1(kk,:,:) = py1;
%     end
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

% function [dr1 xx0 yy0 xx1 yy1] = calcEval(frame, dr, num, dt)
%     dx = real(dr);
%     dy = imag(dr);
%     px = (frame.px{num(1)});
%     py = (frame.py{num(1)});
%     pxL = (frame.pxL{1});
%     pyL = (frame.pyL{1});
%     px1 = pxL + dx;
%     py1 = pyL + dy;
%     px0 = pxL;
%     py0 = pyL;
% %     dt = 1/frame.freq;
%     kk = 0;
%     for ii = num
%         kk = kk + 1; 
%         vx = (frame.vx{ii});
%         vy = (frame.vy{ii});
%         vx(find(isnan(vx))) = 0;
%         vy(find(isnan(vy))) = 0;
%         vx1 = interp2(px, py, vx, px1, py1);
%         vy1 = interp2(px, py, vy, px1, py1);
%         px1 = px1 + vx1*dt;
%         py1 = py1 + vy1*dt;
%         vx0 = interp2(px, py, vx, px0, py0);
%         vy0 = interp2(px, py, vy, px0, py0);
%         px0 = px0 + vx0*dt;
%         py0 = py0 + vy0*dt;
%         dx1 = px1-px0;
%         dy1 = py1-py0;
% %         dr1(kk,:,:) = abs(dx1 + 1i*dy1);
% %         xx0(kk,:,:) = px0;
% %         yy0(kk,:,:) = py0;
% %         xx1(kk,:,:) = px1;
% %         yy1(kk,:,:) = py1;
%     end
%     dr1 = abs(dx1 + 1i*dy1);
%     
% end
