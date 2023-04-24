function [px1 py1] = moveVirtualParticles(frame, px1, py1, num, dt)
    vx = cat(3, frame.vx{num});
    vy = cat(3, frame.vy{num});
    px = (frame.px{num(1)});
    py = (frame.py{num(1)});
    px = px - mean(px(:));
    py = py - mean(py(:));
    ttOrig = single(frame.tt(num));
    ttFull = ttOrig(1)+dt/2:dt:ttOrig(end)-dt/2;
    tt = repmat(ttOrig, [1, size(px)]);  
    tt = permute(tt, [2,3,1]);
    vx(find(isnan(vx))) = 0;
    vy(find(isnan(vy))) = 0;    
    for ii = 1:numel(ttFull)
        ttCurr = ttFull(ii);%*ones(size(px1));
        [num1 alpha] = findTtNum(ttOrig, ttCurr(1));
%         if numel(alpha) == 48
%             num1 = num1;
%         end
%         alpha
        vx1 = interp2(px, py, sum(vx(:,:,num1).*alpha,3), px1, py1, 'linear');
        vy1 = interp2(px, py, sum(vy(:,:,num1).*alpha,3), px1, py1, 'linear');
        px1 = px1 + vx1*dt;
        py1 = py1 + vy1*dt;
    end
    
end

function [num alpha] = findTtNum(tt, ttCurr)
    if tt(2)>tt(1)
        indMin = max(find(tt<=ttCurr));
        indMax = min(find(tt>=ttCurr));
        num = indMin:indMax;
    else
        indMin = max(find(tt>=ttCurr));
        indMax = min(find(tt<=ttCurr));
        num = indMin:indMax;
        
    end
    if numel(num) == 1
        alpha = 1;
    else
        deltaAll = abs(diff(tt(num)));
        delta1 = abs(tt(num(1))-ttCurr);
        delta2 = abs(tt(num(2))-ttCurr);
        alpha(1,1,:) = [delta2/deltaAll delta1/deltaAll];
    end
    
end

