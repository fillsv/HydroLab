function [xx yy] = showTrackSpot(frame, px0, py0, num, dt, rMax, numPoints, color)
    px = double(frame.px{num(1)});
    py = double(frame.py{num(1)});
%     px0 = px0 + mean(px(1,:));
%     py0 = py0 + mean(py(:,1));    
%     rMax = .3;
    dr = rMax*sqrt(pi/numPoints);
    [px1 py1] = meshgrid(px0 + (-rMax:dr:rMax),...
        py0 + (-rMax:dr:rMax));
    ind = find(abs(px1-px0+1i*(py1-py0))<rMax);
    px1 = px1(ind);
    py1 = py1(ind);
    kk = 1;
    [px1, py1] = moveVirtualParticles(frame, px1, py1, num, dt);
%     xx(kk,:) = px1 - mean(px(1,:));
%     yy(kk,:) = py1 - mean(py(:,1));
%     for ii = num
%         kk = kk + 1; 
%         vx = double(frame.vx{ii});
%         vy = double(frame.vy{ii});
%         vx(find(isnan(vx))) = 0;
%         vy(find(isnan(vy))) = 0;
%         vx1 = interp2(px, py, vx, px1, py1);
%         vy1 = interp2(px, py, vy, px1, py1);
%         if ii == num(end)
%             continue;
%         end
%         px1 = px1 + vx1*dt;
%         py1 = py1 + vy1*dt;
% %         xx(kk,:) = px1 - mean(px(1,:));
% %         yy(kk,:) = py1 - mean(py(:,1));
%     end
    xx = px1;% - mean(px(1,:));
    yy = py1;% - mean(py(:,1));
    zz = 100*ones(size(xx));
    hold on
    s = plot3(xx, yy, zz, ['.' color],...
        'MarkerSize', 3);
%     alpha(s, 1e-10)
    hold off
end
