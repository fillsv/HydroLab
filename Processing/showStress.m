function showStress(frame, px0, py0, num)

    if ~exist('frame', 'var')
        disp('function showStream(frame, num, step, distance)');
        return
    end
    
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end    
    px = cat(3, frame.px{num});
    py = cat(3, frame.py{num});
    vx = cat(3, frame.vx{num});
    vy = cat(3, frame.vy{num});
    px = px - mean(px, 'all');
    py = py - mean(py, 'all');
    
    rx = px-px0;
    ry = py-py0;     
    v_phi = (-ry.*vx+rx.*vy)./abs(rx+i*ry);
    v_r = (rx.*vx+ry.*vy)./abs(rx+i*ry);
    r = abs(rx+i*ry);
    dr = .5;
    [r ind] = sort(r(:));
    v_phi = v_phi(ind);
    v_r = v_r(ind);
    for ii = 1:floor(max(r(:))/dr)
        ind = find((r>=(ii-1)*dr)&(r<ii*dr));
        st(ii) = nanmean(v_phi(ind).*v_r(ind));
        ar(ii) = (ii-.5)*dr;
    end
   
    plot(ar, st);
    xlabel('Radius, cm');
    ylabel('Reynolds tensor, cm^2/s^2');
