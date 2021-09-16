function showTensors(frame, px0, py0, num)

    if ~exist('frame', 'var')
        disp('function showTensors(frame, px0, py0, num)');
        return
    end
    
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end    
    px = cat(3, frame.px{num});
    py = cat(3, frame.py{num});
    vx = cat(3, frame.vx{num});
    vy = cat(3, frame.vy{num});
    dx = px(1,2)-px(1,1);
    dy = py(2,1)-py(1,1);    
    mult = 20;
    [px1 py1] = meshgrid(px(1):dx/mult:px(end), py(1):dy/mult:py(end));
    vx = interp2(px, py, vx, px1, py1);
%     fprintf('.')
    vy = interp2(px, py, vy, px1, py1);
%     fprintf('.')
%     sig = 10;
%     vx = imgaussfilt(vx, sig, 'padding', 'symmetric');
%     vy = imgaussfilt(vy, sig, 'padding', 'symmetric');
    px = px1;
    py = py1;

   
    px = px - mean(px, 'all');
    py = py - mean(py, 'all');
    rx = px-px0;
    ry = py-py0;     

    [omega, cav] = curl(px, py, vx, vy);
    v_phi = (-ry.*vx+rx.*vy)./abs(rx+i*ry);
    v_r = (rx.*vx+ry.*vy)./abs(rx+i*ry);
    
    r = abs(rx+i*ry);

    dr = frame.Lx/250;
    r0 = r;
    [r r_ind] = sort(r(:));
    v_phi = v_phi(r_ind);
    r0 = r0(r_ind);
    v_r = v_r(r_ind);
    omega = omega(r_ind);


    Lx = frame.Lx;
    Ly = frame.Ly;

    numk = double(floor(max(r(:))/dr));
%     ar = ((1:numk)-.5)*dr;
 
    u = zeros(1,numk);
    ar = zeros(1,numk);
    ar1 = zeros(1,numk);
    num_step = 2*numel(r)/((numk+1)^2+numk+1);
% 
    for ii = 1:numk
%         ind = r_ind(round(num_step*ii*(ii-1)/2+1):round(num_step*ii*(ii+1)/2));
        ind = (round(num_step*ii*(ii-1)/2+1):round(num_step*ii*(ii+1)/2));

        u(ii) = nanmean(v_phi(ind));
        ar(ii) = nanmean(r0(ind));
        ov(ii) = nanmean(omega(ind).*v_r(ind));
        o(ii) = nanmean(omega(ind));
        st(ii) = nanmean(v_phi(ind).*v_r(ind));
        
    end
    do = diff(o)./diff(ar);
    arn = (ar(1:end-1)+ar(2:end))/2;

    ur = u.*ar;
    dur = diff(ur)./diff(ar);
    omega = dur./arn;
   
    ur = u./ar;
    dur = diff(ur)./diff(ar);
    sigma = dur.*arn;

%     plot(ar, u);
%     plot(ar, st, arn, sigma*.01)
    yyaxis left
    plot(ar, st, arn, sigma*.01)
    ylim([-1 1]*1.1*max(abs([st, sigma*.01])));
    yyaxis right
    plot(ar, u)
    ylim([-1 1]*1.1*max(abs(u)));
    legend('Reynolds', '\Sigma*\nu', 'u_\phi')
    
%     plot(ar, st)
%       plot(arn, omega, ar, ov)
%       legend('\Omega', '\Omega*u_r')
%      plot(arn, omega)%, ar, o);
%      ind = find(ov<0);
%      plot(ar, ov,'.b', ar(ind), ov(ind), '.r')%, ar, dmo);
%      plot(ar, ov,'.b', ar(ind), ov(ind), '.r')%, ar, dmo);
     
%     ylim([-1 1]*3e-4)
     xlabel('Radius, cm');
%     ylabel('Omega, 1/s');
%     ylim([-1 1]*1);
    
%--------------------------    
    
%    
%     plot((ar(1:end-1)+.25), sigma);
% %     plot(ar, u);
%     xlabel('Radius, cm');
%     ylabel('Sigma, 1/s');

    
%     nyu = 1e-2;
%     
% %     ur = mov.*ar;
% %     dur = (ur(2:end)-ur(1:end-1))/dr;
% %     
% %     omega = dur./(ar(1:end-1)+.25);
% %    
% %     plot((ar(1:end-1)+.25), omega);
%     plot(ar, mov)%, ar, dmo);
%     xlabel('Radius, cm');
%     ylabel('<\Omega U_r>, cm/s^2');
