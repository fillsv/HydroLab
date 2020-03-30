function frameVort = calcVort(frame, num)
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end
    px = frame.px{num(1)};
    py = frame.py{num(1)};
    for ii = 1:numel(num)
        [omega, cav] = curl(frame.px{num(ii)}, frame.py{num(ii)}, frame.vx{num(ii)}, frame.vy{num(ii)});

%         omega(find(isnan(omega))) = 0;

        frameVort.omega{ii} = omega;
        frameVort.px{ii} = px;
        frameVort.py{ii} = py;
%         mo(ii) = nanmean(omega(:).^2).^.5;
    end
    frameVort.tt(1:numel(num)) = frame.tt(num);
    frameVort.Lx = frame.Lx;    
    frameVort.Ly = frame.Ly; 
    frameVort.freq = frame.freq; 
    
