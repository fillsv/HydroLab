function frameDiver = calcDiver(frame, num)
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end
    px = frame.px{num(1)};
    py = frame.py{num(1)};
    for ii = 1:numel(num)
        [diver] = divergence(frame.px{num(ii)}, frame.py{num(ii)}, frame.vx{num(ii)}, frame.vy{num(ii)});

%         omega(find(isnan(omega))) = 0;

        frameDiver.diver{ii} = diver;
        frameDiver.px{ii} = px;
        frameDiver.py{ii} = py;
%         mo(ii) = nanmean(omega(:).^2).^.5;
    end
    frameDiver.tt(1:numel(num)) = frame.tt(num);
    frameDiver.Lx = frame.Lx;    
    frameDiver.Ly = frame.Ly; 
    frameDiver.freq = frame.freq; 
    
