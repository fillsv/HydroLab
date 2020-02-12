function frameO = calcVort(frame, num)
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end
    px = frame.px{1}(1:end-1, 1:end-1);
    py = frame.py{1}(1:end-1, 1:end-1);
    for ii = num
        [omega, cav] = curl(frame.px{ii}, frame.py{ii}, frame.vx{ii}, frame.vy{ii});

        omega(find(isnan(omega))) = 0;

        frameO.omega{ii} = omega;
        frameO.px{ii} = px;
        frameO.py{ii} = py;
        
%         mo(ii) = nanmean(omega(:).^2).^.5;
    end
    

    
