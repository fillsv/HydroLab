function frame = sumVel(frame, numstep, num)
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end
    px = frame.px{1};
    py = frame.py{1};
    numNewFrame = floor(numel(num)/numstep);
    vx = cat(3, frame.vx{num(1:numNewFrame*numstep)});
    vy = cat(3, frame.vy{num(1:numNewFrame*numstep)});
    frame.tt = mean(reshape(frame.tt(num(1:numNewFrame*numstep)), [numstep, numNewFrame]),1)';
    vx = reshape(vx, [size(vx,1) size(vx,2), numstep, numNewFrame]);
    vy = reshape(vy, [size(vx,1) size(vx,2), numstep, numNewFrame]);
    vx = squeeze(nanmean(vx, 3));
    vy = squeeze(nanmean(vy, 3));
    
    frame.vx = {};
    frame.vy = {};
    frame.px = {};
    frame.py = {};    
    for ii = 1:numNewFrame
        frame.vx{ii,1} = vx(:,:,ii);
        frame.vy{ii,1} = vy(:,:,ii);
        frame.px{ii,1} = px;
        frame.py{ii,1} = py;
        
    end
    frame.freq = frame.freq/numstep;
end