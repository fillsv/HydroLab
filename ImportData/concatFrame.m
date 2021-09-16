function frame1 = concatFrame(frame1, frame2)
    if ~exist('frame1', 'var') || ~exist('frame2', 'var')
        disp('frame = concatFrame(frame1, frame2)');
        return
    end
    for ii = 1:numel(frame2.px)
        frame1.px{end+1,1} = frame2.px{ii};
        frame1.py{end+1,1} = frame2.py{ii};
        frame1.vx{end+1,1} = frame2.vx{ii};
        frame1.vy{end+1,1} = frame2.vy{ii};
        frame1.tt(end+1,1) = frame2.tt(ii);
    end
end