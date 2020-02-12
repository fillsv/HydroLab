function num = genNumFrame(frame, timeStep, timeWindow)
    numelFrame = numel(frame.px);
    dt = 1/frame.freq;
    numStep = round(timeStep/dt);
    numWindow = round(timeWindow/dt);
    numFrame = floor((numelFrame-numWindow)/numStep)+1;
    kk = 0;
    jj = 0;
    while numelFrame >= jj + numWindow + numStep
%         disp(jj)
        kk = kk + 1;
        num{kk} = jj + (1:numWindow);
        jj = jj + numStep;
    end
