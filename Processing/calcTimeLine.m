clear tl
kk = 0;
timeStep = .2; %s
timeWindow = 1; %s
for ii = 1:numelInfoMat()
    disp(ii)
    frame = loadMat(ii);
    frameOmega = calcVort(frame);
    for num = genFrameNum(frame, timeStep, timeWindow);
        kk = kk + 1;
        num = num{1};
        omega = cat(3, frameOmega.omega{num});
        tl.mo(kk) = nanmean(abs(omega(:)));
        tl.tt(kk) = mean(frame.tt(num));
    end
    loglog(tl.tt, tl.mo)
    xlabel('Time, s');
    drawnow;
end
save('tl', 'tl')
