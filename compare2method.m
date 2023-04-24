freq = 6;
kk = 0;
for t1 = 1:150
    kk = kk + 1;
    num = 1:250+t1*12;
    fr = calcWaves(frame, 16, freq, num);
    fr1 = calcWaves(frame1, 16, freq, num);
    kxm = findMax(fr);
    kxm1 = findMax(fr1);
    rel(kk)=kxm1./kxm;
    tt(kk) = t1;
    plot(tt, rel)
    drawnow;
end

function kxm = findMax(fr)
    ind1 = find(fr.kx>0);
    ind = find(max(abs(fr.fft2vox(ind1)))==abs(fr.fft2vox(ind1)));
    ind = ind1(ind(1));
    [c, kx_p, ky_p] = findMaxFFT2(fr.fft2vox, fr.kx, fr.ky, fr.padFactor, fr.kx(ind), 0);
    kxm = kx_p;
end