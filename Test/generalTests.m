genTestSignal
createInfo('Lx=66;Ly=66;freq=12;')
cd testsignal;
frame = loadMat(1);
cd ..

fprintf('Vort test')
% Ao = 1;
frO = calcVort(frame);
num = 1:numel(frame.vx);
omega = cat(3, frO.omega{num});

if abs(max(omega(:))-Ao) < .05*Ao
    fprintf(' passed!\n')
else
    fprintf(' failed!\n')
end
padFactor = 4;
frOF = calcVortFFT(frame, padFactor);
fomega = frOF.fomega;
passed = 1;
fprintf('VortFFT test')
for ii = 1:20;
    fprintf('.');
    nn = round(rand(1)*numel(num));
    [c1, kx1, ky1] = findMaxFFT2(frOF.fomega{nn}, frOF.kx, frOF.ky, frOF.padFactor, k, k);
    [c2, kx2, ky2] = findMaxFFT2(frOF.fomega{nn}, frOF.kx, frOF.ky, frOF.padFactor, k, -k);

    Aof = Ao/4;

    if ~((abs(abs(c1)-Aof)<0.05*Aof)&...
        (abs(abs(c2)-Aof)<0.05*Aof)&...
        (abs(abs(kx1)-k)<0.05*k)&...
        (abs(abs(ky1)-k)<0.05*k)&...
        (abs(abs(kx2)-k)<0.05*k)&...
        (abs(abs(ky2)-k)<0.05*k))
        passed = passed -1;
    end
end
if passed 
    fprintf(' passed!\n')
else
    fprintf(' %d failed!\n', abs(passed)+1)
end



padFactor = 2;
fprintf('Waves test')
frW = calcWaves(frame, padFactor, nyu);
passed = 1;

kxref = [k, -k, 0, 0];
kyref = [0, 0, k, -k];
A = [Ax/2 Ax/2 Ay/2 Ay/2];
for ii = 1:4
    fprintf('.');
     
    if ii <= 2
        f = frW.fft2vox;
    else
        f = frW.fft2voy;
    end

    [c, kx1, ky1] = findMaxFFT2(f, frW.kx, frW.ky, frW.padFactor,...
    kxref(ii)+rand(1)*1e-2*k, kyref(ii)+rand(1)*1e-2*k);
    if ~((abs(abs(c)-A(ii))<0.05*A(ii))&...
     (abs(abs(kx1-kxref(ii)))<0.05*k)&...
     (abs(abs(ky1-kyref(ii)))<0.05*k))
        passed = 0;
    end
end
if passed 
    fprintf(' passed!\n')
else
    fprintf(' failed!\n')
end


fprintf('WavesSpectra test')
frWS = calcWavesSpectra(frame, 2);
c = frWS.vsp(find(frWS.nyu==nyu));
cref = ((Ax/2)^2+(Ay/2)^2)^.5;
if (abs(c-cref)<0.05*cref)
     fprintf(' passed!\n')
else
    fprintf(' failed!\n')
end
fprintf('EnergySpectra test')

fr1 =  sumVel(frame, 12);
vx = cat(3, fr1.vx{:});
vy = cat(3, fr1.vx{:});
cref = mean(vx(:).^2+vy(:).^2)/2;
for padFactor = [1 4 8 16]
    fprintf('.');
    frES = calcEnergySpectra(calcEnergy(fr1, padFactor));
    dk = frES.k(2)-frES.k(1);
    c = sum(frES.Ek)*dk;
    if ~(abs(c-cref)<0.05*cref)
        passed = 0
    end
end
if passed 
    fprintf(' passed!\n')
else
    fprintf(' failed!\n')
end