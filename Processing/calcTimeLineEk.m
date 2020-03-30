clear tl
cd /media/fill/Exp/bunkerlqc/mat/Exppool/BigVortexDepth
b = dir(); b = b(3:end); b = b([b.isdir]);
b  = selectNReg(b, '10cm');
b  = selectNReg(b, 'DeepWater_19cm');

kk = 0;
for oo = 1:numel(b)
    cd(b(oo).name)
    a = dir(); a = a(3:end); a = a([a.isdir]);
 %   createInfo;
    
    timeStep = 20; %s
    timeWindow = 20; %s
    for jj = 3;%1:numel(a)
        cd(a(jj).name)
        kk = kk + 1;
        ll = 0;
        for ii = 1:9;%numelInfoMat()
            disp([oo jj ii])
            frame = loadMat(ii);
            frame = sumVel(frame, 240);
            for num = genFrameNum(frame, timeStep, timeWindow);
                fprintf('.')
                ll = ll + 1;
                num = num{1};
%                 omega = cat(3, frameOmega.omega{num});
%                 tl.mo(kk) = nanmean(abs(omega(:)));
                tl.tt(kk, ll) = mean(frame.tt(num));
%                   num = numel(frame.tt)-50:numel(frame.tt);
                frEk = calcEnergySpectr(calcEnergy(frame, 4, num));
                tl.Ek(kk, ll, :) = frEk.Ek;
                tl.k(kk, ll, :) = frEk.k;
                tl.Ebig(kk, ll) = mean(frEk.Ek(find(frEk.k<.2)));
                tl.kmax(kk, ll) = frEk.k(find(max(frEk.Ek)==frEk.Ek));
                tl.Epump(kk, ll) = mean(frEk.Ek(find((frEk.k>.4)&(frEk.k<2))));
            end            
            fprintf('\n')

            tl.mV(kk) = findmV(cd);
            tl.depth(kk) = findDepth(cd);
            plot(tl.tt(kk,:), tl.Ebig(kk,:), '.-')
            xlabel('Time, s');
%             xlabel('Depth, cm');
            drawnow;
        end
        cd ..
    end
    cd ..
end
save('tl', 'tl')
