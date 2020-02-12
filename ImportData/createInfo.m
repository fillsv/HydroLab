function createInfo
    b = dir; b = b(3:end); b = b([b.isdir]);

    for kk = 1:numel(b)
%         disp(b(kk).name)
        cd(b(kk).name)
        pathm = 'mat/';
        freq = 24;
        Lx = 66;
        Ly = 66;
        a = dir(pathm); a = a(3:end); a = a([a.isdir]);
        for ii = 1:numel(a)
            path1 = [a(ii).folder filesep a(ii).name filesep];
            disp(path1)
            load(namemat(1), 'x');
            matInfo.freq = freq;
            matInfo.Lx = Lx;
            matInfo.Ly = Ly;
            for jj = 1:numel(x)
                matInfo.tt(jj,1) = 180*(ii-1)+(jj-1)/freq;
            end
            save([path1 'matInfo.mat'], 'matInfo')
        end
        cd ..
    end
