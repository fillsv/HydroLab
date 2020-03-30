function createInfo(param)
    dirName = dir; 
    dirName = dirName(3:end); 
    dirName = dirName([dirName.isdir]);
    varNames = {'Lx', 'Ly', 'freq', 'sigma', 'rho', 'visc', 'h'};
    for kk = 1:numel(dirName)
%         disp(b(kk).name)
        cd(dirName(kk).name)
        pathm = 'mat/';
%         freq = 24;
%         Lx = 66;
%         Ly = 66;
        freq = 1012;
        Lx = 2.1;
        Ly = 2.1;
        tv_threshold = 0.14;
        eval(param)
        a = dir(pathm); a = a(3:end); a = a([a.isdir]);
        for ii = 1:numel(a)
            path1 = [a(ii).folder filesep a(ii).name filesep];
            disp(path1)
            load(nameMat(ii,1,pathm), 'x');
            for varName = varNames
                varName = varName{1};
                if exist(varName,'var')
                    eval(['matInfo.' varName ' = ' varName ';']);
                end
            end
%             matInfo.freq = freq;
%             matInfo.Lx = Lx;
%             matInfo.Ly = Ly;
            matInfo.tv_threshold = tv_threshold; 
            matInfo.tt = [];
            for jj = 1:numel(x)
                matInfo.tt(jj,1) = 180*(ii-1)+(jj-1)/freq;
            end
            save([path1 'matInfo.mat'], 'matInfo')
        end
        cd ..
    end
