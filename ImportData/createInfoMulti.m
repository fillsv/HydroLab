function createInfoMulti(param)
    dirName = dir; 
    dirName = dirName(3:end); 
    dirName = dirName([dirName.isdir]);
    varNames = {'Lx', 'Ly', 'freq', 'sigma', 'rho', 'visc', 'h', 'n'};
%     for kk = 1:numel(dirName)
%         disp(b(kk).name)
%         cd(dirName(kk).name)
    pathm = 'mat/';
%         freq = 24;
%         Lx = 66;
%         Ly = 66;
    freq = 400;
    Lx = 7;
    Ly = 7;
    n = 1.33;
%     h = 1;
    sigma = 73;
    tv_threshold = 0.14;
    if exist('param', 'var')
        eval(param)
    end
    a = dir(pathm); a = a(3:end); a = a([a.isdir]);
    for ii = 1:numel(a)
        pp = 0;

        path1 = [a(ii).folder filesep a(ii).name filesep];
        disp(path1)
        [ nb] = numelmat(ii);
        for kk = 1:nb
            if nameMat(ii,kk,pathm) == -1 
                continue;
            end
            load(nameMat(ii,kk,pathm), 'x');
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
                pp = pp + 1;
                matInfo.tt(jj,1) = 180*(ii-1)+(pp-1)/freq;
            end
            tmp_name = namemat(ii, kk);
            ind = find(tmp_name == filesep);
            tmp_name = tmp_name(ind(end)+1:end);
            nameout = ['matInfo_' tmp_name];
            save([path1 nameout], 'matInfo')
        end
    end
