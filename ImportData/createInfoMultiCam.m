function createInfoMultiCam(param)
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
%    a = dir(pathm); a = a(3:end); a = a([a.isdir]);
    if ~exist('num', 'var') num = 1:numel(a); end
    for ii = num
        strnum = num2str(ii, '%02d_cam*');
        a = dir([pathm strnum]);
        a = a(find([a.isdir]));

        pp = 0;

        path1 = [a(1).folder filesep a(1).name filesep];
        disp(path1)
        [ nb] = numelmatCam(ii);
        for kk = 1:nb
            disp(kk)
            if nameMatCam(ii,kk,pathm) == -1 
                continue;
            end
            s = load(nameMatCam(ii,kk,pathm), 'x');
            if ~isfield(s, 'x')
                fileInfo = whos('-file', nameMatCam(ii,kk,pathm));
            else
                fileInfo(1).size = size(s.x);
            end
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
            for jj = 1:max(fileInfo(1).size)
                pp = pp + 1;
%                 matInfo.tt(jj,1) = 180*(ii-1)+(pp-1)/freq;
                matInfo.tt(jj,1) = (pp-1)/freq;
            end
            tmp_name = nameMatCam(ii, kk);
            ind = find(tmp_name == filesep);
            tmp_name = tmp_name(ind(end)+1:end);
            nameout = ['matInfo_' tmp_name];
            save([path1 nameout], 'matInfo')
        end
    end
