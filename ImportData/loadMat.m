function frame = loadMat(num, nn, pathm, numres)
    if exist('pathm', 'var')
        if isempty(pathm)
            clear pathm
        end
    end


    if exist('pathm', 'var')
        nameMatFile = nameMat(num, nn, pathm);
        nameInfoFile = nameInfo(num, nn, pathm);
    else
        if exist('nn', 'var')
            nameMatFile = nameMat(num, nn);
            nameInfoFile = nameInfo(num, nn);
        else
            nameMatFile = nameMat(num);
            nameInfoFile = nameInfo(num);
        end
    end
    
    if nameMatFile == -1; frame = []; return; end
    if nameInfoFile == -1; frame = []; return; end
    load(nameMatFile)
    load(nameInfoFile)
    if ~exist('numres', 'var')
        numres = numel(results{1}.x);
    end
    frame = [];
    if exist('results', 'var')
        
        for ii = 1:numel(results)
%             disp(ii)
            res = filterResults(results{ii},10);
%             res = results{ii};
            x{ii} = res.x{numres};
            y{ii} = res.y{numres};
            u{ii} = res.u{numres};    
            v{ii} = res.v{numres};    
%             x{ii} = res.x{3};
%             y{ii} = res.y{3};
%             u{ii} = res.u{3};    
%             v{ii} = res.v{3};    
        end 
    end
    if exist('u', 'var')
        frame = importMat(x,y,u,v, matInfo);
    end
