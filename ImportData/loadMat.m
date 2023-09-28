function frame = loadMat(num, nn, pathm, numres)
    if exist('pathm', 'var')
        if isempty(pathm)
            clear pathm
        end
    end


    if exist('pathm', 'var')
        nameMatFile = nameMatCam(num, nn, pathm);
        nameInfoFile = nameInfoCam(num, nn, pathm);
    else
        if exist('nn', 'var')
            nameMatFile = nameMatCam(num, nn);
            nameInfoFile = nameInfoCam(num, nn);
        else
            nameMatFile = nameMatCam(num);
            nameInfoFile = nameInfoCam(num);
        end
    end
    
    if nameMatFile == -1; frame = []; return; end
    if nameInfoFile == -1; frame = []; return; end
    load(nameMatFile)
    load(nameInfoFile)
    if exist('results', 'var')
        if ~exist('numres', 'var')
            numres = numel(results{1}.x);
        end
        frame = [];
        if exist('results', 'var')
            for ii = 1:numel(results)
%                 res = filterResults(results{ii},10);
                res = results{ii};
                x{ii} = res.x{numres};
                y{ii} = res.y{numres};
                u{ii} = res.u{numres};    
                v{ii} = res.v{numres};    
            end 
        end
    end    
    if exist('u', 'var')
        frame = importMat(x,y,u,v, matInfo);
    end
