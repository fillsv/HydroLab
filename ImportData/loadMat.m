function frame = loadMat(num, nn, pathm)
    if exist('pathm', 'var')
        nameMatFile = namemat(num, nn, pathm);
        nameInfoFile = nameInfo(num, pathm);
    else
        if exist('nn', 'var')
            nameMatFile = namemat(num, nn);
            nameInfoFile = nameInfo(num);
        else
            nameMatFile = namemat(num);
            nameInfoFile = nameInfo(num);
        end
    end
    
    if nameMatFile == -1; frame = []; return; end
    if nameInfoFile == -1; frame = []; return; end
    load(nameMatFile)
    load(nameInfoFile)

    frame = importMat(x,y,u,v, matInfo);