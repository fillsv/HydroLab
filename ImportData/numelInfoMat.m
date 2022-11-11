function ii = numelInfoMat(matNumber)
    if ~exist('matNumber', 'var')
        ii = 1;
        while nameInfo(ii)~=-1
            ii = ii + 1;
        end
        ii = ii - 1;
    else
        ii = 1;
        while nameInfo(matNumber, ii)~=-1
            ii = ii + 1;
        end
        ii = ii - 1;
    end
