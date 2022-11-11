function frame = loadFragment(expNumber, tBeg, deltaT)
% tBeg = 35;
tEnd = tBeg+deltaT;
% expNumber = 3;
num = numelInfoMat(expNumber);
tt = {};
tt_all = [];
for ii = 1:num
%     disp(ii)
    load(nameInfo(expNumber, ii))
    tt{ii} = matInfo.tt;
    if size(matInfo.tt,2) == 1
        matInfo.tt = matInfo.tt';
    end
    tt_all = [tt_all matInfo.tt];
end
ttBeg = min(tt_all(find(tt_all>=tBeg)));
ttEnd = max(tt_all(find(tt_all<=tEnd)));
for ii = 1:num
    if ~isempty(find(tt{ii}==ttBeg))
        matBeg = ii;
        iiBeg = find(tt{ii}==ttBeg);
        jjBeg = find(tt_all==ttBeg);
    end
    if ~isempty(find(tt{ii}==ttEnd))
        matEnd = ii;
        iiEnd = find(tt{ii}==ttEnd);
        jjEnd = find(tt_all==ttEnd);
    end
end
frame = [];

for ii = matBeg:matEnd
    frame = concatFrame(frame,sumVel(loadMat(expNumber,ii), 1));
end
frame = sumVel(frame,1, iiBeg+(0:jjEnd-jjBeg));
    