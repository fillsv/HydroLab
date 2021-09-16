function name = nameInfo(num, nn, pathm)
if(exist('nn', 'var') == 0) nn = 1; end;
if(exist('pathm', 'var') == 0) pathm = 'mat/'; end;
a = dir(pathm);
a = a(3:end);
a = a(find([a.isdir]));
name = -1;
if(numel(a)>=num)
    path1 = [a(num).name '/'];
    b = dir([pathm path1 'matInfo*.mat']); 
    if numel(b)
        name = [pathm path1 b(nn).name];
    end
end
