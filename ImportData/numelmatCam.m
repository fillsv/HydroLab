function numela = numelmat(num, nn, pathm)
    if(exist('nn', 'var') == 0) nn = 1; end;
    if(exist('pathm', 'var') == 0) pathm = 'mat/'; end;
%    a = a(3:end);
    strnum = num2str(num, '%02d_cam*');
    a = dir([pathm strnum]);

    a = a(find([a.isdir]));
    numela = numel(a);
    
    numelb = 0;
    if numela > 0 & exist('num', 'var')
        path1 = [a(1).name '/'];
        b = dir([pathm path1 'MVI*.mat']); 
        numela = numel(b);
        b = dir([pathm path1 'cam*.mat']); 
        numela = numela + numel(b);
    end
    
end