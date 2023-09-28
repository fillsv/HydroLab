function [name, nameDir] = nameMatCam(num, nn, pathm)
    if(exist('nn', 'var') == 0) nn = 1; end;
    if(exist('pathm', 'var') == 0) pathm = 'mat/'; end;
    strnum = num2str(num, '%02d_cam*');
    if pathm(end)~=filesep pathm(end+1) = filesep; end
    a = dir([pathm strnum]);
%     a = a(3:end);
    a = a(find([a.isdir]));
    name = -1;
    nameDir = -1;
    if(numel(a)>=1)
        path1 = [a(1).name filesep];
% 
%         b = dir([pathm path1 'MVI*.mat']); 
%         if isempty(b)
        b = dir([pathm path1 'cam*.mat']); 
%         end
%         if isempty(b)
%             b = dir([pathm path1 'tif*.mat']); 
%         end
        if(numel(b)>=nn)
            name = [pathm path1 b(nn).name];
            nameDir = path1;
        end
    end
end
