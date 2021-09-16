function openFolder(mask1, mask2, mask3, mask4)
    if ~exist('mask1', 'var') 
        [str_temp rez_temp] = system(['nautilus -w ' cd]);
        return;
    else
        mask1 = ['*' mask1 '*'];
    end
    if ~exist('mask2', 'var') 
        mask2 = '';
    else
        mask2 = [mask2 '*'];
    end
    if ~exist('mask3', 'var') 
        mask3 = '';
    else
        mask3 = [mask3 '*'];
    end
    if ~exist('mask4', 'var') 
        mask4 = '';
    else
        mask4 = [mask4 '*'];
    end
    a = dir([mask1 mask2 mask3 mask4]);
    a = a(find([a.isdir]));
%     a
    for ii = 1:numel(a)
        fprintf([a(ii).name '\n'])
    end
    if numel(a) == 1
        [str_temp rez_temp] = system(['nautilus -w ' cd '/' a(1).name]);
    end
end