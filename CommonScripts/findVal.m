function value = findVal(template, str)
    if ~exist('str', 'var') str = cd; end
    inde = strfind(str, template)-1;
    str1 = str;
    str1(find(str1 == 'p')) = '0';
    indb = inde;
    for oo = inde:-1:1
        if isempty(str2num(str1(oo)))
            break;
        end
        indb = oo;        
    end
    str(find(str == 'p')) = '.';
    
    value = str2num(str(indb:inde));
end