function freq = findGl(str)
    inde = strfind(str, 'Gl')+2;
    str1 = str;
    str1(find(str1 == 'p')) = '0';
    indb = inde;
    for oo = indb:numel(str)
        if isempty(str2num(str1(oo)))
            break;
        end
        inde = oo;        
    end
    str(find(str == 'p')) = '.';
    
    freq = str2num(str(indb:inde));
end