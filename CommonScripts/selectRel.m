function dd_out = selectRel(dd, minRel, maxRel)
%     dd_out = [];
if exist('maxRel', 'var') == 0 maxRel = minRel; end;
kk = 0;
for ii = 1:numel(dd)
    rel = findRel(dd(ii).name);
    if ~isempty(rel)
        if (rel>=minRel)&&(rel<=maxRel)
            kk = kk + 1;
            dd_out(kk) = dd(ii);
        end
    end
end
if exist('dd_out', 'var') == 0; dd_out = []; end    

