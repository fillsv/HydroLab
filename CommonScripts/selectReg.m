function dd_out = selectReg(dd, reg)
kk = 0;
for ii = 1:numel(dd)
    if ~isempty(strfind(dd(ii).name, reg))
        kk = kk + 1;
        dd_out(kk) = dd(ii);
    end
end
if exist('dd_out', 'var') == 0; dd_out = []; end    
