function dd_out = selectHz(dd, minHz, maxHz)
if exist('maxHz', 'var') == 0 maxHz = minHz; end;
kk = 0;
for ii = 1:numel(dd)
    Hz = findHz(dd(ii).name);
    if ~isempty(Hz)
        if (Hz>=minHz)&&(Hz<=maxHz)
            kk = kk + 1;
            dd_out(kk) = dd(ii);
        end
    end
end
if exist('dd_out', 'var') == 0; dd_out = []; end    
