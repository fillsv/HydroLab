function dd_out = selectmV(dd, minmV, maxmV)
if exist('maxmV', 'var') == 0 maxmV = minmV; end;
kk = 0;
for ii = 1:numel(dd)
    mV = findmV(dd(ii).name);
    if ~isempty(mV)
        if (mV>=minmV)&&(mV<=maxmV)
            kk = kk + 1;
            dd_out(kk) = dd(ii);
        end
    end
end
if exist('dd_out', 'var') == 0; dd_out = []; end    
