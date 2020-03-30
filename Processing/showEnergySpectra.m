function showEnergySpectr(frame, padFactor, num)
    %
    if exist('padFactor', 'var')==0 
        padFactor = 4;
    end 
    if ~isfield(frame, 'Ek')
        if isfield(frame, 'vx')
            if exist('num', 'var')==0 
                num = 1:numel(frame.vx);
            end 
            frame = calcEnergy(frame, padFactor, num);
            num = 1:numel(frame.E);

        end
        if isfield(frame, 'E')
            if exist('num', 'var')==0 
                num = 1:numel(frame.E);
            end 
            frame = calcEnergySpectra(frame, num)
        end
    end    
    if (~isfield(frame, 'Ek'))||(~isfield(frame, 'k'))
        warning('frame is incorrect!')
        return
    end
    
    p = loglog(frame.k, frame.Ek, 'b');
    xlabel('Wave vector, cm^{-1}');
    ylabel('Energy spectrum, cm^3/s^2');

    xlim([min(frame.k) max(frame.k)])
    ylim([min(frame.Ek)/2 max(frame.Ek)*2])
    set(gcf, 'PaperPosition', [0 0 17 10]);
end