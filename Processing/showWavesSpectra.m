function [vsp nyu] = showWavesSpectra(frame, padFactor, num)
    %
    if ~exist('frame', 'var')
        disp('function showWavesSpectra(frame, padFactor, num)');
        return
    end    
    if exist('padFactor', 'var')==0 
        padFactor = 4;
    end 
    if ~isfield(frame, 'vsp')
        if isfield(frame, 'vx')
            if exist('num', 'var')==0 
                num = 1:numel(frame.vx);
            end 
            frame = calcWavesSpectra(frame, padFactor, num);
%             num = 1:numel(frame.E);

        end
    end    
    vsp = frame.vsp;
    nyu = frame.nyu;
    semilogy(frame.nyu, frame.vsp);
%     plot(frame.nyu, frame.vsp./disper(frame.nyu)', 'b');
    xlabel('Frequency, Hz');
    ylabel('RMS velocity, cm/s');

    xlim([min(frame.nyu) max(frame.nyu)])
    ylim([min(frame.vsp)/1.1 max(frame.vsp)*1.1])
    set(gcf, 'PaperPosition', [0 0 17 10]);
end