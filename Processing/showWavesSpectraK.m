function [vsp nyu] = showWavesSpectraK(frame, padFactor, num)
    %
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
%     vsp = abs(frame.vsp./frame.coef);
    vsp = frame.vsp./frame.coef;
    nyu = frame.nyu;
%     k = disper(frame.nyu, frame.sigma);
%     coef = (frame.n-1)/frame.n*frame.h*k;
    plot(frame.nyu, vsp);
%     plot(frame.nyu, frame.vsp./disper(frame.nyu)', 'b');
    xlabel('Frequency, Hz');
    ylabel('RMS velocity, cm/s');

    xlim([min(frame.nyu) max(frame.nyu)])
    ylim([min(vsp)/1.1 max(vsp)*1.1])
    set(gcf, 'PaperPosition', [0 0 17 10]);
end