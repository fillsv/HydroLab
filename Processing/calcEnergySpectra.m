function frameEnergySpectra = calcEnergySpectra(frame, num)
    %
    if exist('num', 'var')==0 
        num = 1:numel(frame.E);
    end 
%     num
    E = cat(3, frame.E{num(:)});
    E = mean(E, 3);
    
    kx = frame.kx;
    ky = frame.ky;
    
    kmax = max(kx(:));

    kr = (kx.^2+ky.^2).^.5;
    [k_temp, k_ind] = sort(kr(:));
    Lx = frame.Lx;
    Ly = frame.Ly;
    kstep = 2/Lx*2/sqrt(2)/2;
    dkx = kx(1,2)-kx(1,1);
    dky = ky(2,1)-ky(1,1);
    numk = floor(kmax/kstep);
%     k = ((1:numk)-.5)*kstep;
    Ek = zeros(numk,1);
    num_step = 2*pi*numel(kr)/(numk^2+numk)/4;

    for ii = 1:numk
        in = k_ind(round(num_step*ii*(ii-1)/2+1):round(num_step*ii*(ii+1)/2));
        k(ii) = mean(kr(in));

        Ek(ii) = 1/2/Lx/Ly/kstep*sum(abs(2*E(in)))*dkx*dky; % intergal!
    
    end
    frameEnergySpectra.Ek = Ek;
    frameEnergySpectra.k = k';
    frameEnergySpectra.Lx = frame.Lx;
    frameEnergySpectra.Ly = frame.Ly;
    frameEnergySpectra.tt = mean(frame.tt(num));
% 
%     if ~exist('name', 'var') name = 1; end
%     if name(1)~=0
%         p = loglog(k, Ek, 'b');
%         xlabel('Wave vector, cm^{-1}');
%         ylabel('Energy spectrum, cm^3/s^2');
% 
%         xlim([min(k) max(k)])
%         % xlim([0.05 max(k)])
%         ylim([min(Ek)/2 max(Ek)*2])
%         set(gcf, 'PaperPosition', [0 0 17 10]);
%     end
end