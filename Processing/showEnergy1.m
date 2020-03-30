function showEnergy1(frameEnergy, maxk, padFactor,  num)
%     frameEnergy
    if exist('padFactor', 'var')==0 
        num = 1:numel(frameEnergy.E);
    end 

    if ~isfield(frameEnergy, 'E')
        if isfield(frameEnergy, 'vx')
            if exist('num', 'var')==0 
                num = 1:numel(frameEnergy.vx);
            end 
            frameEnergy = calcEnergy(frameEnergy, padFactor, num);
        else
            warning('frame is incorrect!')
            return
        end
    end
    if exist('num', 'var')==0 
        num = 1:numel(frameEnergy.E);
    end 

    if exist('maxk', 'var')==0 
        maxk = 2;
    end 

    kx = frameEnergy.kx;
    ky = frameEnergy.ky;
    
    E = cat(3, frameEnergy.E{:});
    E = mean(E, 3);

    ixy = find(abs(kx)<maxk & abs(ky)<maxk);
    pr_right = 14;
    pr_top = 12;  
    newxsize = length(find(abs(kx(1,:))<maxk));
    newysize = length(find(abs(ky(:,1))<maxk));
    newkx = reshape(kx(ixy),[newysize, newxsize]);
    newky = reshape(ky(ixy),[newysize, newxsize]);
    newE = reshape(E(ixy),[newysize, newxsize]);
    surf(newkx, newky, log10(newE), 'LineStyle', 'None'); view(2);
    set(gca, 'fontsize', 20);
    set(gcf, 'PaperPosition', [0, 0, pr_right, pr_top]);
    colorbar( 'fontsize', 20);
    xlim([-maxk maxk]);
    ylim([-maxk maxk]);   
    xlabel('k_x, cm^{-1}', 'fontsize', 20);
    ylabel('k_y, cm^{-1}', 'fontsize', 20);
        
