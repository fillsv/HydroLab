function [reza phi] = angleLine(rez, kx, ky, dph, k, dk)
    kk = 0;
    kz = kx+1i*ky;
%     tic
    indk = find((abs(abs(kz)-k)<=dk));
    
%     toc
    [akz indk1] = sort(angle(kz(indk)));
    rez = rez(indk(indk1));
    phi = zeros(1, numel(-pi:dph:pi-dph));
    reza = zeros(size(phi));
    for kk = 1:numel(-pi:dph:pi-dph);
        ph = -pi+(kk-1)*dph;
        phi(kk) = ph+dph/2;
        inda = find((akz>=ph)&(akz<ph+dph));
%         inda1 = findnear(akz, ph);
%         inda2 = findnear(akz, ph+dph);
        if isempty(inda)
            reza(kk) = NaN;
        else
            reza(kk) = mean(abs(rez(inda)));
%             reza(kk) = mean(abs(rez(inda1:inda2)));
        end
    end
        

end