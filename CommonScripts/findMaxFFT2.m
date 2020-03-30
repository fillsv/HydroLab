function [c, kx_p, ky_p] = findMaxFFT2(fomega, kx, ky, padFactor, fkx1, fky1)
    
    data = ifft2(ifftshift(fomega));
    data = data(1:size(data,1)/padFactor,1:size(data,2)/padFactor);  
    skx = ifftshift(kx);
    sky = ifftshift(ky);
    dinx = .1;
    diny = .1;

    iny = (interp1(sky(:,1), 1:size(sky,1), fky1)-1)/padFactor;
    inx = (interp1(skx(1,:), 1:size(skx,1), fkx1)-1)/padFactor;
    coefx = .9;
    coefy = .9;
    for uu = 1:1000
        ff1 = goertzel1(data, 1+iny);
        ff2 = goertzel1(ff1, 1+inx);
        if(abs(dinx)>=1e-7)
            ffx1 = goertzel1(data, 1+iny);
            ffx2 = abs(goertzel1(ffx1, 1+inx+dinx));
            gx = (ffx2-abs(ff2))/dinx/abs(ff2);
        else
            gx = 0;
        end

        if(abs(diny)>=1e-7)
            ffy1 = goertzel1(data, 1+iny+diny);
            ffy2 = abs(goertzel1(ffy1, 1+inx));
            gy = (ffy2-abs(ff2))/diny/abs(ff2);
        else
            gy = 0;
        end
        if((abs(dinx)<1e-13)&(abs(diny)<1e-13)) 
            break; 
        end
       
        if(gx*dinx<0) coefx=coefx/2.1; end
        if(gy*diny<0) coefy=coefy/2.1; end

        dinx = coefx*gx;
        diny = coefy*gy;
        inx = inx + dinx;
        iny = iny + diny;

    end 
    dkx = kx(1,2)-kx(1,1);
    dky = ky(2,1)-ky(1,1);
    Nx = size(kx,2)/padFactor;
    Ny = size(ky,1)/padFactor;
    
    fbx = double(dkx/2*padFactor*(Nx));
    fby = double(dky/2*padFactor*(Ny));

    kx_p = double(inx*padFactor*dkx);
    ky_p = double(iny*padFactor*dky);
    
    if (kx_p > fbx) kx_p = kx_p - 2*fbx; end
    if (ky_p > fby) ky_p = ky_p - 2*fby; end
    
    c = abs(ff2);
end