function [xtable ytable utable vtable typevector] = piv_FFTmulti_s_GPU1 (image1,image2,interrogationarea, step, subpixfinder, mask_inpt, roi_inpt, passes,int2,int3,int4,imdeform)


%profile on
%this funtion performs the  PIV analysis.
t = now;
image1 = single(image1);
image2 = single(image2);
interrogationarea = single(interrogationarea);
step = single(step);
subpixfinder = single(subpixfinder);
mask_inpt = single(mask_inpt);
roi_inpt = single(roi_inpt);
passes = single(passes);
int2 = single(int2);
int3 = single(int3);
int4 = single(int4);
% imdeform = single(imdeform);


% warning off %#ok<*WNOFF> %MATLAB:log:logOfZero
if numel(roi_inpt)>0
    xroi=roi_inpt(1);
    yroi=roi_inpt(2);
    widthroi=roi_inpt(3);
    heightroi=roi_inpt(4);
    image1_roi=image1(yroi:yroi+heightroi,xroi:xroi+widthroi);
    image2_roi=image2(yroi:yroi+heightroi,xroi:xroi+widthroi);
else
    xroi=0;
    yroi=0;
    image1_roi=image1;
    image2_roi=image2;
end
clear image1 image2
gen_image1_roi = image1_roi;
gen_image2_roi = image2_roi;

if numel(mask_inpt)>0
    cellmask=mask_inpt;
    mask=single(zeros(size(image1_roi)));
    for ii=1:size(cellmask,1);
        masklayerx=cellmask{ii,1};
        masklayery=cellmask{ii,2};
        mask = mask + single(poly2mask(masklayerx-xroi,masklayery-yroi,size(image1_roi,1),size(image1_roi,2))); %kleineres eingangsbild und maske geshiftet
    end
else
    mask=single(zeros(size(image1_roi)));
end
mask(mask>1)=1;
gen_mask = mask;

miniy=1+(ceil(interrogationarea/2));
minix=1+(ceil(interrogationarea/2));
maxiy=step*(floor(size(image1_roi,1)/step))-(interrogationarea-1)+(ceil(interrogationarea/2)); %statt size deltax von ROI nehmen
maxix=step*(floor(size(image1_roi,2)/step))-(interrogationarea-1)+(ceil(interrogationarea/2));

numelementsy=floor((maxiy-miniy)/step+1);
numelementsx=floor((maxix-minix)/step+1);

LAy=miniy;
LAx=minix;
LUy=size(image1_roi,1)-maxiy;
LUx=size(image1_roi,2)-maxix;
shift4centery=round((LUy-LAy)/2);
shift4centerx=round((LUx-LAx)/2);
if shift4centery<0 %shift4center will be negative if in the unshifted case the left border is bigger than the right border. the vectormatrix is hence not centered on the image. the matrix cannot be shifted more towards the left border because then image2_crop would have a negative index. The only way to center the matrix would be to remove a column of vectors on the right side. but then we weould have less data....
    shift4centery=0;
end
if shift4centerx<0 %shift4center will be negative if in the unshifted case the left border is bigger than the right border. the vectormatrix is hence not centered on the image. the matrix cannot be shifted more towards the left border because then image2_crop would have a negative index. The only way to center the matrix would be to remove a column of vectors on the right side. but then we weould have less data....
    shift4centerx=0;
end
miniy=miniy+shift4centery;
minix=minix+shift4centerx;
maxix=maxix+shift4centerx;
maxiy=maxiy+shift4centery;

image1_roi=padarray(image1_roi, double([ceil(interrogationarea/2) ceil(interrogationarea/2)]), min(min(image1_roi)));
image2_roi=padarray(image2_roi, double([ceil(interrogationarea/2) ceil(interrogationarea/2)]), min(min(image1_roi)));
mask=padarray(mask, double([ceil(interrogationarea/2) ceil(interrogationarea/2)]), 0);

if (rem(interrogationarea,2) == 0) %for the subpixel displacement measurement
    SubPixOffset=1;
else
    SubPixOffset=0.5;
end
% xtable=zeros(numelementsy,numelementsx);
% ytable=xtable;
% utable=xtable;
% vtable=xtable;
typevector=single(ones(numelementsy,numelementsx));

%% MAINLOOP
try %check if used from GUI
    handles=guihandles(getappdata(0,'hgui'));
    GUI_avail=1;
catch %#ok<CTCH>
    GUI_avail=0;
end

% divide images by small pictures
% new index for image1_roi and image2_roi
s0 = (repmat((miniy:step:maxiy)'-1, 1,numelementsx) + repmat(((minix:step:maxix)-1)*size(image1_roi, 1), numelementsy,1))'; 
s0 = permute(s0(:), [2 3 1]);
s1 = repmat((1:interrogationarea)',1,interrogationarea) + repmat(((1:interrogationarea)-1)*size(image1_roi, 1),interrogationarea,1);
ss1 = repmat(s1, 1, 1, size(s0,3)) + repmat(s0, interrogationarea, interrogationarea, 1);
clear s0 s1
image1_cut = image1_roi(ss1);
image2_cut = image2_roi(ss1);
% clear image1_roi image2_roi

%do fft2
result_conv = fftshift(fftshift(real(ifft2(conj(fft2(image1_cut)).*fft2(image2_cut))), 1), 2);
% clear image1_cut image2_cut

minres = permute(repmat(squeeze(min(min(result_conv))), 1, size(result_conv, 1), size(result_conv, 2)), [2 3 1]);
result_conv = result_conv-minres;
clear minres 
maxres = permute(repmat(squeeze(max(max(result_conv))), 1, size(result_conv, 1), size(result_conv, 2)), [2 3 1]);
result_conv = (result_conv./maxres)*255;
clear maxres
%apply mask
ii = find(mask(ss1(round(interrogationarea/2+1), round(interrogationarea/2+1), :)));
jj = find(mask((miniy:step:maxiy)+round(interrogationarea/2), (minix:step:maxix)+round(interrogationarea/2)));
typevector(jj) = 0;
result_conv(:,:, ii) = 0;
clear ii jj
[y, x, z] = ind2sub(size(result_conv), find(result_conv==255));
y = single(y);
x = single(x);
z = single(z);
 % we need only one peak from each couple pictures
[z1, zi] = sort(z);
zi = single(zi);
dz1 = [z1(1); diff(z1)];
i0 = find(dz1~=0);
x1 = x(zi(i0));
y1 = y(zi(i0));
z1 = z(zi(i0));
% clear x y z dz1 i0 zi

xtable = repmat((minix:step:maxix)+interrogationarea/2, length(miniy:step:maxiy), 1);
ytable = repmat(((miniy:step:maxiy)+interrogationarea/2)', 1, length(minix:step:maxix));

if subpixfinder==1
    [vector] = SUBPIXGAUSS (result_conv,interrogationarea, x1, y1, z1, SubPixOffset);
elseif subpixfinder==2
    [vector] = SUBPIX2DGAUSS (result_conv,interrogationarea, x1, y1, z1, SubPixOffset);
end
% clear result_conv x1 y1 z1
vector = permute(reshape(vector, [size(xtable') 2]), [2 1 3]);

utable = vector(:,:,1);
vtable = vector(:,:,2);
% clear vector


%assignin('base','corr_results',corr_results);


%multipass
%feststellen wie viele passes
%wenn intarea=0 dann keinen pass.
for multipass=1:passes-1

    if GUI_avail==1
        set(handles.progress, 'string' , ['Frame progress: ' int2str(j/maxiy*100/passes+((multipass-1)*(100/passes))) '%' sprintf('\n') 'Validating velocity field']);drawnow;
     else
        fprintf('.');
    end
    %multipass validation, smoothing
    %stdev test
    utable_orig=utable;
    vtable_orig=vtable;
    stdthresh=4;
    utable(find(isinf(utable))) = NaN;
    vtable(find(isinf(vtable))) = NaN;     
    meanu=nanmean(nanmean(utable));
    meanv=nanmean(nanmean(vtable));
    std2u=nanstd(reshape(utable,size(utable,1)*size(utable,2),1));
    std2v=nanstd(reshape(vtable,size(vtable,1)*size(vtable,2),1));
    minvalu=meanu-stdthresh*std2u;
    maxvalu=meanu+stdthresh*std2u;
    minvalv=meanv-stdthresh*std2v;
    maxvalv=meanv+stdthresh*std2v;
    utable(utable<minvalu)=NaN;
    utable(utable>maxvalu)=NaN;
    vtable(vtable<minvalv)=NaN;
    vtable(vtable>maxvalv)=NaN;
    
    %median test
    %info1=[];
    epsilon=0.02;
    thresh=2;
    [J,I]=size(utable);
    %medianres=zeros(J,I);
    normfluct=zeros(J,I,2);
    b=1;
    %eps=0.1;
    for c=1:2
        if c==1;
            velcomp=utable;
        else
            velcomp=vtable;
        end
        clear neigh
        for ii = -b:b;
            for jj = -b:b;
                neigh(:, :, ii+2*b, jj+2*b)=velcomp((1+b:end-b)+ii, (1+b:end-b)+jj);
            end
        end

        neighcol = reshape(neigh, size(neigh,1), size(neigh,2), (2*b+1)^2);
        clear neigh
        neighcol2= neighcol(:,:, [(1:(2*b+1)*b+b) ((2*b+1)*b+b+2:(2*b+1)^2)]);
        clear neighcol
        neighcol2 = permute(neighcol2, [3, 1, 2]);
        med=median(neighcol2);
        velcomp = velcomp((1+b:end-b), (1+b:end-b));
        fluct=velcomp-permute(med, [2 3 1]);
        res=neighcol2-repmat(med, (2*b+1)^2-1, 1,1);
        clear neighcol2 med velcomp
        medianres=permute(median(abs(res)), [2 3 1]);
        normfluct((1+b:end-b), (1+b:end-b), c)=abs(fluct./(medianres+epsilon));
    end
    
        
    info1=(sqrt(normfluct(:,:,1).^2+normfluct(:,:,2).^2)>thresh);
    clear res fluct normfluct
    utable(info1==1) = NaN;
    vtable(info1==1) = NaN;
    clear info1
    %find typevector...
    %maskedpoints=numel(find((typevector)==0));
    %amountnans=numel(find(isnan(utable)==1))-maskedpoints;
    %discarded=amountnans/(size(utable,1)*size(utable,2))*100;
    %disp(['Discarded: ' num2str(amountnans) ' vectors = ' num2str(discarded) ' %'])
    
    if GUI_avail==1
        if verLessThan('matlab','8.4')
            delete (findobj(getappdata(0,'hgui'),'type', 'hggroup'))
        else
              delete (findobj(getappdata(0,'hgui'),'type', 'quiver'))
        end
        hold on;
        vecscale=str2double(get(handles.vectorscale,'string'));
        %Problem: wenn colorbar an, z�hlt das auch als aexes...
        colorbar('off')
        quiver ((findobj(getappdata(0,'hgui'),'type', 'axes')),xtable(isnan(utable)==0)+xroi-interrogationarea/2,ytable(isnan(utable)==0)+yroi-interrogationarea/2,utable_orig(isnan(utable)==0)*vecscale,vtable_orig(isnan(utable)==0)*vecscale,'Color', [0.15 0.7 0.15],'autoscale','off')
        quiver ((findobj(getappdata(0,'hgui'),'type', 'axes')),xtable(isnan(utable)==1)+xroi-interrogationarea/2,ytable(isnan(utable)==1)+yroi-interrogationarea/2,utable_orig(isnan(utable)==1)*vecscale,vtable_orig(isnan(utable)==1)*vecscale,'Color',[0.7 0.15 0.15], 'autoscale','off')
        drawnow
        hold off
    end
    
    %replace nans
    utable = single(inpaint_nans(double(utable),4));
    vtable = single(inpaint_nans(double(vtable),4));
    %smooth predictor
    try
        if multipass<passes-1
            utable = single(smoothn(utable,0.6)); %stronger smoothing for first passes
            vtable = single(smoothn(vtable,0.6));
        else
            utable = single(smoothn(utable)); %weaker smoothing for last pass
            vtable = single(smoothn(vtable));
        end
    catch
        
        %old matlab versions: gaussian kernel
        h=fspecial('gaussian',5,1);
        utable=single(imfilter(utable,h,'replicate'));
        vtable=single(imfilter(vtable,h,'replicate'));
    end

    if multipass==1
        interrogationarea=round(int2/2)*2;
    end
    if multipass==2
        interrogationarea=round(int3/2)*2;
    end
    if multipass==3
        interrogationarea=round(int4/2)*2;
    end
    step=interrogationarea/2;
    
    %bildkoordinaten neu errechnen:
    %roi=[];

    image1_roi = gen_image1_roi;
    image2_roi = gen_image2_roi;
    mask = gen_mask;
    
    
    miniy=1+(ceil(interrogationarea/2));
    minix=1+(ceil(interrogationarea/2));
    maxiy=step*(floor(size(image1_roi,1)/step))-(interrogationarea-1)+(ceil(interrogationarea/2)); %statt size deltax von ROI nehmen
    maxix=step*(floor(size(image1_roi,2)/step))-(interrogationarea-1)+(ceil(interrogationarea/2));
    
    numelementsy=floor((maxiy-miniy)/step+1);
    numelementsx=floor((maxix-minix)/step+1);
    
    LAy=miniy;
    LAx=minix;
    LUy=size(image1_roi,1)-maxiy;
    LUx=size(image1_roi,2)-maxix;
    shift4centery=round((LUy-LAy)/2);
    shift4centerx=round((LUx-LAx)/2);
    if shift4centery<0 %shift4center will be negative if in the unshifted case the left border is bigger than the right border. the vectormatrix is hence not centered on the image. the matrix cannot be shifted more towards the left border because then image2_crop would have a negative index. The only way to center the matrix would be to remove a column of vectors on the right side. but then we weould have less data....
        shift4centery=0;
    end
    if shift4centerx<0 %shift4center will be negative if in the unshifted case the left border is bigger than the right border. the vectormatrix is hence not centered on the image. the matrix cannot be shifted more towards the left border because then image2_crop would have a negative index. The only way to center the matrix would be to remove a column of vectors on the right side. but then we weould have less data....
        shift4centerx=0;
    end
    miniy=miniy+shift4centery;
    minix=minix+shift4centerx;
    maxix=maxix+shift4centerx;
    maxiy=maxiy+shift4centery;
    
    image1_roi=padarray(image1_roi, double([ceil(interrogationarea/2) ceil(interrogationarea/2)]), min(min(image1_roi)));
    image2_roi=padarray(image2_roi, double([ceil(interrogationarea/2) ceil(interrogationarea/2)]) , min(min(image1_roi)));
    mask=padarray(mask, double([ceil(interrogationarea/2) ceil(interrogationarea/2)]), 0);
    if (rem(interrogationarea,2) == 0) %for the subpixel displacement measurement
        SubPixOffset=1;
    else
        SubPixOffset=0.5;
    end
    
    xtable_old = xtable;
    ytable_old = ytable;
    typevector = single(ones(numelementsy,numelementsx));
    xtable = repmat((minix:step:maxix), numelementsy, 1) + interrogationarea/2;
    ytable = repmat((miniy:step:maxiy)', 1, numelementsx) + interrogationarea/2;

    %xtable alt und neu geben koordinaten wo die vektoren herkommen.
    %d.h. u und v auf die gew�nschte gr��e bringen+interpolieren
    if GUI_avail==1
        set(handles.progress, 'string' , ['Frame progress: ' int2str(j/maxiy*100/passes+((multipass-1)*(100/passes))) '%' sprintf('\n') 'Interpolating velocity field']);drawnow;
        %set(handles.progress, 'string' , 'Interpolating velocity field');drawnow;
    else
        fprintf('.');
    end

    utable = interp2(xtable_old,ytable_old,utable,xtable,ytable,'*spline');
    vtable = interp2(xtable_old,ytable_old,vtable,xtable,ytable,'*spline');
    clear xtable_old ytable_old
   
    %add 1 line around image for border regions... linear extrap
    
    firstlinex = xtable(1,:);
    firstlinex_intp=interp1(1:1:size(firstlinex,2),firstlinex,0:1:size(firstlinex,2)+1,'linear','extrap');
    xtable_1=repmat(firstlinex_intp,size(xtable,1)+2,1);
    clear firstlinex firstlinex_intp
    
    firstliney=ytable(:,1);
    firstliney_intp=interp1(1:1:size(firstliney,1),firstliney,0:1:size(firstliney,1)+1,'linear','extrap')';
    ytable_1=repmat(firstliney_intp,1,size(ytable,2)+2);

    clear firstliney firstliney_intp
    
    X = xtable_1; %original locations of vectors in whole image
    Y = ytable_1;
    
    U = padarray(utable, [1,1], 'replicate'); %interesting portion of u
    V = padarray(vtable, [1,1], 'replicate'); % "" of v
    
    X1=single(X(1,1):1:X(1,end)-1); 
    Y1=single((Y(1,1):1:Y(end,1)-1))';
    X1=repmat(X1,size(Y1, 1),1);
    Y1=repmat(Y1,1,size(X1, 2));

    U1 = interp2(X,Y,U,X1,Y1,'*linear');
    V1 = interp2(X,Y,V,X1,Y1,'*linear');
    clear X Y U V
    
    image2_crop_i1 = interp2(1:size(image2_roi,2), (1:size(image2_roi,1))', image2_roi, X1+U1, Y1+V1, imdeform); %linear is 3x faster and looks ok...
    clear U1 V1
    xb = find(X1(1,:) == xtable_1(1,1));
    yb = find(Y1(:,1) == ytable_1(1,1));
    clear X1 Y1 ytable_1 xtable_1
    % divide images by small pictures
    % new index for image1_roi
    s0 = (repmat((miniy:step:maxiy)'-1, 1,numelementsx) + repmat(((minix:step:maxix)-1)*size(image1_roi, 1), numelementsy,1))'; 
    s0 = permute(s0(:), [2 3 1]);
    s1 = repmat((1:interrogationarea)',1,interrogationarea) + repmat(((1:interrogationarea)-1)*size(image1_roi, 1),interrogationarea,1);
    ss1 = repmat(s1, 1, 1, size(s0,3)) + repmat(s0, interrogationarea, interrogationarea, 1);
    clear s1
    % new index for image2_crop_i1
    s0 = (repmat(yb-step+step*(1:numelementsy)'-1, 1,numelementsx) + repmat((xb-step+step*(1:numelementsx)-1)*size(image2_crop_i1, 1), numelementsy,1))'; 
    s0 = permute(s0(:), [2 3 1]) - s0(1);
    s2 = repmat((1:2*step)',1,2*step) + repmat(((1:2*step)-1)*size(image2_crop_i1, 1),2*step,1);
    ss2 = repmat(s2, 1, 1, size(s0,3)) + repmat(s0, interrogationarea, interrogationarea, 1);
    clear s0 s2
    image1_cut = image1_roi(ss1);
    image2_cut = image2_crop_i1(ss2);
    clear ss2

    %do fft2
    result_conv = fftshift(fftshift(real(ifft2(conj(fft2(image1_cut)).*fft2(image2_cut))), 1), 2);
    minres = permute(repmat(squeeze(min(min(result_conv))), 1, size(result_conv, 1), size(result_conv, 2)), [2 3 1]);
    result_conv = result_conv-minres;
    clear minres 
    maxres = permute(repmat(squeeze(max(max(result_conv))), 1, size(result_conv, 1), size(result_conv, 2)), [2 3 1]);
    result_conv = (result_conv./maxres)*255;
    clear maxres
    
   

    %apply mask
    ii = single(find(mask(ss1(round(interrogationarea/2+1), round(interrogationarea/2+1), :))));
    jj = single(find(mask((miniy:step:maxiy)+round(interrogationarea/2), (minix:step:maxix)+round(interrogationarea/2))));
    typevector(jj) = 0;
    result_conv(:,:, ii) = 0;
    clear ss1

    [y, x, z] = ind2sub(size(result_conv), find(result_conv==255));
    x = single(x);
    y = single(y);
    z = single(z);
    
    [z1, zi] = sort(z);
    zi = single(zi);
    % we need only one peak from each couple pictures
    dz1 = [z1(1); diff(z1)];
    i0 = single(find(dz1~=0));
    x1 = x(zi(i0));
    y1 = y(zi(i0));
    z1 = z(zi(i0));
    clear x y z
    %new xtable and ytable
    xtable = repmat((minix:step:maxix)+interrogationarea/2, length(miniy:step:maxiy), 1);
    ytable = repmat(((miniy:step:maxiy)+interrogationarea/2)', 1, length(minix:step:maxix));

    if subpixfinder==1
        [vector] = SUBPIXGAUSS (result_conv,interrogationarea, x1, y1, z1,SubPixOffset);
    elseif subpixfinder==2
        [vector] = SUBPIX2DGAUSS (result_conv,interrogationarea, x1, y1, z1,SubPixOffset);
    end
    vector = permute(reshape(vector, [size(xtable') 2]), [2 1 3]);
    clear x1 y1 z1 result_conv
    utable = utable+vector(:,:,1);
    vtable = vtable+vector(:,:,2);
    clear vector

end

%assignin('base','pass_result',pass_result);
%__________________________________________________________________________


xtable=xtable-ceil(interrogationarea/2);
ytable=ytable-ceil(interrogationarea/2);

xtable=xtable+xroi;
ytable=ytable+yroi;
fprintf('\n');
disp((now - t)*24*3600);

end
%profile viewer
%p = profile('info');
%profsave(p,'profile_results')

function [vector] = SUBPIXGAUSS(result_conv, interrogationarea, x, y, z, SubPixOffset)
    xi = single(find(~((x <= (size(result_conv,2)-1)) & (y <= (size(result_conv,1)-1)) & (x >= 2) & (y >= 2))));
    x(xi) = [];
    y(xi) = [];
    z(xi) = [];
    xmax = size(result_conv, 2);
    vector=NaN*(single(ones(size(result_conv,3),2)));
    if(numel(x)~=0)
        ip = sub2ind(size(result_conv), y, x, z);
%         the following 8 lines are copyright (c) 1998, Uri Shavit, Roi Gurka, Alex Liberzon, Technion � Israel Institute of Technology
%         http://urapiv.wordpress.com
        f0 = log(result_conv(ip));
        f1 = log(result_conv(ip-1));
        f2 = log(result_conv(ip+1));
        peaky = y + (f1-f2)./(2*f1-4*f0+2*f2);
        f0 = log(result_conv(ip));
        f1 = log(result_conv(ip-xmax));
        f2 = log(result_conv(ip+xmax));
        peakx = x + (f1-f2)./(2*f1-4*f0+2*f2);

        SubpixelX=peakx-(interrogationarea/2)-SubPixOffset;
        SubpixelY=peaky-(interrogationarea/2)-SubPixOffset;
        vector(z, :) = [SubpixelX, SubpixelY];  
    end
end
    
function [vector] = SUBPIX2DGAUSS(result_conv, interrogationarea, x, y, z, SubPixOffset)
    xi = single(find(~((x <= (size(result_conv,2)-1)) & (y <= (size(result_conv,1)-1)) & (x >= 2) & (y >= 2))));
    x(xi) = [];
    y(xi) = [];
    z(xi) = [];
    xmax = size(result_conv, 2);
    vector = single(NaN*(ones(size(result_conv,3),2)));
    if(numel(x)~=0)
        c10 = single(zeros(3,3, length(z)));
        c01 = c10;
        c11 = c10;
        c20 = c10;
        c02 = c10;
        ip = sub2ind(size(result_conv), y, x, z);

        for ii = -1:1
            for j = -1:1
                %following 15 lines based on
                %H. Nobach � M. Honkanen (2005)
                %Two-dimensional Gaussian regression for sub-pixel displacement
                %estimation in particle image velocimetry or particle position
                %estimation in particle tracking velocimetry
                %Experiments in Fluids (2005) 38: 511�515
                c10(j+2,ii+2, :) = ii*log(result_conv(ip+xmax*ii+j));
                c01(j+2,ii+2, :) = j*log(result_conv(ip+xmax*ii+j));
                c11(j+2,ii+2, :) = ii*j*log(result_conv(ip+xmax*ii+j));
                c20(j+2,ii+2, :) = (3*ii^2-2)*log(result_conv(ip+xmax*ii+j));
                c02(j+2,ii+2, :) = (3*j^2-2)*log(result_conv(ip+xmax*ii+j));
                %c00(j+2,i+2)=(5-3*i^2-3*j^2)*log(result_conv_norm(maxY+j, maxX+i));
            end
        end
        c10 = (1/6)*sum(sum(c10));
        c01 = (1/6)*sum(sum(c01));
        c11 = (1/4)*sum(sum(c11));
        c20 = (1/6)*sum(sum(c20));
        c02 = (1/6)*sum(sum(c02));
        %c00=(1/9)*sum(sum(c00));

        deltax = squeeze((c11.*c01-2*c10.*c02)./(4*c20.*c02-c11.^2));
        deltay = squeeze((c11.*c10-2*c01.*c20)./(4*c20.*c02-c11.^2));
        peakx = x+deltax;
        peaky = y+deltay;

        SubpixelX = peakx-(interrogationarea/2)-SubPixOffset;
        SubpixelY = peaky-(interrogationarea/2)-SubPixOffset;

        vector(z, :) = [SubpixelX, SubpixelY];
    end
end
    
function out = PIVlab_preproc (in,roirect,clahe, clahesize,highp,highpsize,intenscap,wienerwurst,wienerwurstsize)
if size(in,3)>1
    in(:,:,2:end)=[];
end
%this function preprocesses the images
if numel(roirect)>0
    x=roirect(1);
    y=roirect(2);
    width=roirect(3);
    height=roirect(4);
else
    x=1;
    y=1;
    width=size(in,2)-1;
    height=size(in,1)-1;
end
%roi (x,y,width,height)
in_roi=in(y:y+height,x:x+width);

if intenscap == 1
    %Intensity Capping: a simple method to improve cross-correlation PIV results
    %Uri Shavit � Ryan J. Lowe � Jonah V. Steinbuck
    n = 2; 
    up_lim_im_1 = median(double(in_roi(:))) + n*std2(in_roi); % upper limit for image 1
    brightspots_im_1 = find(in_roi > up_lim_im_1); % bright spots in image 1
    capped_im_1 = in_roi; capped_im_1(brightspots_im_1) = up_lim_im_1; % capped image 1
    in_roi=capped_im_1;
end
if clahe == 1
    numberoftiles1=round(size(in_roi,1)/clahesize);
    numberoftiles2=round(size(in_roi,2)/clahesize);
    if numberoftiles1 < 2
    numberoftiles1=2;
    end
    if numberoftiles2 < 2
    numberoftiles2=2;
    end
    in_roi=adapthisteq(in_roi, 'NumTiles',[numberoftiles1 numberoftiles2], 'ClipLimit', 0.01, 'NBins', 256, 'Range', 'full', 'Distribution', 'uniform');
end

if highp == 1
    h = fspecial('gaussian',highpsize,highpsize);
    in_roi=double(in_roi-(imfilter(in_roi,h,'replicate')));
    in_roi=in_roi/max(max(in_roi))*255;
end

if wienerwurst == 1
    in_roi=wiener2(in_roi,[wienerwurstsize wienerwurstsize]);
end

out=in;
out(y:y+height,x:x+width)=in_roi;
out=uint8(out);
end

function B = inpaint_nans(A,method)

% Author: John D'Errico
% e-mail address: woodchips@rochester.rr.com
% Release: 2
% Release date: 4/15/06

% INPAINT_NANS: in-paints over nans in an array
% usage: B=INPAINT_NANS(A)          % default method
% usage: B=INPAINT_NANS(A,method)   % specify method used
%
% Solves approximation to one of several pdes to
% interpolate and extrapolate holes in an array
%
% arguments (input):
%   A - nxm array with some NaNs to be filled in
%
%   method - (OPTIONAL) scalar numeric flag - specifies
%       which approach (or physical metaphor to use
%       for the interpolation.) All methods are capable
%       of extrapolation, some are better than others.
%       There are also speed differences, as well as
%       accuracy differences for smooth surfaces.
%
%       methods {0,1,2} use a simple plate metaphor.
%       method  3 uses a better plate equation,
%                 but may be much slower and uses
%                 more memory.
%       method  4 uses a spring metaphor.
%       method  5 is an 8 neighbor average, with no
%                 rationale behind it compared to the
%                 other methods. I do not recommend
%                 its use.
%
%       method == 0 --> (DEFAULT) see method 1, but
%         this method does not build as large of a
%         linear system in the case of only a few
%         NaNs in a large array.
%         Extrapolation behavior is linear.
%         
%       method == 1 --> simple approach, applies del^2
%         over the entire array, then drops those parts
%         of the array which do not have any contact with
%         NaNs. Uses a least squares approach, but it
%         does not modify known values.
%         In the case of small arrays, this method is
%         quite fast as it does very little extra work.
%         Extrapolation behavior is linear.
%         
%       method == 2 --> uses del^2, but solving a direct
%         linear system of equations for nan elements.
%         This method will be the fastest possible for
%         large systems since it uses the sparsest
%         possible system of equations. Not a least
%         squares approach, so it may be least robust
%         to noise on the boundaries of any holes.
%         This method will also be least able to
%         interpolate accurately for smooth surfaces.
%         Extrapolation behavior is linear.
%         
%       method == 3 --+ See method 0, but uses del^4 for
%         the interpolating operator. This may result
%         in more accurate interpolations, at some cost
%         in speed.
%         
%       method == 4 --+ Uses a spring metaphor. Assumes
%         springs (with a nominal length of zero)
%         connect each node with every neighbor
%         (horizontally, vertically and diagonally)
%         Since each node tries to be like its neighbors,
%         extrapolation is as a constant function where
%         this is consistent with the neighboring nodes.
%
%       method == 5 --+ See method 2, but use an average
%         of the 8 nearest neighbors to any element.
%         This method is NOT recommended for use.
%
%
% arguments (output):
%   B - nxm array with NaNs replaced
%
%
% Example:
%  [x,y] = meshgrid(0:.01:1);
%  z0 = exp(x+y);
%  znan = z0;
%  znan(20:50,40:70) = NaN;
%  znan(30:90,5:10) = NaN;
%  znan(70:75,40:90) = NaN;
%
%  z = inpaint_nans(znan);
%
%
% See also: griddata, interp1
%
% Author: John D'Errico
% e-mail address: woodchips@rochester.rr.com
% Release: 2
% Release date: 4/15/06


% I always need to know which elements are NaN,
% and what size the array is for any method
[n,m]=size(A);
A=A(:);
nm=n*m;
k=isnan(A(:));

% list the nodes which are known, and which will
% be interpolated
nan_list=find(k);
known_list=find(~k);

% how many nans overall
nan_count=length(nan_list);

% convert NaN indices to (r,c) form
% nan_list==find(k) are the unrolled (linear) indices
% (row,column) form
[nr,nc]=ind2sub([n,m],nan_list);

% both forms of index in one array:
% column 1 == unrolled index
% column 2 == row index
% column 3 == column index
nan_list=[nan_list,nr,nc];

% supply default method
if (nargin<2) || isempty(method)
  method = 0;
elseif ~ismember(method,0:5)
  error 'If supplied, method must be one of: {0,1,2,3,4,5}.'
end

% for different methods
switch method
 case 0
  % The same as method == 1, except only work on those
  % elements which are NaN, or at least touch a NaN.
  
  % horizontal and vertical neighbors only
  talks_to = [-1 0;0 -1;1 0;0 1];
  neighbors_list=identify_neighbors(n,m,nan_list,talks_to);
  
  % list of all nodes we have identified
  all_list=[nan_list;neighbors_list];
  
  % generate sparse array with second partials on row
  % variable for each element in either list, but only
  % for those nodes which have a row index > 1 or < n
  L = find((all_list(:,2) > 1) & (all_list(:,2) < n)); 
  nl=length(L);
  if nl>0
    fda=sparse(repmat(all_list(L,1),1,3), ...
      repmat(all_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
      repmat([1 -2 1],nl,1),nm,nm);
  else
    fda=spalloc(n*m,n*m,size(all_list,1)*5);
  end
  
  % 2nd partials on column index
  L = find((all_list(:,3) > 1) & (all_list(:,3) < m)); 
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(all_list(L,1),1,3), ...
      repmat(all_list(L,1),1,3)+repmat([-n 0 n],nl,1), ...
      repmat([1 -2 1],nl,1),nm,nm);
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  k=find(any(fda(:,nan_list(:,1)),2));
  
  % and solve...
  B=A;
  B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
  
 case 1
  % least squares approach with del^2. Build system
  % for every array element as an unknown, and then
  % eliminate those which are knowns.

  % Build sparse matrix approximating del^2 for
  % every element in A.
  % Compute finite difference for second partials
  % on row variable first
  [ii,j]=ndgrid(2:(n-1),1:m);
  ind=ii(:)+(j(:)-1)*n;
  np=(n-2)*m;
  fda=sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
      repmat([1 -2 1],np,1),n*m,n*m);
  
  % now second partials on column variable
  [ii,j]=ndgrid(1:n,2:(m-1));
  ind=ii(:)+(j(:)-1)*n;
  np=n*(m-2);
  fda=fda+sparse(repmat(ind,1,3),[ind-n,ind,ind+n], ...
      repmat([1 -2 1],np,1),nm,nm);
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  k=find(any(fda(:,nan_list),2));
  
  % and solve...
  B=A;
  B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
  
 case 2
  % Direct solve for del^2 BVP across holes

  % generate sparse array with second partials on row
  % variable for each nan element, only for those nodes
  % which have a row index > 1 or < n
  L = find((nan_list(:,2) > 1) & (nan_list(:,2) < n)); 
  nl=length(L);
  if nl>0
    fda=sparse(repmat(nan_list(L,1),1,3), ...
      repmat(nan_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
      repmat([1 -2 1],nl,1),n*m,n*m);
  else
    fda=spalloc(n*m,n*m,size(nan_list,1)*5);
  end
  
  % 2nd partials on column index
  L = find((nan_list(:,3) > 1) & (nan_list(:,3) < m)); 
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,3), ...
      repmat(nan_list(L,1),1,3)+repmat([-n 0 n],nl,1), ...
      repmat([1 -2 1],nl,1),n*m,n*m);
  end
  
  % fix boundary conditions at extreme corners
  % of the array in case there were nans there
  if ismember(1,nan_list(:,1))
    fda(1,[1 2 n+1])=[-2 1 1];
  end
  if ismember(n,nan_list(:,1))
    fda(n,[n, n-1,n+n])=[-2 1 1];
  end
  if ismember(nm-n+1,nan_list(:,1))
    fda(nm-n+1,[nm-n+1,nm-n+2,nm-n])=[-2 1 1];
  end
  if ismember(nm,nan_list(:,1))
    fda(nm,[nm,nm-1,nm-n])=[-2 1 1];
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  
  % and solve...
  B=A;
  k=nan_list(:,1);
  B(k)=fda(k,k)\rhs(k);
  
 case 3
  % The same as method == 0, except uses del^4 as the
  % interpolating operator.
  
  % del^4 template of neighbors
  talks_to = [-2 0;-1 -1;-1 0;-1 1;0 -2;0 -1; ...
      0 1;0 2;1 -1;1 0;1 1;2 0];
  neighbors_list=identify_neighbors(n,m,nan_list,talks_to);
  
  % list of all nodes we have identified
  all_list=[nan_list;neighbors_list];
  
  % generate sparse array with del^4, but only
  % for those nodes which have a row & column index
  % >= 3 or <= n-2
  L = find( (all_list(:,2) >= 3) & ...
            (all_list(:,2) <= (n-2)) & ...
            (all_list(:,3) >= 3) & ...
            (all_list(:,3) <= (m-2)));
  nl=length(L);
  if nl>0
    % do the entire template at once
    fda=sparse(repmat(all_list(L,1),1,13), ...
        repmat(all_list(L,1),1,13) + ...
        repmat([-2*n,-n-1,-n,-n+1,-2,-1,0,1,2,n-1,n,n+1,2*n],nl,1), ...
        repmat([1 2 -8 2 1 -8 20 -8 1 2 -8 2 1],nl,1),nm,nm);
  else
    fda=spalloc(n*m,n*m,size(all_list,1)*5);
  end
  
  % on the boundaries, reduce the order around the edges
  L = find((((all_list(:,2) == 2) | ...
             (all_list(:,2) == (n-1))) & ...
            (all_list(:,3) >= 2) & ...
            (all_list(:,3) <= (m-1))) | ...
           (((all_list(:,3) == 2) | ...
             (all_list(:,3) == (m-1))) & ...
            (all_list(:,2) >= 2) & ...
            (all_list(:,2) <= (n-1))));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(all_list(L,1),1,5), ...
      repmat(all_list(L,1),1,5) + ...
        repmat([-n,-1,0,+1,n],nl,1), ...
      repmat([1 1 -4 1 1],nl,1),nm,nm);
  end
  
  L = find( ((all_list(:,2) == 1) | ...
             (all_list(:,2) == n)) & ...
            (all_list(:,3) >= 2) & ...
            (all_list(:,3) <= (m-1)));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(all_list(L,1),1,3), ...
      repmat(all_list(L,1),1,3) + ...
        repmat([-n,0,n],nl,1), ...
      repmat([1 -2 1],nl,1),nm,nm);
  end
  
  L = find( ((all_list(:,3) == 1) | ...
             (all_list(:,3) == m)) & ...
            (all_list(:,2) >= 2) & ...
            (all_list(:,2) <= (n-1)));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(all_list(L,1),1,3), ...
      repmat(all_list(L,1),1,3) + ...
        repmat([-1,0,1],nl,1), ...
      repmat([1 -2 1],nl,1),nm,nm);
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  k=find(any(fda(:,nan_list(:,1)),2));
  
  % and solve...
  B=A;
  B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
  
 case 4
  % Spring analogy
  % interpolating operator.
  
  % list of all springs between a node and a horizontal
  % or vertical neighbor
  hv_list=[-1 -1 0;1 1 0;-n 0 -1;n 0 1];
  hv_springs=[];
  for ii=1:4
    hvs=nan_list+repmat(hv_list(ii,:),nan_count,1);
    k=(hvs(:,2)>=1) & (hvs(:,2)<=n) & (hvs(:,3)>=1) & (hvs(:,3)<=m);
    hv_springs=[hv_springs;[nan_list(k,1),hvs(k,1)]];
  end

  % delete replicate springs
  hv_springs=unique(sort(hv_springs,2),'rows');
  
  % build sparse matrix of connections, springs
  % connecting diagonal neighbors are weaker than
  % the horizontal and vertical springs
  nhv=size(hv_springs,1);
  springs=sparse(repmat((1:nhv)',1,2),hv_springs, ...
     repmat([1 -1],nhv,1),nhv,nm);
  
  % eliminate knowns
  rhs=-springs(:,known_list)*A(known_list);
  
  % and solve...
  B=A;
  B(nan_list(:,1))=springs(:,nan_list(:,1))\rhs;
  
 case 5
  % Average of 8 nearest neighbors
  
  % generate sparse array to average 8 nearest neighbors
  % for each nan element, be careful around edges
  fda=spalloc(n*m,n*m,size(nan_list,1)*9);
  
  % -1,-1
  L = find((nan_list(:,2) > 1) & (nan_list(:,3) > 1)); 
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-n-1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end
  
  % 0,-1
  L = find(nan_list(:,3) > 1);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-n, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % +1,-1
  L = find((nan_list(:,2) < n) & (nan_list(:,3) > 1));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-n+1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % -1,0
  L = find(nan_list(:,2) > 1);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % +1,0
  L = find(nan_list(:,2) < n);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % -1,+1
  L = find((nan_list(:,2) > 1) & (nan_list(:,3) < m)); 
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([n-1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end
  
  % 0,+1
  L = find(nan_list(:,3) < m);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([n, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % +1,+1
  L = find((nan_list(:,2) < n) & (nan_list(:,3) < m));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([n+1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  
  % and solve...
  B=A;
  k=nan_list(:,1);
  B(k)=fda(k,k)\rhs(k);
  
end

% all done, make sure that B is the same shape as
% A was when we came in.
B=reshape(B,n,m);
end


function neighbors_list=identify_neighbors(n,m,nan_list,talks_to)
% identify_neighbors: identifies all the neighbors of
%   those nodes in nan_list, not including the nans
%   themselves
%
% arguments (input):
%  n,m - scalar - [n,m]=size(A), where A is the
%      array to be interpolated
%  nan_list - array - list of every nan element in A
%      nan_list(i,1) == linear index of i'th nan element
%      nan_list(i,2) == row index of i'th nan element
%      nan_list(i,3) == column index of i'th nan element
%  talks_to - px2 array - defines which nodes communicate
%      with each other, i.e., which nodes are neighbors.
%
%      talks_to(i,1) - defines the offset in the row
%                      dimension of a neighbor
%      talks_to(i,2) - defines the offset in the column
%                      dimension of a neighbor
%      
%      For example, talks_to = [-1 0;0 -1;1 0;0 1]
%      means that each node talks only to its immediate
%      neighbors horizontally and vertically.
% 
% arguments(output):
%  neighbors_list - array - list of all neighbors of
%      all the nodes in nan_list

if ~isempty(nan_list)
  % use the definition of a neighbor in talks_to
  nan_count=size(nan_list,1);
  talk_count=size(talks_to,1);
  
  nn=zeros(nan_count*talk_count,2);
  j=[1,nan_count];
  for ii=1:talk_count
    nn(j(1):j(2),:)=nan_list(:,2:3) + ...
        repmat(talks_to(ii,:),nan_count,1);
    j=j+nan_count;
  end
  
  % drop those nodes which fall outside the bounds of the
  % original array
  L = (nn(:,1)<1)|(nn(:,1)>n)|(nn(:,2)<1)|(nn(:,2)>m); 
  nn(L,:)=[];
  
  % form the same format 3 column array as nan_list
  neighbors_list=[sub2ind([n,m],nn(:,1),nn(:,2)),nn];
  
  % delete replicates in the neighbors list
  neighbors_list=unique(neighbors_list,'rows');
  
  % and delete those which are also in the list of NaNs.
  neighbors_list=setdiff(neighbors_list,nan_list,'rows');
  
else
  neighbors_list=[];
end

end










