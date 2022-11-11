function [frameCoh, frameIncoh] = calcCoh_Incoh(frameOmega, num, keep, deep, wavelette )
%Рассчет когерентной и некогерентной части DWT
    if exist('num', 'var')==0
        num = 1:numel(frameOmega.px);
    end
    if exist('deep', 'var')==0
        deep = 3; % По гафикам показывает лучшее показания
    end
    if exist('wavelette', 'var')==0
        wavelette = 'coif2';
    end
    px = frameOmega.px{num(1)};
    py = frameOmega.py{num(1)};
    for ii = 1:numel(num)
        omega = cat(3, frameOmega.omega{num(ii)});
        omega(find(isnan(omega))) = 0;
        
%         omega1 = zeros(256,256);
%         omega1(1:209, 1:209) = omega1(1:209, 1:209)+omega;
%         
%         omega = omega1;
%         
        
         [C,S] = wavedec2(omega, deep, wavelette);
%        [C,S] = wavedec2(omega, 7, wavelette);
        Coeff_sort =sort(abs(C(:)));
        
%         if exist('keep', 'var')==0
        if(keep==0)
%             keep = 0.05;
%             keep = calcN0_from_Entropy(omega(:))/length(omega(:));
             keep = calcN0_from_Entropy(C(:))/length(C(:));
            disp("Keep " + keep*100 + "%");
%             break;
        end
        
        
        thresh = Coeff_sort(floor((1-keep)*length(Coeff_sort)));
        index = abs(C)>thresh;
        C_filter = C.*index;
        BB = waverec2(C_filter, S, wavelette);
        
        
%         frameCoh.omega{ii} = BB(1:209, 1:209);
%         frameIncoh.omega{ii} = omega - BB(1:209, 1:209);
        frameCoh.omega{ii} = BB;
        frameIncoh.omega{ii} = omega - BB;
        
        frameCoh.px{ii} = px;
        frameCoh.py{ii} = py;
        frameIncoh.px{ii} = px;
        frameIncoh.py{ii} = py;

    end
    frameCoh.tt(1:numel(num)) = frameOmega.tt(num);
    frameCoh.Lx = frameOmega.Lx;
    frameCoh.Ly = frameOmega.Ly;
    frameCoh.freq = frameOmega.freq;
    
    frameIncoh.tt(1:numel(num)) = frameOmega.tt(num);
    frameIncoh.Lx = frameOmega.Lx;
    frameIncoh.Ly = frameOmega.Ly;
    frameIncoh.freq = frameOmega.freq; 

%@D
end