function image1 = loadFrameBackground(namein, width, height, frameNumber, cutoff, background)

    f = fopen(namein, 'r');
    fseek(f, height*width*frameNumber, 'bof');
    image1 = fread(f, [width height], 'uint8');
    fclose(f);
    
    if ~exist('cutoff', 'var')
        cutoff = 5e-3;
    end
    if ~exist('background', 'var')
        name_background = [namein(1:end-4) '_background.mat'];
        if exist(name_background, 'file');
            load(name_background)
        else
            background = 0;
        end
    end
%     background = 0;
    image1 = double(image1)-double(background);
%     image1 = image1-min(image1(:));

    minimage1 = min(image1(:));
    maximage1 = max(image1(:));
%     h = hist([minimage1-1; image1(:)], minimage1-1:maximage1);
% h = 1:100;

%     lightmin = max(find(cumsum(h)/sum(h)<cutoff));
%     if isempty(lightmin) lightmin = 0; end
%     lightmax = min(find(cumsum(h)/sum(h)>1-cutoff));
    lightmax = maximage1;
%     image1(find(image1<lightmin)) = lightmin;
%     image1(find(image1>lightmax)) = lightmax;
%     image1 = image1-lightmin;
    image1 = image1/(lightmax)*255;
    image1 = uint8(image1);

