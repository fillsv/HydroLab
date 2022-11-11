function image = loadFrame(namein, width, height, frameNumber)
    if ~exist('frameNumber', 'var')
        disp('function image = loadFrame(namein, width, height, frameNumber)')
        return
    end
    f = fopen(namein, 'r');
    fseek(f, height*width*frameNumber, 'bof');
    image = uint8(fread(f, [width height], 'uint8'));
    fclose(f);