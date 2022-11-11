function frameFiltred = filterFrame(frame, num)

    if ~exist('frame', 'var')
        disp('function filterFrame(frame, num)');
        return
    end
    
    if exist('num', 'var')==0 
        num = 1:numel(frame.px);
    end    
    kk = 0;
    frameFiltred = frame;
    frameFiltred.px = {};
    frameFiltred.py = {};
    frameFiltred.vx = {};
    frameFiltred.vy = {};        
    for jj = num
        frameFiltred.px{end+1} = frame.px{jj};
        frameFiltred.py{end+1} = frame.py{jj};
%        sigma = 5;
%        frameFiltred.vx{end+1} = imgaussfilt(frame.vx{jj}, sigma, 'padding', 'symmetric');    
%        frameFiltred.vy{end+1} = imgaussfilt(frame.vy{jj}, sigma, 'padding', 'symmetric');    
         [frameFiltred.vx{end+1} frameFiltred.vy{end+1}] = ...
             filterV(frame.vx{jj}, frame.vy{jj}, 3.5);
    end
    
    

        
        
        
        
        
        