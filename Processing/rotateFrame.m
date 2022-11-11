function frame = rotateFrame(frame, angle)

if exist('frame', 'var')==0 
    disp('function frame = rotateFrame(frame, angle)');
    return
end

vx = cat(3, frame.vx{:});
vy = cat(3, frame.vy{:});
px = frame.px{1};
py = frame.py{1};
dx = (px(1,2)-px(1,1));
dy = (py(2,1)-py(1,1));
dr = abs(dx+i*dy)/sqrt(2);
for ii = 1:size(vx,3)
    vx_new(:,:,ii) = imrotate(vy(:,:,ii), angle);
    vy_new(:,:,ii) = imrotate(vx(:,:,ii), angle);
end
% angle = angle*2*pi/180
% dx1 = dx*cos(angle)+dy*sin(angle);
% dy1 = -dx*sin(angle)+dy*cos(angle);
[px_new py_new] = meshgrid(dr*(1:size(vx_new,2)), dr*(1:size(vx_new,1)) );
 
for ii = 1:size(vx,3)
    frame.vx{ii,1} = vx_new(:,:,ii);
    frame.vy{ii,1} = vy_new(:,:,ii);
    frame.px{ii,1} = px_new;
    frame.py{ii,1} = py_new;
end