function frame = rotateFrame(frame, angle)

if exist('frame', 'var')==0 
    disp('function frame = rotateFrame(frame, angle)');
    return
end
angle = angle/180*pi;
vx = cat(3, frame.vx{:});
vy = cat(3, frame.vy{:});
px = frame.px{1};
py = frame.py{1};
dx = (px(1,2)-px(1,1));
dy = (py(2,1)-py(1,1));
dr = abs(dx+i*dy)/sqrt(2);
vx_new = zeros(size(vx));
vy_new = zeros(size(vx));
px_new = px*cos(angle)-py*sin(angle);
py_new = px*sin(angle)+py*cos(angle);
[px1, py1] = meshgrid(min(px_new(:)):dr:max(px_new(:)), min(py_new(:)):dr:max(py_new(:)) );
vx1 = zeros(size(px));
vy1 = zeros(size(px));
px1_old = px1*cos(-angle)-py1*sin(-angle);
py1_old = px1*sin(-angle)+py1*cos(-angle);

% px1 = 
for ii = 1:size(vx,3)
%     vx_new(:,:,ii) = imrotate(vx(:,:,ii), angle);
%     vy_new(:,:,ii) = imrotate(vy(:,:,ii), angle);
    vx_new(:,:,ii) = vx(:,:,ii)*cos(angle)-vy(:,:,ii)*sin(angle);
    vy_new(:,:,ii) = vx(:,:,ii)*sin(angle)+vy(:,:,ii)*cos(angle);
%     vy_new(:,:,ii) = imrotate(vy(:,:,ii), angle);
    vx1(:,:,ii) = interp2(px, py, vx_new(:,:,ii), px1_old, py1_old);
    vy1(:,:,ii) = interp2(px, py, vy_new(:,:,ii), px1_old, py1_old);    
end

%
% dx1 = dx*cos(angle)+dy*sin(angle);
% dy1 = -dx*sin(angle)+dy*cos(angle);
frame.Lx = max(px1(1,:))-min(px1(1,:));
frame.Ly = max(py1(:,1))-min(py1(:,1));
 
for ii = 1:size(vx,3)
    frame.vx{ii,1} = vx1(:,:,ii);
    frame.vy{ii,1} = vy1(:,:,ii);
    frame.px{ii,1} = px1;
    frame.py{ii,1} = py1;
end