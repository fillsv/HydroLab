function frame = flipFrameY(frame)

if exist('frame', 'var')==0 
    disp('function frame = rlipFrameY(frame)');
    return
end
vx = cat(3, frame.vx{:});
vy = cat(3, frame.vy{:});
px = frame.px{1};
py = frame.py{1};
dx = (px(1,2)-px(1,1));
dy = (py(2,1)-py(1,1));
dr = abs(dx+i*dy)/sqrt(2);
vx1 = zeros(size(vx));
vy1 = zeros(size(vx));

px1 = px;
py1 = frame.Ly-py(end:-1:1,:);

% px1 = 
% for ii = 1:size(vx,3)
vx1(:,:,:) = vx(end:-1:1,:,:);
vy1(:,:,:) = -vy(end:-1:1,:,:);    
% end

%
% dx1 = dx*cos(angle)+dy*sin(angle);
% dy1 = -dx*sin(angle)+dy*cos(angle);
% frame.Lx = max(px1(1,:))-min(px1(1,:));
% frame.Ly = max(py1(:,1))-min(py1(:,1));
 
for ii = 1:size(vx,3)
    frame.vx{ii,1} = vx1(:,:,ii);
    frame.vy{ii,1} = vy1(:,:,ii);
    frame.px{ii,1} = px1;
    frame.py{ii,1} = py1;
end