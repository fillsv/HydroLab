function frame = frame3Filter(frame, filterLimit, sigma)

if exist('frame', 'var')==0 
    disp('frame = frame3Filter(frame, filterLimit, sigma)');
    return
end

vx = cat(3, frame.vx{:});
vy = cat(3, frame.vy{:});
if ~exist('sigma','var') sigma = 3; end
if numel(sigma) == 1 sigma = sigma*[1 1 1]; end
if numel(sigma) == 2 sigma = [sigma(1) sigma(1) sigma(2)]; end
vx(isnan(vx)) = 0;
vy(isnan(vy)) = 0;
fvx = imgaussfilt3(vx,sigma);
fvy = imgaussfilt3(vy,sigma);
if ~exist('filterLimit','var') filterLimit = .5; end
if isempty(filterLimit) filterLimit = .5; end
% filterLimit
filterLimitX = mean(abs(fvx),'all')*filterLimit;
filterLimitY = mean(abs(fvy),'all')*filterLimit;

ind = find((abs(fvx-vx)>filterLimitX)|(abs(fvy-vy)>filterLimitY));
vx(ind) = fvx(ind);
vy(ind) = fvy(ind);

for ii = 1:size(vx,3)
    frame.vx{ii,1} = vx(:,:,ii);
    frame.vy{ii,1} = vy(:,:,ii);
end