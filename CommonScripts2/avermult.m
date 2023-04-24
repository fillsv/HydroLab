function d = avermult(a, mult)

b = reshape(a, mult, size(a,1)/mult, mult, size(a,2)/mult);
c = permute(b, [2,4,1,3]);
c = c(:,:,:);
d = nanmean(c,3);

