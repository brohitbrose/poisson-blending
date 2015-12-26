function [fore_al, mask_al] = alignSource(fore, mask, back)
% im_s2 = alignSource(im_s, mask, im_t);
% After interactively requesting the approximate bottom-center location of
% the paste, returns a pseudo-mask fore_al containing a rectangular crop of
% the foreground image (slightly more generous than mask) on a black
% background, and a true mask mask_al. Both of these are now aligned with
% the back and thus share its dimensions.

figure(1), hold off, imagesc(fore), axis image
figure(2), hold off, imagesc(back), axis image
[y, x] = find(mask);
y1 = min(y)-1; y2 = max(y)+1; x1 = min(x)-1; x2 = max(x)+1;
fore_al = zeros(size(back));
disp(['Select the background location for the bottom-center of the ' ...
    'foreground image']);
[tx, ty] = ginput(1);

yind = (y1:y2);
yind2 = yind - max(y) + round(ty);
xind = (x1:x2);
xind2 = xind - round(mean(x)) + round(tx);

y = y - max(y) + round(ty);
x = x - round(mean(x)) + round(tx);
ind = y + (x-1)*size(back, 1);
mask_al = false(size(back, 1), size(back, 2));
mask_al(ind) = true;

fore_al(yind2, xind2, :) = fore(yind, xind, :);
back(repmat(mask_al, [1 1 3])) = fore_al(repmat(mask_al, [1 1 3]));

figure(1), hold off, imagesc(fore_al), axis image
figure(2), hold off, imagesc(back), axis image
drawnow;

end