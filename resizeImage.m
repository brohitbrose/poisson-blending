function [im_b, im_o] = resizeImage(im_b, im_o)

disp(['HOLD SHIFT and select a (square) region approximately ' ...
    'corresponding to desired BACKGROUND region.']);
figure(1);
imshow(im_b);
rect_b = getrect(1);
if ~(abs(rect_b(3) - rect_b(4)) < 5)
    error('Please HOLD SHIFT when selecting region.');
end
disp(['HOLD SHIFT and select a (square) region approximately ' ...
    'corresponding to desired FOREGROUND region.']);
figure(2);
imshow(im_o);
rect_o = getrect(2);
if ~(abs(rect_o(3) - rect_o(4)) < 5)
    error('Please HOLD SHIFT when selecting region.');
end
if rect_b(3) < rect_o(3)
    im_o = imresize(im_o, rect_b(3)/rect_o(3), 'bilinear'); 
elseif rect_o(3) < rect_b(3)
    im_b = imresize(im_b, rect_o(3)/rect_b(3), 'bilinear');
end

end