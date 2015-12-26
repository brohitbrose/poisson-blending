% Master script for using and testing the image blending algorithms.

RECONSTRUCT = false;
SEAMLESS = false;
MIXED  = true;

if RECONSTRUCT 
    % Display original image for reference.
    toyim = im2double(imread('./samples/dooba.jpg'));
    imshow(toyim);
    title('Original')
    
    % Create x- and y-gradients for this image. Reconstruct the image
    % using these and one original pixel value alone.
    im_out = toy_reconstruct(toyim);
    
    % Display reconstructed image as well as its accumulated error from the
    % original image. There exists a bug if the error isn't near 0.
    figure(9);
    imshow(im_out);
    title('Reconstructed');
    disp(['Error: ' num2str(sqrt(sum((toyim(:)-im_out(:)).^2)))])
end

if SEAMLESS
    % Roughly resize images so that desired overlapping regions are of
    % approximately the same size.
    background = im2double(imread('./samples/monalisa.jpg'));
    foreground = im2double(imread('./samples/john-cena.jpg'));
    [background, foreground] = resizeImage(background, foreground);
    close all;
    
    % Interactively create a binary mask of the foreground image. True only
    % inside the user-inputted polygon.
    disp(['Choose vertices corresponding to the polygon containing ' ...
        'foreground region of interest.']);
    objmask = roipoly(foreground);
    
    % See alignSource.m for details regarding this function.
     [fore, mask] = alignSource(foreground, objmask, background);

    % Perform seamless cloning of foreground onto background.
    im_blend = blend(fore, mask, background, 'poisson');
    figure(10), hold off, imshow(im_blend)
end

if MIXED
    % Roughly resize images so that overlapping regions are of
    % approximately the same size.
    background = im2double(imread('./samples/monalisa.jpg'));
    foreground = im2double(imread('./samples/john-cena.jpg'));
    [background, foreground] = resizeImage(background, foreground);
    close all;
    
    % Interactively create a binary mask of the foreground image. True only
    % inside the user-inputted polygon.
    disp(['Choose vertices corresponding to the polygon containing ' ...
        'foreground region of interest.']);
    objmask = roipoly(foreground);
    
    % See alignSource.m for details regarding this call.
    [fore, mask] = alignSource(foreground, objmask, background);

    % Perform mixed seamless cloning of foreground onto background.
    im_blend = blend(fore, mask, background, 'mixed');
    
    % Set to true if only the blended region is desired.
    texture = false;
    if (texture)
        im_blend = cat(3, im_blend(:,:,1) .* mask, ...
            im_blend(:,:,2) .* mask, ...
            im_blend(:,:,3) .* mask);
    end
    
    figure(11), hold off, imshow(im_blend);
end
