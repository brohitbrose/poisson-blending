close all;
clear;
clc;

% Master script for using and testing the image blending algorithms.

% Set the value corresponding to the desired test(s) to true.
RECONSTRUCT = false;
SEAMLESS = true;
MIXED  = false;

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
    tic
    im_blend = blend(fore, mask, background, 'poisson');
    toc
    
    % Repeated testing showed that argument passing involving large
    % matrices in MATLAB has nontrivial overhead. This function is now
    % obsolete performance-wise because blend.m will always blend images
    % more quickly. However, blend.m is written in terrible style, with
    % several dozens of lines directly copy/pasted between the various
    % "get_Ab"s. Uncomment the lines below to test the modularized
    % function's performace. Ideally, we would want to keep modularized
    % code that somehow circumvents the argument passing overhead.
    
    % tic
    % im_blend = blend_modularized(fore, mask, background, 'poisson');
    % toc
    figure(10), hold off, imshow(im_blend)
end

if MIXED
    % Roughly resize images so that overlapping regions are of
    % approximately the same size.
    background = im2double(imread('./samples/wood.jpg'));
    foreground = im2double(imread('./samples/charizard.jpg'));
    [background, foreground] = resizeImage(background, foreground);
    close all;
    
    % Interactively create a binary mask of the foreground image. True only
    % inside the user-inputted polygon.
    % TODO: due to the nature of mixed-blending images, give an option for
    % rectangular selection in addition to the default free-form select.
    disp(['Choose vertices corresponding to the polygon containing ' ...
        'foreground region of interest.']);
    objmask = roipoly(foreground);
    
    % See alignSource.m for details regarding this call.
    [fore, mask] = alignSource(foreground, objmask, background);

    % Perform mixed seamless cloning of foreground onto background.
    tic
    im_blend = blend(fore, mask, background, 'mixed');
    toc
    
    % Repeated testing showed that argument passing involving large
    % matrices in MATLAB has nontrivial overhead. This function is now
    % obsolete performance-wise because blend.m will always blend images
    % more quickly. However, blend.m is written in terrible style, with
    % several dozens of lines directly copy/pasted between the various
    % "get_Ab"s. Uncomment the lines below to test the modularized
    % function's performace. Ideally, we would want to keep modularized
    % code that somehow circumvents the argument passing overhead.
    
    % tic
    % im_blend = blend_modularized(fore, mask, background, 'mixed');
    % toc
    
    % Set to true if only the blended region is desired.
    texture = false;
    if (texture)
        im_blend = cat(3, im_blend(:,:,1) .* mask, ...
            im_blend(:,:,2) .* mask, ...
            im_blend(:,:,3) .* mask);
    end
    
    figure(11), hold off, imshow(im_blend);
end
