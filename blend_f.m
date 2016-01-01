function im_blend = blend_f(fore, mask, background, method)
% im_blend = blend(im_s, mask_s, im_background, method);
% Merges im_s with im_background in accordance with mask_s into im_blend.
% The nature of the boundary conditions is determined by method.

% Some error identification...
if ~(strcmp(method, 'poisson') || strcmp(method, 'mixed'))
    error('Invalid blending method specified');
end
[bh, bw, bcolor] = size(background);
[~, ~, imcolor] = size(fore);
if ~(bcolor == imcolor)
    error('Blending color with grayscale images currently unsupported');
end

% If color images, perform the logic for each color channel. Otherwise,
% performing it once on the grayscales will suffice.
if (bcolor == 3)

    % Solve a least-squares problem Ax=b for each of the color channels.
    back_red = background(:,:,1);
    fore_red = fore(:,:,1);
    back_green = background(:,:,2);
    fore_green = fore(:,:,2);
    back_blue = background(:,:,3);
    fore_blue = fore(:,:,3);
    
    if (strcmp(method, 'poisson'))
        [A_red, b_red] = get_Ab_poisson(fore_red, mask, back_red);
        red_out = reshape(A_red\b_red, [bh bw]);
        [A_green, b_green] = get_Ab_poisson(fore_green, mask, back_green);
        green_out = reshape(A_green\b_green, [bh bw]);    
        [A_blue, b_blue] = get_Ab_poisson(fore_blue, mask, back_blue);
        blue_out = reshape(A_blue\b_blue, [bh bw]);
    elseif (strcmp(method, 'mixed'))
        [A_red, b_red] = get_Ab_mixed(fore_red, mask, back_red);
        red_out = reshape(A_red\b_red, [bh bw]);
        [A_green, b_green] = get_Ab_mixed(fore_green, mask, back_green);
        green_out = reshape(A_green\b_green, [bh bw]);    
        [A_blue, b_blue] = get_Ab_mixed(fore_blue, mask, back_blue);
        blue_out = reshape(A_blue\b_blue, [bh bw]);
    end
    
    im_blend = cat(3, red_out, green_out, blue_out);

else
    [A, b] = get_Ab(fore, mask, background, method);
    im_blend = reshape(A\b, [bh bw]);
end

end

function [A,b] = get_Ab_poisson(fore, mask, back)
% [A,b] = get_Ab(im_s, mask_s, im_t, method);
% Get the appropriate A and b matrices for least-squares problem Ax=b.
% Details can be found in report.pdf.

% For future convenience, create a mapping from array to linear indices.
[bh, bw] = size(back);
im2var = zeros(bh, bw);
im2var(1:bh*bw) = 1:bh*bw;

% The resulting A matrix will be extremely sparse. To preserve memory while
% optimizing speed, we determine the indices and values of nonzero A
% entries in advance and set them with a single call to sparse().

% Equations that directly copy background pixels to output pixels.
A_normal_rows = zeros(1,bh*bw);
A_normal_cols = zeros(1,bh*bw);
A_normal_vals = zeros(1,bh*bw);
b_normal = zeros(1,bh*bw);
norm_index = 1;

% Equations corresponding to pixels inside mask_s with a right- or down-
% neighbor also inside mask_s.
A_source_rows = zeros(1,2*bh*bw);
A_source_cols = zeros(1,2*bh*bw);
A_source_vals = zeros(1,2*bh*bw);
b_source = zeros(1,bh*bw);
source_index = 1;

% % Equations corresponding to pixels inside mask_s with a right- or down-
% % neighbor outside mask_s.
% A_border_rows = zeros(1,bh*bw);
% A_border_cols = zeros(1,bh*bw);
% A_border_vals = zeros(1,bh*bw);
% b_border = zeros(1,bh*bw);
% border_index = 1;

% Equations corresponding to pixels inside mask_s with an up- or a left-
% neighbor outside mask_s.
A_outer_rows = zeros(1,bh*bw);
A_outer_cols = zeros(1,bh*bw);
A_outer_vals = zeros(1,bh*bw);
b_outer = zeros(1,bh*bw);
outer_index = 1;

for j=1:bw
    for i=1:bh
        if (~mask(i,j))
            if ((j+1) <= bw)
                if (~mask(i,j+1))
                    t = (source_index + 1) ./ 2;
                    A_source_rows(source_index) = t;
                    A_source_cols(source_index) = im2var(i,j);
                    A_source_vals(source_index) = 1;
                    d = back(i,j)-back(i,j+1);
                    b_source(t) = d;
                    t = source_index + 1;               
                    A_source_rows(t) = A_source_rows(source_index);
                    A_source_cols(t) = im2var(i,j+1);
                    A_source_vals(t) = -1;
                    source_index = t + 1;
                else
                    A_outer_rows(outer_index) = outer_index;
                    A_outer_cols(outer_index) = im2var(i,j);
                    A_outer_vals(outer_index) = 1;
                    d = back(i,j)-back(i,j+1);
                    b_outer(outer_index) = d + fore(i,j+1);
                    outer_index = outer_index + 1;     
                end
            end
            if ((i+1) <= bh)
                if (~mask(i+1,j))
                    t = (source_index + 1) ./ 2;
                    A_source_rows(source_index) = t;
                    A_source_cols(source_index) = im2var(i,j);
                    A_source_vals(source_index) = 1;
                    d = back(i,j)-back(i+1,j);
                    b_source(t) = d;
                    t = source_index + 1;
                    A_source_rows(t) = A_source_rows(source_index);
                    A_source_cols(t) = im2var(i+1,j);
                    A_source_vals(t) = -1;
                    source_index = t + 1;
                else
                    A_outer_rows(outer_index) = outer_index;
                    A_outer_cols(outer_index) = im2var(i,j);
                    A_outer_vals(outer_index) = 1;
                    d = back(i,j)-back(i+1,j);
                    b_outer(outer_index) = d + fore(i+1,j);
                    outer_index = outer_index + 1;
                end
            end
        else
            A_normal_rows(norm_index) = norm_index;
            A_normal_cols(norm_index) = im2var(i,j);
            A_normal_vals(norm_index) = 1;
            b_normal(norm_index) = fore(i,j);
            norm_index = norm_index + 1;
            % Then check whether neighbors in the positive x- and y-
            % directions lie inside mask_s.
            if ((j+1 <= bw) && ~mask(i,j+1))
                A_outer_rows(outer_index) = outer_index;
                A_outer_cols(outer_index) = im2var(i,j+1);
                A_outer_vals(outer_index) = 1;
                d = back(i,j+1)-back(i,j);
                b_outer(outer_index) = d + fore(i,j);
                outer_index = outer_index + 1;
            end
            if ((i+1 <= bh) && ~mask(i+1,j))
                A_outer_rows(outer_index) = outer_index;
                A_outer_cols(outer_index) = im2var(i+1,j);
                A_outer_vals(outer_index) = 1;
                d = back(i+1,j)-back(i,j);
                b_outer(outer_index) = d + fore(i,j);
                outer_index = outer_index + 1;
            end
        end
    end
end

% Having generously initialized arrays in order to prevent resizing inside
% memory, we now strip trailing zeros.
n = norm_index-1;
A_normal_rows = A_normal_rows(1:n);
A_normal_cols = A_normal_cols(1:n);
A_normal_vals = A_normal_vals(1:n);
A_normal = sparse(A_normal_rows, A_normal_cols, A_normal_vals, ...
    n, bh*bw);
b_normal = b_normal(1:n);

s = source_index-1;
A_source_rows = A_source_rows(1:s);
A_source_cols = A_source_cols(1:s);
A_source_vals = A_source_vals(1:s);
t = s./2;
A_source = sparse(A_source_rows, A_source_cols, A_source_vals, ...
    t, bh*bw);
b_source = b_source(1:t);

% bi = border_index-1;
% A_border_rows = A_border_rows(1:bi);
% A_border_cols = A_border_cols(1:bi);
% A_border_vals = A_border_vals(1:bi);
% A_border = sparse(A_border_rows, A_border_cols, A_border_vals, ...
%     bi, bh*bw);
% b_border = b_border(1:bi);

o = outer_index-1;
A_outer_rows = A_outer_rows(1:o);
A_outer_cols = A_outer_cols(1:o);
A_outer_vals = A_outer_vals(1:o);
A_outer = sparse(A_outer_rows, A_outer_cols, A_outer_vals, ...
    o, bh*bw);
b_outer = b_outer(1:o);

A = cat(1,A_normal,A_source,A_outer);
b = cat(2,b_normal,b_source,b_outer);
b = b';

end


% INCOMPLETE 
function [A,b] = get_Ab_mixed(fore, mask, back)
% [A,b] = get_Ab(im_s, mask_s, im_t, method);
% Get the appropriate A and b matrices for least-squares problem Ax=b.
% Details can be found in report.pdf.

% For future convenience, create a mapping from array to linear indices.
[bh, bw] = size(back);
im2var = zeros(bh, bw);
im2var(1:bh*bw) = 1:bh*bw;

% The resulting A matrix will be extremely sparse. To preserve memory while
% optimizing speed, we determine the indices and values of nonzero A
% entries in advance and set them with a single call to sparse().

% Equations that directly copy background pixels to output pixels.
A_normal_rows = zeros(1,bh*bw);
A_normal_cols = zeros(1,bh*bw);
A_normal_vals = zeros(1,bh*bw);
b_normal = zeros(1,bh*bw);
norm_index = 1;

% Equations corresponding to pixels inside mask_s with a right- or down-
% neighbor also inside mask_s.
A_source_rows = zeros(1,2*bh*bw);
A_source_cols = zeros(1,2*bh*bw);
A_source_vals = zeros(1,2*bh*bw);
b_source = zeros(1,bh*bw);
source_index = 1;

% Equations corresponding to pixels inside mask_s with a right- or down-
% neighbor outside mask_s.
A_border_rows = zeros(1,bh*bw);
A_border_cols = zeros(1,bh*bw);
A_border_vals = zeros(1,bh*bw);
b_border = zeros(1,bh*bw);
border_index = 1;

% Equations corresponding to pixels inside mask_s with an up- or a left-
% neighbor outside mask_s.
A_outer_rows = zeros(1,bh*bw);
A_outer_cols = zeros(1,bh*bw);
A_outer_vals = zeros(1,bh*bw);
b_outer = zeros(1,bh*bw);
outer_index = 1;

for j=1:bw
    for i=1:bh
        % If outside mask_s, set final values to background values.
        if (~mask(i,j))
            A_normal_rows(norm_index) = norm_index;
            A_normal_cols(norm_index) = im2var(i,j);
            A_normal_vals(norm_index) = 1;
            b_normal(norm_index) = back(i,j);
            norm_index = norm_index + 1;
            % Then check whether neighbors in the positive x- and y-
            % directions lie inside mask_s.
            if ((j+1 <= bw) && mask(i,j+1))
                A_outer_rows(outer_index) = outer_index;
                A_outer_cols(outer_index) = im2var(i,j+1);
                A_outer_vals(outer_index) = 1;
                d = getMaxMag(fore(i,j), fore(i,j+1), back(i,j), ...
                    back(i,j+1));
                b_outer(outer_index) = d + back(i,j);
                outer_index = outer_index + 1;
            end
            if ((i+1 <= bh) && mask(i+1,j))
                A_outer_rows(outer_index) = outer_index;
                A_outer_cols(outer_index) = im2var(i+1,j);
                A_outer_vals(outer_index) = 1;
                d = getMaxMag(fore(i,j), fore(i+1,j), back(i,j), ...
                    back(i+1,j)); 
                b_outer(outer_index) = d + back(i,j);
                outer_index = outer_index + 1;
            end
        
        % Otherwise, we must differentiate between (inner) border and
        % foreground pixels in both directions.
        else
            % Handle x-gradient
            if (mask(i,j+1))
                t = (source_index + 1) ./ 2;
                A_source_rows(source_index) = t;
                A_source_cols(source_index) = im2var(i,j);
                A_source_vals(source_index) = 1;
                d = getMaxMag(fore(i,j+1),fore(i,j), ...
                    back(i,j+1),back(i,j));
                b_source(t) = d;
                t = source_index + 1;               
                A_source_rows(t) = A_source_rows(source_index);
                A_source_cols(t) = im2var(i,j+1);
                A_source_vals(t) = -1;
                source_index = t + 1;
            else
                A_border_rows(border_index) = border_index;
                A_border_cols(border_index) = im2var(i,j);
                A_border_vals(border_index) = 1;
                d = getMaxMag(fore(i,j+1),fore(i,j), ...
                    back(i,j+1),back(i,j));
                b_border(border_index) = d+back(i,j+1);
                border_index = border_index + 1;
            end
            
            % Handle y-gradient
            if (mask(i+1,j))
                t = (source_index + 1) ./ 2;
                A_source_rows(source_index) = t;
                A_source_cols(source_index) = im2var(i,j);
                A_source_vals(source_index) = 1;
                d = getMaxMag(fore(i+1,j),fore(i,j), ...
                    back(i+1,j),back(i,j));
                b_source(t) = d;
                t = source_index + 1;
                A_source_rows(t) = A_source_rows(source_index);
                A_source_cols(t) = im2var(i+1,j);
                A_source_vals(t) = -1;
                source_index = t + 1;
            else
                A_border_rows(border_index) = border_index;
                A_border_cols(border_index) = im2var(i,j);
                A_border_vals(border_index) = 1;
                d = getMaxMag(fore(i+1,j),fore(i,j), ...
                    back(i+1,j),back(i,j));
                b_border(border_index) = d+back(i+1,j);
                border_index = border_index + 1; 
            end
        end
    end
end

% Having generously initialized arrays in order to prevent resizing inside
% memory, we now strip trailing zeros.
n = norm_index-1;
A_normal_rows = A_normal_rows(1:n);
A_normal_cols = A_normal_cols(1:n);
A_normal_vals = A_normal_vals(1:n);
A_normal = sparse(A_normal_rows, A_normal_cols, A_normal_vals, ...
    n, bh*bw);
b_normal = b_normal(1:n);

s = source_index-1;
A_source_rows = A_source_rows(1:s);
A_source_cols = A_source_cols(1:s);
A_source_vals = A_source_vals(1:s);
t = s./2;
A_source = sparse(A_source_rows, A_source_cols, A_source_vals, ...
    t, bh*bw);
b_source = b_source(1:t);

bi = border_index-1;
A_border_rows = A_border_rows(1:bi);
A_border_cols = A_border_cols(1:bi);
A_border_vals = A_border_vals(1:bi);
A_border = sparse(A_border_rows, A_border_cols, A_border_vals, ...
    bi, bh*bw);
b_border = b_border(1:bi);

o = outer_index-1;
A_outer_rows = A_outer_rows(1:o);
A_outer_cols = A_outer_cols(1:o);
A_outer_vals = A_outer_vals(1:o);
A_outer = sparse(A_outer_rows, A_outer_cols, A_outer_vals, ...
    o, bh*bw);
b_outer = b_outer(1:o);

A = cat(1,A_normal,A_source,A_border,A_outer);
b = cat(2,b_normal,b_source,b_border,b_outer);
b = b';

end

function d = getMaxMag(a1,a2,b1,b2)
d1 = a2 - a1;
d2 = b2 - b1;
if abs(d1) > abs(d2)
    d = d1;
else
    d = d2;
end

end