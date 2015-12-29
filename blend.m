function im_blend = blend(fore, mask, background, method)
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
            [A_normal_rows, A_normal_cols, A_normal_vals, b_normal] = ...
                addNormal(A_normal_rows, A_normal_cols, A_normal_vals, ...
                b_normal,norm_index,im2var(i,j),back(i,j));
            norm_index = norm_index + 1;
            % Then check whether neighbors in the positive x- and y-
            % directions lie inside mask_s.
            if ((j+1 <= bw) && mask(i,j+1))
                d = fore(i,j+1)-fore(i,j);
                [A_outer_rows, A_outer_cols, A_outer_vals, b_outer] = ...
                    addOuter(A_outer_rows, A_outer_cols, ... 
                    A_outer_vals, b_outer, outer_index, im2var(i,j+1), ... 
                    back(i,j), d);
                outer_index = outer_index + 1;
            end
            if ((i+1 <= bh) && mask(i+1,j))
                d = fore(i+1,j)-fore(i,j);
                [A_outer_rows, A_outer_cols, A_outer_vals, b_outer] = ...
                    addOuter(A_outer_rows, A_outer_cols, ... 
                    A_outer_vals, b_outer, outer_index, im2var(i+1,j), ...
                    back(i,j), d);
                outer_index = outer_index + 1;
            end
        
        % Otherwise, we must differentiate between (inner) border and
        % foreground pixels in both directions.
        else
            % Handle x-gradient
            if (mask(i,j+1))
                 d = fore(i,j)-fore(i,j+1);
                [A_source_rows, A_source_cols, A_source_vals, b_source] = ...
                    addSource(A_source_rows, A_source_cols, ...
                    A_source_vals,b_source,source_index,im2var(i,j),...
                    im2var(i,j+1),d);
                source_index = source_index + 2;
            else
                d = fore(i,j) - fore(i,j+1);
                [A_border_rows, A_border_cols, A_border_vals, b_border] = ...
                    addBorder(A_border_rows, A_border_cols, ... 
                    A_border_vals,b_border,border_index,im2var(i,j),...
                    back(i,j+1),d);
                border_index = border_index + 1;
            end
            
            % Handle y-gradient
            if (mask(i+1,j))
                d = fore(i,j)-fore(i+1,j);
                [A_source_rows, A_source_cols, A_source_vals, b_source] = ...
                    addSource(A_source_rows, A_source_cols, ...
                    A_source_vals,b_source,source_index,im2var(i,j),...
                    im2var(i+1,j),d);
               source_index = source_index + 2;
            else
                d = fore(i,j) - fore(i+1,j);
                [A_border_rows, A_border_cols, A_border_vals, b_border] = ...
                    addBorder(A_border_rows, A_border_cols, ... 
                    A_border_vals,b_border,border_index,im2var(i,j),...
                    back(i+1,j),d);
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

function [A,b] = get_Ab_mixed(fore, mask, back)
% [A,b] = get_Ab(im_s, mask_s, im_t, method);
% Get the appropriate A and b matrices for least-squares problem Ax=b.
% Details can be found in report.pdf.

% For future convenience, create a mapping from array to linear indices.
[bh, bw] = size(back);
im2var = zeros(bh, bw);
im2var(1:bh*bw) = 1:bh*bw;

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
            [A_normal_rows, A_normal_cols, A_normal_vals, b_normal] = ...
                addNormal(A_normal_rows, A_normal_cols, A_normal_vals, ...
                b_normal,norm_index,im2var(i,j),back(i,j));
            norm_index = norm_index + 1;
            % Then check whether neighbors in the positive x- and y-
            % directions lie inside mask_s.
            if ((j+1 <= bw) && mask(i,j+1))
                d = getMaxMag(fore(i,j), fore(i,j+1), back(i,j), ... 
                    back(i,j+1));
                [A_outer_rows, A_outer_cols, A_outer_vals, b_outer] = ...
                    addOuter(A_outer_rows, A_outer_cols, ... 
                    A_outer_vals, b_outer, outer_index, im2var(i,j+1), ...
                    back(i,j), d);
                outer_index = outer_index + 1;
            end
            if ((i+1 <= bh) && mask(i+1,j))
                d = getMaxMag(fore(i,j), fore(i+1,j), back(i,j), ...
                    back(i+1,j));
                [A_outer_rows, A_outer_cols, A_outer_vals, b_outer] = ...
                    addOuter(A_outer_rows, A_outer_cols, ... 
                    A_outer_vals, b_outer, outer_index, im2var(i+1,j), ...
                    back(i,j), d);
                outer_index = outer_index + 1;
            end
        
        % Otherwise, we must differentiate between (inner) border and
        % foreground pixels in both directions.
        else
            % Handle x-gradient
            if (mask(i,j+1))
                 d = getMaxMag(fore(i,j+1), fore(i,j), back(i,j+1),...
                     back(i,j));
                [A_source_rows, A_source_cols, A_source_vals, b_source] = ...
                    addSource(A_source_rows, A_source_cols, ...
                    A_source_vals,b_source,source_index,im2var(i,j),...
                    im2var(i,j+1),d);
                source_index = source_index + 2;
            else
                d = getMaxMag(fore(i,j+1), fore(i,j), back(i,j+1), ...
                    back(i,j));
                [A_border_rows, A_border_cols, A_border_vals, b_border] = ...
                    addBorder(A_border_rows, A_border_cols, ... 
                    A_border_vals,b_border,border_index,im2var(i,j),...
                    back(i,j+1),d);
                border_index = border_index + 1;
            end
            
            % Handle y-gradient
            if (mask(i+1,j))
                d = getMaxMag(fore(i+1,j), fore(i,j), back(i+1,j), ...
                    back(i,j));
                [A_source_rows, A_source_cols, A_source_vals, b_source] = ...
                    addSource(A_source_rows, A_source_cols, ...
                    A_source_vals,b_source,source_index,im2var(i,j),im2var(i+1,j),d);
               source_index = source_index + 2;
            else
                d = getMaxMag(fore(i+1,j),fore(i,j), back(i+1,j), ...
                    back(i,j));
                [A_border_rows, A_border_cols, A_border_vals, b_border] = ...
                    addBorder(A_border_rows, A_border_cols, ... 
                    A_border_vals,b_border,border_index,im2var(i,j),back(i+1,j),d);
                border_index = border_index + 1; 
            end
        end
    end
end

% Strip trailing zeros.
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

function [nr, nc, nv, nb] = addNormal(nr,nc,nv,nb,ni,col,bval)

nr(ni) = ni;
nc(ni) = col;
nv(ni) = 1;
nb(ni) = bval;

end

function [or, oc, ov, ob] = addOuter(or,oc,ov,ob,oi,col,bval,d)

or(oi) = oi;
oc(oi) = col;
ov(oi) = 1;
ob(oi) = bval+d;

end

function [sr, sc, sv, sb] = addSource(sr,sc,sv,sb,si,col,col1,d)

t = (si+1) ./ 2;
sr(si) = t;
sc(si) = col;
sv(si) = 1;
sb(t) = d;
si1 = si+1;
sr(si1) = sr(si);
sc(si1) = col1;
sv(si1) = -1;

end

function [br,bc,bv,bb] = addBorder(br,bc,bv,bb,bi,col,bval,d)

br(bi) = bi;
bc(bi) = col;
bv(bi) = 1;
bb(bi) = bval+d;

end

function d = getMaxMag(a1,a2,b1,b2)

if abs(a2-a1) > abs(b2-b1)
    d = a2-a1;
else
    d = b2-b1;
end

end