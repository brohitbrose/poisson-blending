function im_out = reconstruct(im)
% im_out = toy_reconstruct(im)
% Constructs a basic "gradient" from input image, computed by differencing
% direct neighbors along x- and y- directions. Reconstructs image from this
% "gradient," and then returns the reconstructed image. Details are present
% in report.pdf.

[imh, imw, color] = size(im);

% If im is a color image, perform the logic for each color channel.
if (color == 3)
    im_red = im(:,:,1);
    [A_red, b_red] = get_Ab(im_red);
    red_out = reshape(A_red\b_red, [imh imw]);
    
    im_green = im(:,:,2);
    [A_green, b_green] = get_Ab(im_green);
    green_out = reshape(A_green\b_green, [imh imw]);
    
    im_blue = im(:,:,3);
    [A_blue, b_blue] = get_Ab(im_blue);
    blue_out = reshape(A_blue\b_blue, [imh imw]);
    
    im_out = cat(3, red_out, green_out, blue_out);
    
else
    [A, b] = get_Ab(im);
    im_out = reshape(A\b, [imh imw]);
end

end


function [A, b] = get_Ab(im)

[imh, imw] = size(im);
im2var = zeros(imh, imw);
im2var(1:imh*imw) = 1:imh*imw;

% Least-squares objective: "x-gradients"
A1 = get_A1(imh, imw);
b1 = zeros(imh*(imw-1),1);
for j = 1:(imw-1)
    for i = 1:imh
        b1(im2var(i,j),1) = im(i,j+1)-im(i,j);
    end
end

xgrad = reshape(b1, [imh imw-1]);
figure();
imshow(mat2gray(xgrad));
title('Single color x-gradient');

% Least-squares objective: "y-gradients"
A2 = get_A2(imh, imw);
b2 = zeros((imh-1)*imw,1);
d = 1;
for j = 1:imw
    for i = 1:imh-1
        b2(d,1) = im(i+1,j)-im(i,j);
        d = d+1;
    end
end

ygrad = reshape(b2, [imh-1 imw]);
figure();
imshow(mat2gray(ygrad));
title('Single color y-gradient');

% Least squares objective: single value for reconstruction
A3 = sparse(1, 1, 1, 1, imh*imw);
b3 = im(1,1);

% Combine the three above into one system...
A = cat(1,A1,A2,A3);
b = cat(1,b1,b2,b3);

end

% TODO: explain why this works
function A1 = get_A1(imh, imw)

rel = imh*(imw-1);
A1_rows = kron(1:rel,[1 1]);
A1_columns = zeros(rel,2);
A1_columns(:,1) = (1:rel)';
A1_columns(:,2) = ((1:rel)+imh)';
A1_columns = reshape(A1_columns',[1 2*rel]);
A1_vals = 1:2*rel;
A1_vals = -round(2 * (mod(A1_vals, 2) - 1/2));
A1 = sparse(A1_rows, A1_columns, A1_vals, rel, imw*imh);

end

% TODO: explain why 4AM coding is an excellent idea
function A2 = get_A2(imh, imw)

rel = (imh-1)*imw;
A2_rows = kron(1:rel,[1 1]);
A2_columns = zeros(imw,2*(imh-1));
for i = 1:imw
    A2_columns = A2set(A2_columns, i, 2*(imh-1), 1 + imh*(i-1));
end
A2_columns = reshape(A2_columns',[1 2*rel]);
A2_vals = 1:2*rel;
A2_vals = -round(2 * (mod(A2_vals, 2) - 1/2)); % lol
A2 = sparse(A2_rows, A2_columns, A2_vals, rel, imw*imh);

end

function A2_columns = A2set(A2_columns, r, len, ival)
A2_columns(r,1) = ival;
v = 1;
for i = 2:2:len-1
    A2_columns(r,i) = ival+v;
    A2_columns(r,i+1) = ival+v;
    v = v+1;
end
A2_columns(r,len) = ival+v;
end