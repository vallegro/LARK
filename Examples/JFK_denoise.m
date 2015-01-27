% A film grain removal example by iterative steering kernel regression
%
% written by hiro, June 17 2007

% load image
img = double(imread('JFKreg.jpg'));
[N,M] = size(img(:,:,1));

% convert the image into the YCrCb channels
[Y, Cr, Cb] = RGB2YCC(img);

% the input image
y = zeros(N, M, 3);
y(:,:,1) = Y;
y(:,:,2) = Cr;
y(:,:,3) = Cb;

% pilot estimation by second order classic kernel regression
% on the luminance channel
h = 0.5;    % the global smoothing parameter
r = 1;      % the upscaling factor
ksize = 5;  % the kernel size
[zc, zx1c, zx2c] = ckr2_regular(y(:,:,1), h, r, ksize);

% iteartive steering kernel regression (second order)
IT = 5;     % the total number of iterations
wsize = 11;  % the size of the local orientation analysis window
lambda = 1;  % the regularization for the elongation parameter
alpha = 0.5; % the structure sensitive parameter
h = 2.3;     % the global smoothing parameter
ksize = 11;  % the kernel size
z = zeros(N, M, 3, IT+1);
zx1 = zeros(N, M, 3, IT+1);
zx2 = zeros(N, M, 3, IT+1);
z(:,:,:,1) = y;
zx1(:,:,1) = zx1c;
zx2(:,:,1) = zx2c;

for i = 2 : IT+1
    % compute steering matrix from the Y channel
    C = steering(zx1(:,:,i-1), zx2(:,:,i-1), ones(size(img)), wsize, lambda, alpha);
    % steering kernel regression
    % Y channel
    [zs, zx1s, zx2s] = skr2_regular(z(:,:,1,i-1), h, C, r, ksize);
    z(:,:,1,i) = zs;
    zx1(:,:,i) = zx1s;
    zx2(:,:,i) = zx2s;
    % Cr channel
    [zs, zx1s, zx2s] = skr2_regular(z(:,:,2,i-1), h, C, r, ksize);
    z(:,:,2,i) = zs;
    % Cb channel
    [zs, zx1s, zx2s] = skr2_regular(z(:,:,3,i-1), h, C, r, ksize);
    z(:,:,3,i) = zs;
end

% display images
zRGB = YCC2RGB(z(:,:,1,4), z(:,:,2,4), z(:,:,3,4));
figure; imagesc(uint8(img)); colormap(gray); axis image;
title('The original image');
figure; imagesc(uint8(zRGB)); colormap(gray); axis image;
title('The denoised image by iterative steering kernel regression, 3 iterations');
gray_reverse = gray;
gray_reverse = gray_reverse(64:-1:1, :);
figure; imagesc(abs(z(:,:,1,1) - z(:,:,1,4))); colormap(gray_reverse); axis image; colorbar;
title('The absolute residual image in the luminance channel');

