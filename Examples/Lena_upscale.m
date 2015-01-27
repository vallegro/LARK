% An upscaling example by steering kernel regression
%
% written by hiro, June 17 2007

% load image
img = double(imread('lena.png'));

% downsample the original image
r = 3;  % the downsampling factor
y = downsample2(img, r);

% pilot estimation by classic kernel regression
h = 0.5;
ksize = 5;
[zc, zx1c, zx2c] = ckr2_regular(y, h, 1, ksize);

% obtain orientation information
wsize = 3;   % the size of local analysis window
lambda = 1;  % the regularization for the elongation parameter
alpha = 0.05; % the structure sensitive parameter
C = steering(zx1c, zx2c, ones(size(y)), wsize, lambda, alpha);

% upscaling by steering kernel regression
h = 0.75;
ksize = 7;
[zs, zx1s, zx2s] = skr2_regular(y, h, C, r, ksize);
% root mean square error
error = img - zs(1:512, 1:512);
rmse = sqrt(mymse(error(:)));

% display the estimated image
figure; imagesc(y); colormap(gray); axis image;
title('The downsampled image');
figure; imagesc(zs(1:512, 1:512)); colormap(gray); axis image;
title(['The upscaled image by second order steering kernel regression, RMSE=', num2str(rmse)]);

