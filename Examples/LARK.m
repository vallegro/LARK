% A denoising example by iterative steering kernel regression
%
% written by hiro, June 17 2007

% load image
%img = double(imread('TestImages/lena.png'));
[N,M] = size(y);

% add white Gaussian noise
%sigma = 25;        % standard deviation
%randn('state', 0); % initialization
%y = round0_255(img + randn(size(img)) * sigma);

% pilot estimation by second order classic kernel regression
h = 0.5;    % the global smoothing parameter
r = 1;      % the upscaling factor
ksize = 5;  % the kernel size
[zc, zx1c, zx2c] = ckr2_regular(y, h, r, ksize);

% iteartive steering kernel regression (second order)
IT = 30;     % the total number of iterations
wsize = 11;  % the size of the local orientation analysis window
lambda = 1;  % the regularization for the elongation parameter
alpha = 0.5; % the structure sensitive parameter
h = 2.4;     % the global smoothing parameter
ksize = 17;  % the kernel size
z = zeros(N, M, IT+1);
zx1 = zeros(N, M, IT+1);
zx2 = zeros(N, M, IT+1);
rmse = zeros(IT+1, 1);
z(:,:,1) = y;
zx1(:,:,1) = zx1c;
zx2(:,:,1) = zx2c;
error = img - y;
rmse(1) = sqrt(mymse(error(:)));

for i = 2 : IT+1
    % compute steering matrix
    C = steering(zx1(:,:,i-1), zx2(:,:,i-1), ones(size(img)), wsize, lambda, alpha);
    % steering kernel regression
    [zs, zx1s, zx2s] = skr2_regular(z(:,:,i-1), h, C, r, ksize);
    z(:,:,i) = zs;
    %imtool(uint8(zs));
    zx1(:,:,i) = zx1s;
    zx2(:,:,i) = zx2s;
    % root mean square error
    error = img - zs;
    rmse(i) = sqrt(mymse(error(:)));
    rmse(i)
    %figure(99); imagesc(zs); colormap(gray); axis image; pause(1);
end

% display images
figure; imagesc(y); colormap(gray); axis image;
title(['The noisy image, STD=25, RMSE=', num2str(rmse(1))]);
figure; imagesc(z(:,:,13)); colormap(gray); axis image;
title(['ISKR, 12 iterations, RMSE=', num2str(rmse(13))]);