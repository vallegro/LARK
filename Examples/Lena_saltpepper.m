% A salt & pepper noise reduction example by 
% steering kernel regression
%
% written by hiro, June 17 2007

% load image
img = double(imread('lena.png'));
img = img(171:290, 171:290);
[N,M] = size(img);

% add salt & pepper noise
noise_density = 20; % [percent]
y = double(imnoise(uint8(img),'salt & pepper', noise_density/100));

% Initial estimation
% apply 3x3 median filter
z_med = medfilt2(y, [3, 3]);
error = img(6:end-5, 6:end-5) - z_med(6:end-5, 6:end-5);
rmse_median = sqrt(mymse(error(:)));
% estimate the derivatives by L2 classic kernel regression
h = 0.8;    % the global smoothing parameter
r = 1;      % the upscaling factor
ksize = 7;  % the support size of the kernel function
zcL2 = ckr2all_regular(z_med, h, r, ksize); % "zxL2" is the initial estimate

% remove salt & pepper noise by L1 classic kernel regression
z_init = zcL2; % initial state
h = 1.5;       % the global smoothing parameter
r = 1;         % the upscaling factor
ksize = 9;     % the support size of the kernel function
IT = 750;      % the number of iteration of the steepest descent method
step = 0.1;    % the step size for the steepest descent update
zcL1 = ckr2L1_regular(y, z_init, h, r, ksize, IT, step);
% compute root mean square error
error = img(6:end-5, 6:end-5) - zcL1(6:end-5, 6:end-5, 1);
rmse_ckrL1 = sqrt(mymse(error(:)));

% obtain orientation information
wsize = 9;   % the size of local analysis window
lambda = 1;  % the regularization for the elongation parameter
alpha = 0.0; % the structure sensitive parameter
C = steering(zcL1(:,:,2), zcL1(:,:,3), ones(size(y)), wsize, lambda, alpha);

% remove salt & pepper noise by L1 steering kernel regression
z_init = zcL1; % initial state
h = 1.75;      % the global smoothing parameter
r = 1;         % the upscaling factor
ksize = 9;     % the support size of the kernel function
IT = 300;      % the number of iteration of the steepest descent method
step = 0.1;    % the step size for the steepest descent update
zsL1 = skr2L1_regular(y, z_init, h, C, r, ksize, IT, step);
% compute root mean square error
error = img(6:end-5, 6:end-5) - zsL1(6:end-5, 6:end-5, 1);
rmse_skrL1 = sqrt(mymse(error(:)));

% display the result
figure; imagesc(y); colormap(gray); axis image;
title('Salt & Pepper noise (20%)');
figure; imagesc(z_med); colormap(gray); axis image;
title(['The denoised image by 3 \times 3 Median filter, RMSE = ', num2str(rmse_median)]);
figure; imagesc(zcL1(:,:,1)); colormap(gray); axis image;
title(['The denoised image by L1 second order classic kernel regression, RMSE = ', num2str(rmse_ckrL1)]);
figure; imagesc(zsL1(:,:,1)); colormap(gray); axis image;
title(['The denoised image by L1 second order steering kernel regression, RMSE = ', num2str(rmse_skrL1)]);


