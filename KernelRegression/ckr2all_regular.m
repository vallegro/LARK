function z = ckr2all_regular(y, h, r, ksize)
% [CKR2ALL_REGULAR]
% The second order classic kernel regression function for regularly sampled
% data.
%
% [USAGE]
% z = ckr2all_regular(y, h, r, ksize)
%
% [RETURNS]
% z     : the estimated image and all the derivatives (first and second)
%
% [PARAMETERS]
% y     : the input image
% h     : the global smoothing parameter
% r     : the upscaling factor ("r" must be an integer number)
% ksize : the size of the kernel (ksize x ksize, and "ksize" must be
%         an odd number)
%
% [HISTORY]
% July, 01, 2007 : created by Hiro
% Apr  14, 2008 : the transpose operator is fixed by Hiro

% Get the oritinal image size
[N, M] = size(y);

% Initialize the return parameters
z = zeros(N*r, M*r, 6);

% Create the equivalent kernels
radius = (ksize - 1) / 2;
[x2, x1] = meshgrid(-radius-(r-1)/r : 1/r : radius, -radius-(r-1)/r : 1/r : radius);
A = zeros(6, ksize^2, r, r);
for i = 1 : r
    for j = 1 : r 
        xx1 = downsample2(x1(r-i+1:end, r-j+1:end), r);
        xx2 = downsample2(x2(r-i+1:end, r-j+1:end), r);
        % The feture matrix
        Xx = [ones(ksize^2,1), xx1(:), xx2(:), xx1(:).^2, xx1(:).*xx2(:), xx2(:).^2];
        % The weight matrix (Gaussian kernel function)
        tt = xx1.^2 + xx2.^2;
        W = exp(-(0.5/h^2) * tt);
        % Equivalent kernel
        Xw = [Xx(:,1).*W(:), Xx(:,2).*W(:), Xx(:,3).*W(:),...
                Xx(:,4).*W(:), Xx(:,5).*W(:), Xx(:,6).*W(:)];
        A(:,:,i,j) = inv(Xx' * Xw) * (Xw');
    end
end

% Mirroring the input image
y = EdgeMirror(y, [radius, radius]);

% Estimate an image and its first gradients with pixel-by-pixel
for n = 1 : N
    for m = 1 : M
        
        % Neighboring samples to be taken account into the estimation
        yp = y(n:n+ksize-1, m:m+ksize-1);
        
        for i = 1 : r
            nn = (n - 1) * r + i;
            for j = 1 : r
                mm = (m - 1) * r + j;
                % Estimate the pixel values at (nn,mm)
                b = A(:,:,i,j) * yp(:);
                z(nn,mm,1) = b(1); % z(x)
                z(nn,mm,2) = b(2); % z_x1(x)
                z(nn,mm,3) = b(3); % z_x2(x)
                z(nn,mm,4) = b(4); % z_x1x1(x)
                z(nn,mm,5) = b(5); % z_x1x2(x)
                z(nn,mm,6) = b(6); % z_x2x2(x)
            end
        end
        
    end
end