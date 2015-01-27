function z = skr2L1_regular(y, z_init, h, C, r, ksize, IT, step)
% [SKR2L1_REGULAR]
% The second order slassic kernel regression function with L1-norm
% for regularly sampled data.
%
% [USAGE]
% z = skr2L1_regular(y, z_init h, C, r, ksize, IT)
%
% [RETURNS]
% z     : the estimated image and all the derivatives (first and second)
%
% [PARAMETERS]
% y      : the input image
% z_init : the initial state for z 
% h      : the global smoothing parameter
% C      : the covariance matrices containing local orientation information
% r      : the upscaling factor ("r" must be an integer number)
% ksize  : the size of the kernel (ksize x ksize, and "ksize" must be
%          an odd number)
% IT     : the number of iterations for the steepest descnet method
% step   : the step size of the steepest descent update
%
% [HISTORY]
% July 01, 2007 : created by Hiro

% Get the oritinal image size
[N, M] = size(y);

% Initialize the return parameters
z = zeros(N*r, M*r, 6);

% Pixel sampling positions
radius = (ksize - 1) / 2;
[x2, x1] = meshgrid(-radius-(r-1)/r : 1/r : radius, -radius-(r-1)/r : 1/r : radius);

% pre-culculation for covariance matrices
C11 = zeros(N, M);
C12 = zeros(N, M);
C22 = zeros(N, M);
sq_detC = zeros(N, M);
for n = 1 : N
    for m = 1 : M
        C11(n,m) = C(1,1,n,m);
        C12(n,m) = C(1,2,n,m);
        C22(n,m) = C(2,2,n,m);
        sq_detC(n,m) = sqrt(det(C(:,:,n,m)));
    end
end

% Mirroring
y = EdgeMirror(y, [radius, radius]);
C11 = EdgeMirror(C11, [radius, radius]);
C12 = EdgeMirror(C12, [radius, radius]);
C22 = EdgeMirror(C22, [radius, radius]);
sq_detC = EdgeMirror(sq_detC, [radius, radius]);

% Estimate an image and its first gradients with pixel-by-pixel
for i = 1 : r
    for j = 1 : r
        xx1 = downsample2(x1(r-i+1:end, r-j+1:end), r);
        xx2 = downsample2(x2(r-i+1:end, r-j+1:end), r);
        
        % The feture matrix
        Xx = [ones(ksize^2,1), xx1(:), xx2(:), xx1(:).^2, xx1(:).*xx2(:), xx2(:).^2];
        
        for n = 1 : N
            nn = (n - 1) * r + i;
            for m = 1 : M
                mm = (m - 1) * r + j;
                
                % Neighboring samples to be taken account into the estimation
                yp = y(n:n+ksize-1, m:m+ksize-1);
                yp = yp(:);
                
                % compute the weight matrix
                tt = xx1 .* (C11(n:n+ksize-1, m:m+ksize-1) .* xx1...
                        + C12(n:n+ksize-1, m:m+ksize-1) .* xx2)...
                        + xx2 .* (C12(n:n+ksize-1, m:m+ksize-1) .* xx1...
                        + C22(n:n+ksize-1, m:m+ksize-1) .* xx2);
                W = exp(-(0.5/h^2) * tt) .* sq_detC(n:n+ksize-1, m:m+ksize-1);
                
                % Xx^T * W
                XwT = [Xx(:,1).*W(:), Xx(:,2).*W(:), Xx(:,3).*W(:),...
                        Xx(:,4).*W(:), Xx(:,5).*W(:), Xx(:,6).*W(:)].';
                
                % steepest descent iterations
                % initialization
                b = z_init(nn,mm,:);
                b = b(:);
                for it = 1 : IT
                    b = b + step .* (XwT * sign(yp - Xx * b));
                end
                
                % Estimate the pixel values at (nn,mm)
                z(nn,mm,1) = b(1);
                z(nn,mm,2) = b(2);
                z(nn,mm,3) = b(3);
                z(nn,mm,4) = b(4); 
                z(nn,mm,5) = b(5);
                z(nn,mm,6) = b(6);
            end
        end
    end
end