function z = ckr2L1_regular(y, z_init, h, r, ksize, IT, step)
% [CKR2L1_REGULAR]
% The second order classic kernel regression function with L1-norm
% for regularly sampled data.
%
% [USAGE]
% z = ckr2L1_regular(y, z_init h, r, ksize, IT)
%
% [RETURNS]
% z     : the estimated image and all the derivatives (first and second)
%
% [PARAMETERS]
% y      : the input image
% z_init : the initial state for z 
% h      : the global smoothing parameter
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

% Mirroring
y = EdgeMirror(y, [radius, radius]);

% Estimate an image and its first gradients with pixel-by-pixel
for i = 1 : r
    for j = 1 : r
        xx1 = downsample2(x1(r-i+1:end, r-j+1:end), r);
        xx2 = downsample2(x2(r-i+1:end, r-j+1:end), r);
        
        % The feture matrix
        Xx = [ones(ksize^2,1), xx1(:), xx2(:), xx1(:).^2, xx1(:).*xx2(:), xx2(:).^2];
                
        % The weight matrix
        W = exp(-(0.5/h^2) * (xx1.^2 + xx2.^2));
        
        % Xx^T * W
        XwT = [Xx(:,1).*W(:), Xx(:,2).*W(:), Xx(:,3).*W(:),...
               Xx(:,4).*W(:), Xx(:,5).*W(:), Xx(:,6).*W(:)].';
        
        for n = 1 : N
            nn = (n - 1) * r + i;
            for m = 1 : M
                mm = (m - 1) * r + j;
                
                % Neighboring samples to be taken account into the estimation
                yp = y(n:n+ksize-1, m:m+ksize-1);
                yp = yp(:);
                
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



