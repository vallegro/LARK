function [z, zx1, zx2] = skr0_irregular(y, I, h, C, ksize)
% [SKR0_IRREGULAR]
% The zeroth order steering kernel regression function for irregularly
% sampled data.
%
% [USAGE]
% z = skr0_irregular(y, I, h, C, ksize)
%
% [RETURNS]
% z     : the estimated image
%
% [PARAMETERS]
% y     : the input image
% I     : The sampling position map (I(n,m) == 0 for no sample at (n,m))
% h     : the global smoothing parameter
% C     : the covariance matrices containing local orientation information
% ksize : the size of the kernel (ksize x ksize, and "ksize" must be
%         an odd number)
%
% [HISTORY]
% June 16, 2007 : created by Hiro
% Apr  14, 2008 : the transpose operator is fixed by Hiro

% Get the oritinal image size
[N, M] = size(y);

% Initialize the return parameters
z = zeros(N, M);
zx1 = zeros(N, M);
zx2 = zeros(N, M);

% Pixel sampling positions
radius = (ksize - 1) / 2;
[x2, x1] = meshgrid(-radius:radius, -radius:radius);

% The feture matrix
Xx = ones(ksize^2,1);

% pre-culculation for covariance matrices
C11 = zeros(N, M);
C12 = zeros(N, M);
C22 = zeros(N, M);
sq_detC = zeros(N, M);
for n = 1 : N
    for m = 1 : M
        if I(n,m) == 0
            continue;
        end
        C11(n,m) = C(1,1,n,m);
        C12(n,m) = C(1,2,n,m);
        C22(n,m) = C(2,2,n,m);
        sq_detC(n,m) = sqrt(det(C(:,:,n,m)));
    end
end

% Mirroring
y = EdgeMirror(y, [radius, radius]);
% I = EdgeMirror(I, [radius, radius]);
C11 = EdgeMirror(C11, [radius, radius]);
C12 = EdgeMirror(C12, [radius, radius]);
C22 = EdgeMirror(C22, [radius, radius]);
sq_detC = EdgeMirror(sq_detC, [radius, radius]);

% Estimate an image and its first gradients with pixel-by-pixel        
for n = 1 : N
    for m = 1 : M
            
        % Neighboring samples to be taken account into the estimation
        yp = y(n:n+ksize-1, m:m+ksize-1);
        
        % compute the weight matrix
        tt = x1 .* (C11(n:n+ksize-1, m:m+ksize-1) .* x1...
                + C12(n:n+ksize-1, m:m+ksize-1) .* x2)...
                + x2 .* (C12(n:n+ksize-1, m:m+ksize-1) .* x1...
                + C22(n:n+ksize-1, m:m+ksize-1) .* x2);
        WI = exp(-(0.5/h^2) * tt) .* sq_detC(n:n+ksize-1, m:m+ksize-1);
            
%         % Zero weights where there is no sample
%         WI = W .* I(n:n+ksize-1, m:m+ksize-1);
                
        % Equivalent kernel
        Xw = WI(:);
        A = inv(Xx' * Xw) * (Xw');
                
        % Estimate the pixel values at (nn,mm)
        z(n,m)   = A(1,:) * yp(:);
        
    end
end