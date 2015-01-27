function z = YCC2RGB(Y, Cr, Cb);

[N M num] = size(Y);
z = zeros(N, M, 3, num);

for i = 1 : num
    z(:,:,1,i) = Y(:,:,i) - 0.00092460 * Cr(:,:,i) + 1.40168676 * Cb(:,:,i);
    z(:,:,2,i) = Y(:,:,i) - 0.34369538 * Cr(:,:,i) - 0.71416904 * Cb(:,:,i);
    z(:,:,3,i) = Y(:,:,i) + 1.77216042 * Cr(:,:,i) + 0.00099022 * Cb(:,:,i);
end