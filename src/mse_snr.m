function [psnr, mse] = mse_snr(X,Xc)

D = abs(X-Xc).^2;
mse  = sum(D(:))/numel(X)

psnr = 10*log10(255*255/mse)