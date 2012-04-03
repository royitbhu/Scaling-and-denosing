%%%%%%%%%%%%%%%%%Performance analysis of zooming algorithms%%%%%%%%%%%%
function [MSE,rmse,SNR,PSNR,CP]= PerformanceAnalysisZoomMethods(I1,I2)
%%%%%%%%%%%%%%%%Performance Analysis of Speckle Reduction Methods%%%%%%%%%%%%
%-----------   PSNR analysis   ----------------
[x y z]=size(I1);
Q = 255;
MSE= sum(sum((I2-I1) .^ 2)) / (x * y);
rmse=sqrt(MSE);

% SNR ( sum of energies of original image/sum of energies of
% original-reconstructed image)
SNR=10*log10(sum(sum((I1).^2)))./(sum(sum((I2-I1) .^ 2)));
PSNR = 10*log10(Q*Q./MSE);
 
 %%%%%%%%%%%%%%%%%%%%Coorelation Parameter   -CP : metric for edge preservation%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% I1 is the original image
%%%%%% I2 is the reconstructed/denoised image
%%%%%%Ih1=high pass filtered version of I1 
%%%%%%Ih2 =high pass filtered version of I2
%%%%%% High pass filtered version of the images are obtained via a 3x3 pixel std  approximation of the Laplacian operator
h = fspecial('laplacian');
Ih1 = imfilter(I1,h,'replicate');
Ih2 = imfilter(I2,h,'replicate');
mI1=mean2(I1);
mI2=mean2(I2);
X1=sum(sum((Ih1-mI1).* (Ih2-mI2)));
X2=sum(sum((Ih1-mI1).*(Ih1-mI1)));
X3=sum(sum((Ih2-mI2).*(Ih2-mI2)));
CP=X1./sqrt(X2.*X3);% Correlation parameter ,it should be close to unity for an optimal effect of edge preservation
out(1) = PSNR(1,1,1);
out(1) = out(1)+PSNR(1,1,2);
out(1) = out(1)+PSNR(1,1,3);
out(2) = MSE(1,1,1);
out(2) = out(2)+MSE(1,1,2);
out(2) = out(2)+MSE(1,1,3);
out(3) = CP(1,1,1);
out(3) = out(3)+CP(1,1,2);
out(3) = out(3)+CP(1,1,3);
out./3
