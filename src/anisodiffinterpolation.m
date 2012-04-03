%%%%%%%%%%%Anisotropic diffusion based interpolation: Guichard et al %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%Author: Rajeev Srivastava, ITBHU%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 							
I=imread('input.jpg');
p=4;%%%Zoom out factor
p1=p./(p*p);%%%Zoom-in factor
%%%%Downsample the image by factor qxq : By default downsample and upsample  
J1=imresize(I,p1,'nearest');
%%%%%%%%%%%Upsample the image to its original size by applysing Bilinear or Cubic interepolation
 J = imresize(I,p,'bicubic');% Magnification factor 2x2, 0.5x0.5
  figure (1),subplot (3,2,1), imshow(uint8(I),[]),title ('Original Image');
  figure (1),subplot (3,2,2), imshow(uint8(J1),[]),title ('Downsampled Image');
 figure (1),subplot (3,2,3), imshow(uint8(J),[]),title ('Bilinear Interpolation');
 %  figure (2),imshow(uint8(J),[]);title ('Bilinear Interpolation /zoom');
  %horg=imhist(uint8(I));
  %hdownsampled=imhist(uint8(J1));
  %hzoom=imhist(uint8(J));
 %figure(2),subplot (3,2,1),plot(horg),title('Histogram of original image');
 %figure(2),subplot (3,2,2),plot(hdownsampled),title('Histogram of downsampled image');
 %figure(2),subplot (3,2,3),plot(hzoom),title('bilinear interpolation');
 %figure(4),plot(horg),title('Histogram of original image');
 
 % Step 2 : Computation of Projection operator :Perform image averaging
%  starting from initial condition.This will be updated in each iteration of pde based filtering
h = fspecial('average',[3 3]);
Iproj= imfilter(J,h);


 %%%%%%%%%%%%Step 3: Using nonlinear complex diffusion process for magnification, noise smoothing and edge forming %%%%%%%%%%%
Iproj=double(Iproj);
J=double(J);
[x y z]=size(J);
dt=0.25; %dt<0.25 for stability of numerical scheme 
im=double(J);
kappa=60;
option=2;% option 1 is associated with less PSNR
t=1;% niter initialized 
% T=4; %Stopping Criteria T=pxp (Magnification size) e.g. T=4 for 2x2
% magnification, T=niter
T=p*p;
% No of iterations= 1 to 4 in increments of 0.20 i.e. 0.2xniter=4 i.e niter=20
%%%%%%%%%%%%%%%%%%%%Anisotropic diff based interpolation code%%%%%%%%
%if ndims(im)==3
  %error('Anisodiff only operates on 2D gray-scale images');
%end
% figure(1);imshow(im); title('Original noisy Image');
im = double(im);
[rows,cols,channels] = size(im);
diff = im;
for t = 1:T
%  fprintf('\rIteration %d',i);

  % Construct diffl which is the same as diff but
  % has an extra padding of zeros around it.
  diffl = zeros(rows+2, cols+2,channels+2);
  diffl(2:rows+1, 2:cols+1,1:channels) = diff;

  % North, South, East and West differences
  deltaN = diffl(1:rows,2:cols+1,1:channels)   - diff;
  deltaS = diffl(3:rows+2,2:cols+1,1:channels) - diff;
  deltaE = diffl(2:rows+1,3:cols+2,1:channels) - diff;
  deltaW = diffl(2:rows+1,1:cols,1:channels)   - diff;

  % Conduction

  if option == 1
    cN = exp(-(deltaN/kappa).^2);
    cS = exp(-(deltaS/kappa).^2);
    cE = exp(-(deltaE/kappa).^2);
    cW = exp(-(deltaW/kappa).^2);
  elseif option == 2
    cN = 1./(1 + (deltaN/kappa).^2);
    cS = 1./(1 + (deltaS/kappa).^2);
    cE = 1./(1 + (deltaE/kappa).^2);
    cW = 1./(1 + (deltaW/kappa).^2);
  end
 diff = diff + dt*((cN.*deltaN + cS.*deltaS + cE.*deltaE + cW.*deltaW)-Iproj+J);
% diff = diff + dt*(cN.*deltaN + cS.*deltaS + cE.*deltaE + cW.*deltaW);
Ianiso=diff;
end
figure (5), imshow(uint8(diff),[]),title ('Anisotropic diff Interpolation:Zoom out 4x4');
imwrite(uint8(diff),'output.jpg','jpg');