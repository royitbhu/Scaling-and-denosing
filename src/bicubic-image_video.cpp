/*
The program has been built using the 'opencv' API on linux platform.
Usage:	
	Compile :g++ program3.cpp -o program -I /usr/local/include/opencv -L /usr/local/lib -lml -lcv -lhighgui -lcvaux -lcxcore 
	Run 	:./a.out image-file-name
Result is saved in a file "result.jpg".
Original image is saved in file "original.jpg"

Approximation : The pixel areas on the boundaries don't have 16 corresponding pixels. For those pixels going out of range, their neighbouring pixels have been taken instead.
Thus, for a pixel(x,y)
If x is beyond width limit, x-1 will be used in its place.
If x is below zero, x+1 (zero) will be used in its place.
Similarly, y part is dealt with.
Pixel positions: 00->lower left;10->lower right; 01-> top left; 11-> top right


Reference : http://en.wikipedia.org/wiki/Bicubic_interpolation
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cv.h>
#include <highgui.h>

double a[4][4];
double calculatep(double x,double y)
{
	double p;
	p=a[0][0] + a[0][1]*y + a[0][2]*y*y + a[0][3]*y*y*y + a[1][0]*x + a[1][1]*x*y + a[1][2]*x*y*y + a[1][3]*x*y*y*y + a[2][0]*x*x + a[2][1]*x*x*y + a[2][2]*x*x*y*y + a[2][3]*x*x*y*y*y + a[3][0]*x*x*x + a[3][1]*x*x*x*y +a[3][2]*x*x*x*y*y +a[3][3]*x*x*x*y*y*y;
	return p;
}
void bicubic(IplImage* img, IplImage* imgfinal, double m, double n)
{
    CvScalar v;
	
	int height,width,step,channels;
	uchar *data;
	int i,j,k,f,g,h;
	height=img->height;
	width=img->width;
	step=img->widthStep;
	channels=img->nChannels;
	data      = (uchar *)img->imageData;	//Get the image data in 'data'
	int ht=height*n;			//new image height
  	int wt=width*m;				//new image width
  	uchar *datafinal =(uchar*) imgfinal->imageData;
  	for(i=0;i<height*n;i++) 
  		for(j=0;j<width*m;j++) 
  			for(k=0;k<channels;k++)
    				datafinal[(int)(i*step*m+j*channels+k)]=0;

  	double w[4],x[4],y[4],z[4],i1,j1;	//w->f00,f10,f01,f11; x->fx00,fx10,fx01,fx11; y->fy00,fy10,fy01,fy11; z->fxy0,fxy10,fxy01,fxy11
  	for(f=0;f<ht;f++)			//Traverse for all the pixels in new image
  	{
  		for(g=0;g<wt;g++)
  		{
  			i1=((double)f)/n;	
  			j1=((double)g)/m;
  			i=floor(i1);		//Map to the original image pixels
  			j=floor(j1);
  			i1=i1-i;
  			j1=j1-j;
  			
  			for(k=0;k<channels;k++)
  			{
  				
  				//Calculate f values
  				w[0]=data[i*step+j*channels+k];			
  				w[1]=data[i*step+(j+1)*channels+k];
  				w[2]=data[(i+1)*step+j*channels+k];
  				w[3]=data[(i+1)*step+(j+1)*channels+k];
  				
  				//Calulate fy values
  				if((i-1)>=0)//Check if value is going beyond the scope of image
  				{
  					y[0] = 0.5*(data[(i+1)*step+j*channels+k] - data[(i-1)*step+j*channels+k]);
  					y[1] = 0.5*(data[(i+1)*step+(j+1)*channels+k] - data[(i-1)*step+(j+1)*channels+k]);
  				}
  				else// Selecting the neighbouring pixels in the place
  				{
  					y[0]=0.5*(data[(i+1)*step+j*channels+k]- data[(i)*step+j*channels+k]);
  					y[1]=0.5*(data[(i+1)*step+(j+1)*channels+k]- data[(i)*step+(j+1)*channels+k]);
  				}
  				if((i+2)<=(height-1))
  				{
  					y[2] = 0.5*(data[(i+2)*step+j*channels+k] - data[(i)*step+j*channels+k]);
  					y[3] = 0.5*(data[(i+2)*step+(j+1)*channels+k] - data[(i)*step+(j+1)*channels+k]);
  				}
  				else
  				{
  					y[2] = 0.5*(data[(i+1)*step+j*channels+k] - data[(i)*step+j*channels+k]);
  					y[3] = 0.5*(data[(i+1)*step+(j+1)*channels+k] - data[(i)*step+(j+1)*channels+k]);
  				}
  				
  				//Calculate fx values
  				if((j-1)>=0)
  				{
  					x[0]=0.5*(data[i*step+(j+1)*channels+k] - data[i*step+(j-1)*channels+k]);
  					x[2]=0.5*(data[(i+1)*step+(j+1)*channels+k] - data[(i+1)*step+(j-1)*channels+k]);
  				}
  				else
  				{
  					x[0]=0.5*(data[i*step+(j+1)*channels+k] - data[i*step+(j)*channels+k]);
  					x[2]=0.5*(data[(i+1)*step+(j+1)*channels+k] - data[(i+1)*step+(j)*channels+k]);
  				}
  				if((j+2)<=(width-1))
  				{
  					x[1]=0.5*(data[i*step+(j+2)*channels+k] - data[(i)*step+j*channels+k]);
  					x[3]=0.5*(data[(i+1)*step+(j+2)*channels+k] - data[(i+1)*step+(j)*channels+k]);
  				}
  				else
  				{
  					x[1]=0.5*(data[i*step+(j+1)*channels+k] - data[(i)*step+j*channels+k]);
  					x[3]=0.5*(data[(i+1)*step+(j+1)*channels+k] - data[(i+1)*step+(j)*channels+k]);
  				}
  				
  				//Calculate fxy values
  				if((i-1)>=0 && (j-1)>=0)
  					z[0]=0.25*(data[(i+1)*step+(j+1)*channels+k] - data[(i-1)*step+(j+1)*channels+k] - data[(i+1)*step+(j-1)*channels+k] + data[(i-1)*step+(j-1)*channels+k]);
  				else if((i-1)>=0)
  				{
  					z[0]=0.25*(data[(i+1)*step+(j+1)*channels+k] - data[(i-1)*step+(j+1)*channels+k] - data[(i+1)*step+(j)*channels+k] + data[(i-1)*step+(j)*channels+k]);
  				}
  				else if((j-1)>=0)
  				{
  					z[0]=0.25*(data[(i+1)*step+(j+1)*channels+k] - data[(i)*step+(j+1)*channels+k] - data[(i+1)*step+(j-1)*channels+k] + data[(i)*step+(j-1)*channels+k]);
  				}
  				else
  				{
  					z[0]=0.25*(data[(i+1)*step+(j+1)*channels+k] - data[(i)*step+(j+1)*channels+k] - data[(i+1)*step+(j)*channels+k] + data[(i)*step+(j)*channels+k]);
  				}
  				if((j+2)<width && (i-1)>=0)
  				{
  					z[1]=0.25*(data[(i+1)*step+(j+2)*channels+k]-data[(i-1)*step+(j+2)*channels+k]-data[(i+1)*step+(j)*channels+k]+data[(i-1)*step+(j)*channels+k]);
  				}
  				else if((j+2)<width)
  				{
  					z[1]=0.25*(data[(i+1)*step+(j+2)*channels+k]-data[(i)*step+(j+2)*channels+k]-data[(i+1)*step+(j)*channels+k]+data[(i)*step+(j)*channels+k]);
  				}
  				else if((i-1)>=0)
  				{
  					z[1]=0.25*(data[(i+1)*step+(j+1)*channels+k]-data[(i-1)*step+(j+1)*channels+k]-data[(i+1)*step+(j)*channels+k]+data[(i-1)*step+(j)*channels+k]);
  				}
  				else
  				{
  					z[1]=0.25*(data[(i+1)*step+(j+1)*channels+k]-data[(i)*step+(j+1)*channels+k]-data[(i+1)*step+(j)*channels+k]+data[(i)*step+(j)*channels+k]);
  				}
  				if((i+2)<height && (j-1)>=0)
  				{
  					z[2]=0.25*(data[(i+2)*step+(j+1)*channels+k] - data[(i)*step+(j+1)*channels+k] - data[(i+2)*step+(j-1)*channels+k] + data[(i)*step+(j-1)*channels+k]);
  				}
  				else if((i+2)<height)
  				{
  					z[2]=0.25*(data[(i+2)*step+(j+1)*channels+k] - data[(i)*step+(j+1)*channels+k] - data[(i+2)*step+(j)*channels+k] + data[(i)*step+(j)*channels+k]);
  				}
  				else if((j-1)>=0)
  				{
  					z[2]=0.25*(data[(i+1)*step+(j+1)*channels+k] - data[(i)*step+(j+1)*channels+k] - data[(i+1)*step+(j-1)*channels+k] + data[(i)*step+(j-1)*channels+k]);
  				}
  				else
  				{
  					z[2]=0.25*(data[(i+1)*step+(j+1)*channels+k] - data[(i)*step+(j+1)*channels+k] - data[(i+1)*step+(j)*channels+k] + data[(i)*step+(j)*channels+k]);
  				}
  				if((j+2)<width && (i+2)<height)
  				{
  					z[3]=0.25*(data[(i+2)*step+(j+2)*channels+k] - data[(i)*step+(j+2)*channels+k] - data[(i+2)*step+(j)*channels+k] + data[i*step+j*channels+k]);
  				}
  				else if((j+2)<width)
  				{
  					z[3]=0.25*(data[(i+1)*step+(j+2)*channels+k] - data[(i)*step+(j+2)*channels+k] - data[(i+1)*step+(j)*channels+k] + data[i*step+j*channels+k]);
  				}
  				else if((i+2)<height)
  				{
  					z[3]=0.25*(data[(i+2)*step+(j+1)*channels+k] - data[(i)*step+(j+1)*channels+k] - data[(i+2)*step+(j)*channels+k] + data[i*step+j*channels+k]);
  				}
  				else
  				{
  					z[3]=0.25*(data[(i+1)*step+(j+1)*channels+k] - data[(i)*step+(j+1)*channels+k] - data[(i+1)*step+(j)*channels+k] + data[i*step+j*channels+k]);
  				}
  				
  				
  				//Find the sixteen coefficients required for the equation to find channel values
  				a[0][0]=w[0];
  				a[1][0]=y[0];
  				a[2][0]=0-(3*w[0])+(3*w[2])-(2*y[0])-y[2];
  				a[3][0]=(2*w[0])-(2*w[2])+y[0]+y[2];
  				a[0][1]=x[0];
  				a[1][1]=z[0];
  				a[2][1]=0-(3*x[0])+(3*x[2])-(2*z[0])-z[2];
  				a[3][1]=(2*x[0])-(2*x[2])+z[0]+z[2];
  				a[0][2]=0-(3*w[0])+(3*w[1])-(2*x[0])-x[1];
  				a[1][2]=0-(3*y[0])+(3*y[1])-(2*z[0])-z[1];
  				a[2][2]=(9*(w[0]-w[1]-w[2]+w[3]))+(6*(x[0]-x[2]+y[0]-y[1]))+(3*(x[1]-x[3]+y[2]-y[3]))+(4*z[0])+(2*z[1])+(2*z[2])+z[3];
  				a[3][2]=(6*(w[1]-w[0]+w[2]-w[3]))+(3*(y[1]-y[0]+y[3]-y[2]))+(2*(x[3]-x[1]-z[0]-z[2]))+(4*x[2])-4*x[0]-z[1]-z[3];
  				a[0][3]=(2*w[0])-(2*w[1])+x[0]+x[1];
  				a[1][3]=(2*y[0])-(2*y[1])+z[0]+z[1];
  				a[2][3]=(6*(w[1]-w[0]+w[2]-w[3]))+(3*(x[2]+x[3]-x[0]-x[1]))+(2*(y[3]-y[2]-z[0]-z[1])) -(4*y[0])+(4*y[1])-z[2]-z[3];
  				a[3][3]=(4*(w[0]-w[1]-w[2]+w[3]))+(2*(x[0]+x[1]-x[2]-x[3]+y[0]-y[1]+y[2]-y[3]))+z[0]+z[1]+z[2]+z[3];
  				
  				//Calculate the values for the k channels in the image
  				v.val[k]=(int)(calculatep(j1,i1));
  				
  				//datafinal[(int)(f*step*m + g*channels+k)]=(int)(calculatep(j1,i1));
  			}
  			//Assign the values to the image
  			cvSet2D(imgfinal,f,g,v);
  		}
  	}	
}

int main()
{
    char input[20], output[20];
    int owidth, oheight;
    double h_scale, v_scale;
    scanf("%s",input);
    scanf("%s%d%d",output,&owidth,&oheight);
    
    IplImage *in_image, *out_image;
    
    CvCapture *capture = cvCreateFileCapture(input);
    CvVideoWriter *writer = 0;
    double fps  = cvGetCaptureProperty(capture, CV_CAP_PROP_FPS);//30
    int frameW  = owidth;
    int frameH  = oheight;
    h_scale = owidth*1.0/cvGetCaptureProperty(capture, CV_CAP_PROP_FRAME_WIDTH);
    v_scale = oheight*1.0/cvGetCaptureProperty(capture, CV_CAP_PROP_FRAME_HEIGHT);
    
    writer=cvCreateVideoWriter(output,cvGetCaptureProperty(capture, CV_CAP_PROP_FOURCC),//CV_FOURCC('P','I','M','1'),
                           fps,cvSize(frameW,frameH));
                           int i=0;
    while(cvGrabFrame(capture))
    {
        in_image=cvRetrieveFrame(capture);  
        out_image = cvCreateImage(cvSize(owidth,oheight),IPL_DEPTH_8U,in_image->nChannels);
        bicubic(in_image, out_image, h_scale, v_scale);
        cvWriteFrame(writer,out_image);  
        printf("%d\n",i++);    
    }
 
    cvReleaseVideoWriter(&writer);
    cvReleaseCapture(&capture);
}
