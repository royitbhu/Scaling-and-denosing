// g++ haar_bilinear.cpp -I /usr/local/include/opencv -L /usr/local/lib -lml -lcv -lhighgui -lcvaux -lcxcore 


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cv.h>
#include <string.h>
#include <highgui.h>

//...................................... Scaling..................................................//

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

//..............................................................................//

void Haar1DRow(double *image, int row, int col, int width, int channels)
{
    double *data = new double [col*channels];
    int i,j,k;
    for(i=0;i<row;i++)
    {
        for(j=0;j<col/2;j++)
        {
            for(k=0;k<channels;k++)
            {
                double a = image[(i*width + 2*j)*channels + k];
                double b = image[(i*width + (2*j+1))*channels + k];
                data[j*channels + k] = (a+b)/sqrt(2);
                data[(col/2+j)*channels + k] = (a-b)/sqrt(2);
            }
        }
        for(j=0;j<col/2;j++)
        {   for(k=0;k<channels;k++)
            {
                image[(i*width+j)*channels + k] = data[j*channels + k];
                image[(i*width+col/2+j)*channels + k] = data[(col/2+j)*channels + k];
            }
        }
    }
    
}

void UnHaar1DRow(double *image, int row, int col, int width, int channels)
{
    double *data = new double [col*channels];
    int i,j,k;
    for(i=0;i<row;i++)
    {
        for(j=0;j<col/2;j++)
        {
            for(k=0;k<channels;k++)
            {
                double a = image[(i*width+j)*channels + k];
                double b = image[(i*width+col/2+j)*channels + k];
                data[j*channels + k] = (a+b)/sqrt(2);
                data[(col/2+j)*channels + k] = (a-b)/sqrt(2);
            }
        }
        for(j=0;j<col/2;j++)
        {   for(k=0;k<channels;k++)
            {
                image[(i*width + 2*j)*channels + k] = data[j*channels + k];
                image[(i*width + (2*j+1))*channels + k] = data[(col/2+j)*channels + k];
            }
        }
    }
    
}

void Haar1DCol(double *image, int row, int col, int width, int channels)
{
    double *data = new double [row*channels];
    int i,j,k;
    for(i=0;i<col;i++)
    {
        for(j=0;j<row/2;j++)
        {
            for(k=0;k<channels;k++)
            {
                double a = image[(2*j*width + i)*channels + k];
                double b = image[((2*j+1)*width + i)*channels + k];
                data[j*channels + k] = (a+b)/sqrt(2);
                data[(row/2+j)*channels + k] = (a-b)/sqrt(2);
            }
        }
        for(j=0;j<row/2;j++)
        {   for(k=0;k<channels;k++)
            {
                image[(j*width+i)*channels + k] = data[j*channels + k];
                image[((row/2+j)*width+i)*channels + k] = data[(row/2+j)*channels + k];
            }
        }
    }
}

void UnHaar1DCol(double *image, int row, int col, int width, int channels)
{
    double *data = new double [row*channels];
    int i,j,k;
    for(i=0;i<col;i++)
    {
        for(j=0;j<row/2;j++)
        {
            for(k=0;k<channels;k++)
            {
                double a = image[(j*width+i)*channels + k];
                double b = image[((row/2+j)*width+i)*channels + k];
                data[j*channels + k] = (a+b)/sqrt(2);
                data[(row/2+j)*channels + k] = (a-b)/sqrt(2);
            }
        }
        for(j=0;j<row/2;j++)
        {   for(k=0;k<channels;k++)
            {
                image[(2*j*width + i)*channels + k] = data[j*channels + k];
                image[((2*j+1)*width + i)*channels + k] = data[(row/2+j)*channels + k];
            }
        }
    }
}

void Haar2D(double *image, int width, int height, int channels)
{
    int level=3, crow = height, ccol = width;
    
    while(level--)
    {
        Haar1DRow(image, crow, ccol, width, channels);
        Haar1DCol(image, crow, ccol, width, channels);
        ccol /= 2;
        crow /= 2;
    }
}

void Haar_Reconstruct(double *image, int width, int height, int channels)
{
    int level=3, crow = height/4, ccol = width/4;
    
    while(level--)
    {
        UnHaar1DCol(image, crow, ccol, width, channels);
        UnHaar1DRow(image, crow, ccol, width, channels);
        ccol *= 2;
        crow *= 2;
    }
}

void Denoise(double *data, int width, int height, int nChannels)
{
    int i,j,k, l;
    //double sigma[3],temp, thr[3];
    double sigx[2000][3], sigy[2000][3], temp, th;
   
    for(l=1;l<=3;l++)
    {
        /*sigma[0] = sigma[1] = sigma[2] = 0.0;
        for(i=0;i<height/pow(2,l);i++)
            for(j=0;j<width/pow(2,l);j++)
                for(k=0;k<nChannels;k++)
                {
                    sigma[k] += abs(data[(i*width + j)*nChannels + k]);
                }
        sigma[0] /= 0.6745*width*height/pow(2,2*l);
        sigma[1] /= 0.6745*width*height/pow(2,2*l);
        sigma[2] /= 0.6745*width*height/pow(2,2*l);
        
        thr[0] = sigma[0]*sqrt(log(width*height/pow(2,2*l))/(width*height/pow(2,2*l)));
        thr[1] = sigma[1]*sqrt(log(width*height/pow(2,2*l))/(width*height/pow(2,2*l)));
        thr[2] = sigma[2]*sqrt(log(width*height/pow(2,2*l))/(width*height/pow(2,2*l)));
        */
        
        for(i=0;i<2000;i++)
            sigx[i][0] = sigx[i][1] = sigx[i][2] = sigy[i][0] = sigy[i][1] = sigy[i][2] = 0.0;
        
        for(i=0;i<height/pow(2,l-1);i++)
            for(j=0;j<width/pow(2,l-1);j++)
                for(k=0;k<3;k++)
                {
                    sigx[j][k] += abs(data[(i*width + j)*nChannels + k]);
                    sigy[i][k] += abs(data[(i*width + j)*nChannels + k]);
                }
        
        for(i=0;i<height/pow(2,l-1);i++)
            for(j=0;j<width/pow(2,l-1);j++)
                for(k=0;k<nChannels;k++)
                {
                    temp = data[(i*width + j)*nChannels + k];
                    
                    double sig = (sigx[j][k]+sigy[i][k])/2;
                    sig /= 0.6745*(width+height)/pow(2,l-1);
                    th = sig*sqrt(log((width+height)/pow(2,l-1))/((width+height)/pow(2,l-1)));
                    
                    if(abs(temp) < th)
                        data[(i*width + j)*nChannels + k] = 0.0;
                    else if(temp > th)
                        data[(i*width + j)*nChannels + k] -= th;
                    else
                        data[(i*width + j)*nChannels + k] += th;
                }
    }
}

void Wavelet_Denoise(IplImage *image)
{
    int i,j,k, width, height, nChannels;
    CvScalar temp;
    width = image->width;
    height = image->height;
    nChannels = image->nChannels;
    double *data = new double [width * height * nChannels];
    
    for(i=0;i<height;i++)
        for(j=0;j<width;j++)
        {
            temp = cvGet2D(image, i, j);
            for(k=0;k<nChannels;k++)
                data[(i*width + j)*nChannels + k] = temp.val[k];
        }
    
    Haar2D(data, width, height, nChannels);
    
    Denoise(data, width, height, nChannels);
    
    Haar_Reconstruct(data, width, height, nChannels);
    
    for(i=0;i<height;i++)
        for(j=0;j<width;j++)
        {
            for(k=0;k<nChannels;k++)
                temp.val[k] = data[(i*width + j)*nChannels + k];
            cvSet2D(image, i, j, temp);
        }
    
    delete [] data;
}

int main()
{
    int o_height, o_width;
    char input[40], output[40];
    scanf("%s", input);
    scanf("%s", output);
    scanf("%d%d", &o_width, &o_height);
    
    IplImage *in_image, *out_image;
    in_image = cvLoadImage(input);
    
    if(!in_image)
    {
        printf("Error!\n");
        return 0;
    }
    
    out_image = cvCreateImage(cvSize(o_width,o_height),IPL_DEPTH_8U,in_image->nChannels);
    
    bicubic(in_image, out_image, o_width*1.0/in_image->width, o_height*1.0/in_image->height);
    Wavelet_Denoise(out_image);
    
    cvSaveImage(output,out_image);
}
