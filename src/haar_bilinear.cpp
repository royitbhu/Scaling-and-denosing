// g++ haar_bilinear.cpp -I /usr/local/include/opencv -L /usr/local/lib -lml -lcv -lhighgui -lcvaux -lcxcore 


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cv.h>
#include <string.h>
#include <highgui.h>

//...................................... Scaling..................................................//

double valueat(double intensity1,double intensity2,double v1,double v2,double v)
{
    if(v1==v2)
        return intensity1;
    else
        return((((v2-v)/(v2-v1))*intensity2) +(((v-v1)/(v2-v1))*intensity1));
}

void bilinear(IplImage* img, IplImage *imgfinal, double horizontal_scale, double vertical_scale)
{
    int orig_height,orig_width,orig_step,orig_channels,new_step,new_channels,new_height,new_width,i,j,k;
    CvScalar v;
    orig_height=img->height;
    orig_width=img->width;
    orig_step=img->widthStep;
    new_step=orig_step*horizontal_scale;
    new_channels=orig_channels=img->nChannels;
    uchar *orig_data = (uchar *)img->imageData;	//Get the image data in 'data'
    new_height=orig_height*vertical_scale;			//new image height
  	new_width=orig_width*horizontal_scale;				//new image width
  	uchar *new_data =(uchar*) imgfinal->imageData;
  	for(i=0;i<new_height;i++) 
  		for(j=0;j<new_width;j++) 
  			for(k=0;k<new_channels;k++)
    				new_data[(int)(i*new_step+j*new_channels+k)]=0;
    for(i=0;i<new_height;i++)
    {
        for(j=0;j<new_width;j++)
        {
            int x1=floor((double)j/(double)horizontal_scale);
            int y1=floor((double)i/(double)vertical_scale);
            int x2=ceil((double)j/(double)horizontal_scale);
            int y2=ceil((double)i/(double)vertical_scale);
            double x = (double)j/(double)horizontal_scale;
            double y = (double)i/(double)vertical_scale;
            for(k=0;k<new_channels;k++)
            {
                
                double r1=valueat((double)orig_data[(int)(y1*orig_step +x1*orig_channels+k)],(double)orig_data[(int)(y1*orig_step +x2*orig_channels+k)],(double)x1,(double)x2,x);
                double r2=valueat((double)orig_data[(int)(y2*orig_step +x1*orig_channels+k)],(double)orig_data[(int)(y2*orig_step +x2*orig_channels+k)],(double)x1,(double)x2,x);
                v.val[k]=valueat(r1,r2,(double)y1,(double)y2,y);                
                //new_data[(int)(i*new_step+j*new_channels+k)]=orig_data[(int)(temp_height*orig_step + temp_width*orig_channels+k)];
            }
            cvSet2D(imgfinal,i,j,v);
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
    char input[20], output[20];
    scanf("%s%s%d%d", input, output, &o_width, &o_height);
    
    IplImage *in_image, *out_image;
    in_image = cvLoadImage(input,-1);
    
    if(!in_image)
    {
        printf("Error!\n");
        return 0;
    }
    
    out_image = cvCreateImage(cvSize(o_width,o_height),IPL_DEPTH_8U,in_image->nChannels);
    
    bilinear(in_image, out_image, o_width*1.0/in_image->width, o_height*1.0/in_image->height);
    Wavelet_Denoise(out_image);
    
    cvSaveImage(output,out_image);
}
