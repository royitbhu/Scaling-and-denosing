#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <cv.h>
#include <highgui.h>

#define pi 3.14159265

int a=3;

double Lanczos(double x)
{
    if(x<0.000001 && x>-0.000001)
        return 1;
    else if(-a < x && a > x)
        return (a*sin(pi*x)*sin(pi*x/a))/(pi*pi*x*x);
    else
        return 0;
}

void Lanczose_Scale(IplImage *in, IplImage *out, int owidth, int oheight)
{
    int channels = in->nChannels;
    int iwidth = in->width;
    int iheight = in->height;
    
    int i,j,k,x,y;
    for(x=0;x<owidth;x++)
        for(y=0;y<oheight;y++)
        {
            double x0 = floor(x*1.0*iwidth/owidth);
            double y0 = floor(y*1.0*iheight/oheight);
            CvScalar temp, temp2;
            for(k=0;k<channels;k++)
                temp2.val[k] = 0;
            for(i=x0-a+1;i<=x0+a;i++)
                for(j=y0-a+1;j<=y0+a;j++)
                    if(i>=0 && i<iwidth && j>=0 && j<iheight)
                    {
                        temp = cvGet2D(in,j,i);
                        for(k=0;k<channels;k++)
                        {
                            temp2.val[k] += temp.val[k]*Lanczos(x0-i)*Lanczos(y0-j);
                        }
                    }
            cvSet2D(out,y,x,temp2);
        }
}

int main()
{
    char input[20], output[20];
    int owidth, oheight;
    scanf("%s",input);
    scanf("%s%d%d",output,&owidth,&oheight);
    
    IplImage *in_image, *out_image;
    in_image = cvLoadImage(input,-1);
    
    if(!in_image)
    {
        printf("Error!\n");
        return 0;
    }
    
    out_image = cvCreateImage(cvSize(owidth,oheight),IPL_DEPTH_8U,in_image->nChannels);
    
    Lanczose_Scale(in_image, out_image, owidth, oheight);
    
    cvSaveImage(output,out_image);
}
