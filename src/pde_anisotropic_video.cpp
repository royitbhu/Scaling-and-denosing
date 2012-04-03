/*
Usage:	
	Compile :g++ pde_anisotropic.cpp -I /usr/local/include/opencv -L /usr/local/lib -lml -lcv -lhighgui -lcvaux -lcxcore 
	Run 	:./a.out image-file-name
Result is saved in a file "result.jpg".

*/

#include<iostream>
#include<stdio.h>
#include <math.h>
#include <cv.h>
#include <highgui.h>
using namespace std;

#define GARBAGE 99999

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

void getrgb(IplImage *I,double *Ir,double *Ig,double *Ib)
{
    int w=I->width;
    int h=I->height;
    int i,j;
    CvScalar temp;
    for(i=0;i<w;i++)
    {
        for(j=0;j<h;j++)
        {
            temp=cvGet2D(I,j,i);
            Ir[j*w+i]=temp.val[0];
            Ig[j*w+i]=temp.val[1];
            Ib[j*w+i]=temp.val[2];
        }
    }
}
void generateImgFrmRGB(IplImage *img,double *Ir,double *Ig, double *Ib,int w,int h)
{
    CvScalar temp;
    int i,j;
    for(i=0;i<h;i++)
    {
        for(j=0;j<w;j++)
        {
            temp.val[0]=Ir[i*w+j];
            temp.val[1]=Ig[i*w+j];
            temp.val[2]=Ib[i*w+j];
            cvSet2D(img,i,j,temp);
        }
    }
}
void average_filter(IplImage *I)
{
    int i,j,k,l1,l2,count=0;
    double *s=(double*)calloc(I->nChannels,sizeof(double));
    CvScalar temp;
    for(i=0;i<I->height;i++)
    {
        for(j=0;j<I->width;j++)
        {
            for(k=0;k<I->nChannels;k++)
                s[k]=0;
            count=0;
            for(l1=-1;l1<=1;l1++)
            {
                for(l2=-1;l2<=1;l2++)
                {
                    if( ((i+l1)<0) || ((i+l1)>=I->height) || ((j+l2)<0) || ((j+l2)>=I->width))
                        continue;
                    count++;
                    temp=cvGet2D(I,i+l1,j+l2);
                    for(k=0;k<I->nChannels;k++)
                    {
                        s[k]=s[k]+temp.val[k];
                    }
                }
            }
            for(k=0;k<I->nChannels;k++)
            {
                s[k]=s[k]/(double)count;
                s[k]=(s[k]>255)?255:s[k];
                s[k]=(s[k]<0)?0:s[k];
                temp.val[k]=(int)s[k];
            }
            cvSet2D(I,i,j,temp);
        }
    }
}
void modify(double *I,int w,int h)
{
    double c=0;
    int i,j;
    for(i=0;i<h;i++)
    {
        for(j=0;j<w;j++)
        {
            //c=1/(1+ pow(I[i*w+j],2)/3600 + 0.0000001);
            c=exp(-pow(I[i*w+j],2)/3600);
            I[i*w+j]*=c;
            I[i*w+j]*=c;
        }
    }
}

void generate_new(double *I,double *IIncN,double *IIncS,double *IIncE,double *IIncW,double *Pimg,double *I0,int w,int h)
{
    int i,j;
    for(i=0;i<h;i++)
    {
        for(j=0;j<w;j++)
        {
            I[i*w+j]+=(0.15*(double)(IIncN[i*w+j]+IIncS[i*w+j]+IIncE[i*w+j]+IIncW[i*w+j]-Pimg[i*w+j]+I0[i*w+j]));
        }
    }
}
void copy_image(IplImage *img1,IplImage *img2)
{
    int i,j;
    CvScalar temp;
    for(i=0;i<img1->width;i++)
    {
        for(j=0;j<img1->height;j++)
        {
            temp=cvGet2D(img1,j,i);
            cvSet2D(img2,j,i,temp);
        }
    }
}

void find_increment(double *I,double *IIncrementN,double *IIncrementS,double *IIncrementE,double *IIncrementW,int owidth,int oheight)
{
    double *I1=(double*)calloc((owidth+2)*(oheight+2),sizeof(double));
    int i,j;
    for(i=1;i<=oheight;i++)
    {
        for(j=1;j<=owidth;j++)
        {
            I1[i*(owidth+2)+j]=I[(i-1)*owidth + (j-1)];
        }
    }
    for(i=1;i<=oheight;i++)
    {
        for(j=1;j<=owidth;j++)
        {
            IIncrementN[(i-1)*owidth+j-1] = I1[(i-1)*(owidth+2) + j]-I1[i*(owidth+2) + j];
            IIncrementS[(i-1)*owidth+j-1] = I1[(i+1)*(owidth+2) + j]-I1[i*(owidth+2) + j];
            IIncrementW[(i-1)*owidth+j-1] = I1[i*(owidth+2) + j-1]-I1[i*(owidth+2) + j];
            IIncrementE[(i-1)*owidth+j-1] = I1[i*(owidth+2) + j+1]-I1[i*(owidth+2) + j];
        }
    }
}
void pde(IplImage* img, IplImage* imgfinal, int owidth, int oheight)
{
    int i,j,k,l;
    int width=img->width;
    int height=img->height;
    double horizontal_scale=(double)(owidth)/(double)(img->width);
	double vertical_scale=(double)(oheight)/(double)(img->height);
    IplImage *img_initial = cvCreateImage(cvSize(owidth,oheight),IPL_DEPTH_8U,img->nChannels);
    IplImage *Pimg = cvCreateImage(cvSize(owidth,oheight),IPL_DEPTH_8U,img->nChannels);
    IplImage *imgtemp = cvCreateImage(cvSize(owidth,oheight),IPL_DEPTH_8U,img->nChannels);
    bilinear(img,img_initial,horizontal_scale,vertical_scale);
    //cvSaveImage("bilinear.jpg",img_initial);
    copy_image(img_initial,imgtemp);
    copy_image(img_initial,Pimg);
    
    average_filter(Pimg);
    double *Ir0=(double*)calloc(owidth*oheight,sizeof(double));
    double *Ig0=(double*)calloc(owidth*oheight,sizeof(double));
    double *Ib0=(double*)calloc(owidth*oheight,sizeof(double));
    double *Pimgr=(double*)calloc(owidth*oheight,sizeof(double));
    double *Pimgg=(double*)calloc(owidth*oheight,sizeof(double));
    double *Pimgb=(double*)calloc(owidth*oheight,sizeof(double));
    double *Ir=(double*)calloc(owidth*oheight,sizeof(double));
    double *Ig=(double*)calloc(owidth*oheight,sizeof(double));
    double *Ib=(double*)calloc(owidth*oheight,sizeof(double));
    getrgb(img_initial,Ir,Ig,Ib);
    getrgb(img_initial,Ir0,Ig0,Ib0);
    getrgb(Pimg,Pimgr,Pimgg,Pimgb);
    
    double *IrIncrementN=(double*)calloc(owidth*oheight,sizeof(double));
    double *IrIncrementS=(double*)calloc(owidth*oheight,sizeof(double));
    double *IrIncrementE=(double*)calloc(owidth*oheight,sizeof(double));
    double *IrIncrementW=(double*)calloc(owidth*oheight,sizeof(double));
    
    double *IgIncrementN=(double*)calloc(owidth*oheight,sizeof(double));
    double *IgIncrementS=(double*)calloc(owidth*oheight,sizeof(double));
    double *IgIncrementE=(double*)calloc(owidth*oheight,sizeof(double));
    double *IgIncrementW=(double*)calloc(owidth*oheight,sizeof(double));
    
    double *IbIncrementN=(double*)calloc(owidth*oheight,sizeof(double));
    double *IbIncrementS=(double*)calloc(owidth*oheight,sizeof(double));
    double *IbIncrementE=(double*)calloc(owidth*oheight,sizeof(double));
    double *IbIncrementW=(double*)calloc(owidth*oheight,sizeof(double));
    
    //single differentiation
 
    
    
    for(i=0;i<20;i++)
    {
        
        find_increment(Ir,IrIncrementN,IrIncrementS,IrIncrementE,IrIncrementW,owidth,oheight);
        find_increment(Ig,IgIncrementN,IgIncrementS,IgIncrementE,IgIncrementW,owidth,oheight);
        find_increment(Ib,IbIncrementN,IbIncrementS,IbIncrementE,IbIncrementW,owidth,oheight);
        
        modify(IrIncrementN,owidth,oheight);
        modify(IrIncrementS,owidth,oheight);
        modify(IrIncrementE,owidth,oheight);
        modify(IrIncrementW,owidth,oheight);
        
        modify(IgIncrementN,owidth,oheight);
        modify(IgIncrementS,owidth,oheight);
        modify(IgIncrementE,owidth,oheight);
        modify(IgIncrementW,owidth,oheight);
        
        modify(IbIncrementN,owidth,oheight);
        modify(IbIncrementS,owidth,oheight);
        modify(IbIncrementE,owidth,oheight);
        modify(IbIncrementW,owidth,oheight);
        
        
        generate_new(Ir,IrIncrementN,IrIncrementS,IrIncrementE,IrIncrementW,Pimgr,Ir0,owidth,oheight);
        generate_new(Ig,IgIncrementN,IgIncrementS,IgIncrementE,IgIncrementW,Pimgg,Ig0,owidth,oheight);
        generate_new(Ib,IbIncrementN,IbIncrementS,IbIncrementE,IbIncrementW,Pimgb,Ib0,owidth,oheight);
    }
    generateImgFrmRGB(imgfinal,Ir,Ig,Ib,owidth,oheight);
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
        pde(in_image,out_image,owidth,oheight);
        cvWriteFrame(writer,out_image);  
        printf("%d\n",i++);    
    }
 
    cvReleaseVideoWriter(&writer);
    cvReleaseCapture(&capture);
}

