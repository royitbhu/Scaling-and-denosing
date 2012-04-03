/*
Usage:	
	Compile :g++ nearest-neighbour.cpp -I /usr/local/include/opencv -L /usr/local/lib -lml -lcv -lhighgui -lcvaux -lcxcore 
	Run 	:./a.out image-file-name
Result is saved in a file "result.jpg".
Original image is saved in file "original.jpg"
*/

#include<iostream>
#include <math.h>
#include <stdio.h>
#include <cv.h>
#include <highgui.h>
using namespace std;

void nearest_neighbour(IplImage * img, IplImage* imgfinal, double horizontal_scale, double vertical_scale)
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
            double x = (double)j/(double)horizontal_scale;
            double y = (double)i/(double)vertical_scale;
            int temp_width,temp_height;
            if((x - (double)floor(x)) > 0.5)
                temp_width=ceil(x);
            else
                temp_width=floor(x);
            if((y - (double)floor(y)) > 0.5)
                temp_height=ceil(y);
            else
                temp_height=floor(y);
            for(k=0;k<new_channels;k++)
            {
                v.val[k]=orig_data[(int)(temp_height*orig_step + temp_width*orig_channels+k)];
                //new_data[(int)(i*new_step+j*new_channels+k)]=orig_data[(int)(temp_height*orig_step + temp_width*orig_channels+k)];
            }
            cvSet2D(imgfinal,i,j,v);
        }
    }
}

int main(int argc,char *argv[])
{
    IplImage* img=0;
	double horizontal_scale,vertical_scale;				
	if(argc<2)				//Check for the format of the command line input
	{
		printf("Usage : ./a.out <image-file-name>\n");
		exit(0);
	}
	img=cvLoadImage(argv[1]);		//Load the given image
  	if(!img){
    	printf("Could not load image file: %s\n",argv[1]);
    	exit(0);
  	}
	printf("Enter the scaling factors(m n) :");
	scanf("%lf%lf",&horizontal_scale,&vertical_scale);
	
  	IplImage * imgfinal;
  	imgfinal = cvCreateImage(cvSize(horizontal_scale*img->width, vertical_scale*img->height),img->depth,img->nChannels);//Create and initialize the new image
  	nearest_neighbour(img, imgfinal, horizontal_scale,vertical_scale);
  	  // show the image
  	cvSaveImage("result.jpg",imgfinal);			//Save the result image
  	cvNamedWindow("Result", CV_WINDOW_AUTOSIZE);
  	cvMoveWindow("Result", 100, 100);
  	cvShowImage("Result", imgfinal);			//Show the result image
	cvNamedWindow("Original", CV_WINDOW_AUTOSIZE);
  	cvMoveWindow("Original", 100, 100);
  	cvShowImage("Original", img);				//Show the original image
  	cvSaveImage("original.jpg",img);			//Save the original image
  	// wait for a key
  	cvWaitKey(0);

  	// release the image
  	cvReleaseImage(&imgfinal);				//Free the memory
  	cvReleaseImage(&img);
}
