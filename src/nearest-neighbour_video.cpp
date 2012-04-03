/*
Usage:	
	Compile :g++ nearest-neighbour.cpp -I /usr/local/include/opencv -L /usr/local/lib -lml -lcv -lhighgui -lcvaux -lcxcore 
	Run 	:./a.out image-file-name
Result is saved in a file "result.jpg".
Original image is saved in file "original.jpg"
*/

#include <math.h>
#include <cv.h>
#include <stdio.h>
#include <highgui.h>

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
        nearest_neighbour(in_image, out_image, h_scale, v_scale);
        cvWriteFrame(writer,out_image);  
        printf("%d\n",i++);    
    }
 
    cvReleaseVideoWriter(&writer);
    cvReleaseCapture(&capture);
}
