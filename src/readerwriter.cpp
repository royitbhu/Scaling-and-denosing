#include "cv.h"
#include "highgui.h"
#include <stdio.h>

int main()
{
    CvCapture *capture = cvCreateFileCapture("test.avi");
    //CvCapture *capture = cvCreateCameraCapture(0);
    CvVideoWriter *writer = 0;
    double fps  = cvGetCaptureProperty(capture, CV_CAP_PROP_FPS);//30
    int frameW  = cvGetCaptureProperty(capture, CV_CAP_PROP_FRAME_WIDTH); // 800
    int frameH  = cvGetCaptureProperty(capture, CV_CAP_PROP_FRAME_HEIGHT); // 600
    
    writer=cvCreateVideoWriter("/home/saket/out.avi",cvGetCaptureProperty(capture, CV_CAP_PROP_FOURCC),//CV_FOURCC('P','I','M','1'),
                           fps,cvSize(frameW,frameH));

    if(capture==NULL)
    {
        printf("cap\n");
        return -1;
    }
    
    if(writer==NULL)
    {
        printf("uffffffffff\n");
       // return -1;
    }
    
    IplImage* img = 0;
    int i=0;
    printf("%lf %d %d\n",fps,frameH, frameW);
    while(cvGrabFrame(capture))
    {
        img=cvRetrieveFrame(capture);  
        cvWriteFrame(writer,img);      
        //cvSaveImage("output.jpg",img);
    }
    
    cvReleaseVideoWriter(&writer);
    cvReleaseCapture(&capture);
    cvReleaseImage(&img);

}

