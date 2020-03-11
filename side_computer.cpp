#define _USE_MATH_DEFINES
#include <cmath>
#include <stdio.h>
#include <unistd.h>
#include <sys/socket.h>
#include <bluetooth/bluetooth.h>
#include <bluetooth/rfcomm.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "opencv2/imgproc.hpp"
#include <opencv2/videoio.hpp>
#include <stdio.h>
#include <Eigen/Dense> 
#include <iostream>
#include <string>  
using namespace std;
using namespace cv;
using namespace std;
using namespace Eigen;

double ramp_function(double x){
    double result;
    if(x<-0.5){
        result=0;
    }else if(x>0.5){
        result=1;
    }else{
        result=x+0.5;
    }
    return result;
}
double jump_function(double x){
    double result;
    if(x<=0){
        result=0;
    }else{
        result=1;
    }
    return result;
}
void reverse(char* str, int len){ 
    int i = 0, j = len - 1, temp; 
    while (i < j) { 
        temp = str[i]; 
        str[i] = str[j]; 
        str[j] = temp; 
        i++; 
        j--; 
    } 
} 
int intToStr(int x, char str[], int d){ 
    int i = 0; 
    while (x) { 
        str[i++] = (x % 10) + '0'; 
        x = x / 10; 
    } 
    while (i < d) 
        str[i++] = '0'; 
    reverse(str, i); 
    str[i] = '\0'; 
    return i; 
}  
void ftoa(float n, char* res, int afterpoint){ 
    int ipart = (int)n; 
    float fpart = n - (float)ipart; 
    int i = intToStr(ipart, res, 0); 
    if (afterpoint != 0) { 
        res[i] = '.'; // add dot 
        fpart = fpart * pow(10, afterpoint); 
        intToStr((int)fpart, res + i + 1, afterpoint); 
    } 
} 

int main( int argc, char** argv )
{
    int red, green, blue;
    int i2, j2;
    int reduced_rows, reduced_cols;
    double comparison_neuron, result_index, result_value;
    int in1, in2;
    int iterator_cols;
    float element;
    char s_element[4];
    unsigned int microseconds = 1000;
    // -----------------------For capturing the image from the camera
	VideoCapture cap = cv::VideoCapture("http://172.23.192.45:4747/mjpegfeed?128x96");
	if(!cap.isOpened()){
		cout << "Error opening video stream or file" << endl;
		return -1;
	}
    // // -----------------------For processing the image
    reduced_rows = 96/4;                               // image 128x96
    reduced_cols = 128/4;                                // image 128x96
    MatrixXd layer_1(reduced_rows,reduced_cols); 
    MatrixXd layer_2(1,reduced_rows);
    MatrixXd layer_3, new_layer_3(1,reduced_cols);
    VectorXd layer_3_in, new_layer_3_in(reduced_cols);
    // // ------------------------For connected to the bluetooth device
    struct sockaddr_rc addr = { 0 };
    int s, status;
    char dest[18] = "98:D3:33:81:14:22";
    s = socket(AF_BLUETOOTH, SOCK_STREAM, BTPROTO_RFCOMM);
    addr.rc_family = AF_BLUETOOTH;
    addr.rc_channel = (uint8_t) 1;
    str2ba( dest, &addr.rc_bdaddr );
    status = connect(s, (struct sockaddr *)&addr, sizeof(addr));
    // //-------------------------For processing the image 
	while(1){
        //------------------------To read the image from the camera
		Mat image;
		cap >> image;
		if (image.empty())
		break;
        //----------------------To process the images
        i2=0;
        for (int i=0; i<image.rows; i++){
            j2=0;
            for (int j=0; j<image.cols; j++){
                if((i2*10)+5==i && (j2*10)+5==j){
                    blue = (uchar) image.at<Vec3b>(i,j)[0];
                    green = (uchar) image.at<Vec3b>(i,j)[1];
                    red = (uchar) image.at<Vec3b>(i,j)[2];
                    layer_1(i2,j2)=ramp_function(-0.002*red + 0.002*green -0.002*blue);       // 1/255
                    j2++;
                    if(j2==reduced_cols){
                        i2++;
                    }
                }
            }
        }
        layer_2 = MatrixXd::Constant(1,reduced_rows,0.04);              // reduccion de columnas
        layer_3 = (layer_2*layer_1);
        layer_3_in = VectorXd::LinSpaced(32, 1, 125)*0.007;             // 1/128
        iterator_cols = reduced_cols/2;
        for(int k=0; k<5; k++){
            for (int i=0; i<iterator_cols; i++){
                in1 = i*2;
                in2 = in1+1;
                comparison_neuron = jump_function(layer_3(in2)-layer_3(in1));
                result_value = ramp_function(ramp_function(layer_3(in1)-2*comparison_neuron-0.5)+ramp_function(2*comparison_neuron+layer_3(in2)-2.5)-0.5);
                result_index = ramp_function(ramp_function(layer_3_in(in1)-2*comparison_neuron-0.5)+ramp_function(2*comparison_neuron+layer_3_in(in2)-2.5)-0.5);
                new_layer_3(i) = result_value;
                new_layer_3_in(i) = result_index;
            }
            layer_3 = new_layer_3;
            layer_3_in = new_layer_3_in;
            iterator_cols = iterator_cols/2;
        }
        element = (layer_3_in(0)*142.85)/128;
        ftoa(element, s_element, 2); 
        printf("point:%f \n",element);
        // //----------------------------For sending the data to the bluetooth device
        status = write(s,s_element,4);
        // status = write(s, "\n", 2);
        //-----------------------------For showing the image
        usleep(microseconds);
		imshow( "Frame", image );
		char c=(char)waitKey(25);
		if(c==27)
		break;
	}
	cap.release();
	destroyAllWindows();
    // close(s);
	return 0;
}

// export PKG_CONFIG_PATH=/usr/lib/x86_64-linux-gnu/pkgconfig
// export PATH=/usr/lib/x86_64-linux-gnu/pkgconfig/opencv.pc:/home/joseph/Apps/Android-Studio/Project/tools/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin
// g++ -I /usr/include/eigen3/ -o union union.cpp -lbluetooth `pkg-config --cflags --libs opencv`