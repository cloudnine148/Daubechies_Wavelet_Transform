#include "wavelet.h"

int main(void)
{
	cv::Mat input_img = cv::imread("test.jpg");
	Daubechies db(input_img);
	db.Daub4b(input_img);
}