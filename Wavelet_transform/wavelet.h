#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

typedef unsigned char BYTE;
class Daubechies
{
public:
	double **m_TempImg;
	double **m_LowImg;
	double C0, C1, C2, C3;
	int nWidth;
	int nHeight;
	cv::Mat LL_Img;
	cv::Mat LH_Img;
	cv::Mat HH_Img;
	cv::Mat HL_Img;
public:
	void change_int(cv::Mat inImage);
	void Daub4b(cv::Mat inImage);
	void transform(int row, int n);
	void transpose(int nn);
	void convert_mat_to_array(cv::Mat inImage);
	Daubechies(cv::Mat input_img);
	~Daubechies();
};