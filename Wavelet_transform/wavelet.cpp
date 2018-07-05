#include "wavelet.h"
#include <math.h>

Daubechies::Daubechies(cv::Mat input_img)
{
	C0 = 0.4829629131445341;
	C1 = 0.8365163037378079;
	C2 = 0.2241438680420134;
	C3 = -0.1294095225512604;
	nWidth = input_img.cols;
	nHeight = input_img.rows;
	m_TempImg = new double*[nHeight];
	m_LowImg = new double*[nHeight];

	for (auto i = 0; i < nHeight; i++)
	{
		m_TempImg[i] = new double[nWidth];
		m_LowImg[i] = new double[nWidth];

		memset(m_TempImg[i], 0, sizeof(double)*nWidth);
		memset(m_LowImg[i], 0, sizeof(double)*nWidth);
	}


	LL_Img = cv::Mat::zeros(cv::Size(nHeight / 2, nWidth / 2), CV_8UC1);
	LH_Img = cv::Mat::zeros(cv::Size(nHeight / 2, nWidth / 2), CV_8UC1);
	HH_Img = cv::Mat::zeros(cv::Size(nHeight / 2, nWidth / 2), CV_8UC1);
	HL_Img = cv::Mat::zeros(cv::Size(nHeight / 2, nWidth / 2), CV_8UC1);
}
Daubechies::~Daubechies()
{
	for (auto i = 0; i < nHeight; i++)
	{
		delete[] m_TempImg[i];
		delete[] m_LowImg[i];
	}
	delete[] m_TempImg;
	delete[] m_LowImg;
}

void Daubechies::Daub4b(cv::Mat inImage)
{
	int i, j, n, k;

	cv::cvtColor(inImage, inImage, CV_BGR2GRAY);
	convert_mat_to_array(inImage);

	for (k = 0; k < nWidth; k++) {
		for (n = nWidth; n >= nWidth; n >>= 1) {
			transform(k, n);
		}
		for (i = 0; i < nWidth; i++) m_TempImg[k][i] = m_LowImg[k][i];
	}

	transpose(nWidth);

	for (k = 0; k < nWidth; k++) {
		for (n = nWidth; n >= nWidth; n >>= 1) {
			transform(k, n);
		}
		for (i = 0; i < nWidth; i++) m_TempImg[k][i] = m_LowImg[k][i];
	}

	transpose(nWidth);

	change_int(inImage);
}

void Daubechies::transform(int k, int n)
{
	int i, j;

	if (n >= 5) {
		int half = n >> 1;

		for (i = 0, j = 0; j < n - 3; j += 2, i++) {
			m_LowImg[k][i] = m_TempImg[k][j] * C0 + m_TempImg[k][j + 1] * C1 + m_TempImg[k][j + 2] * C2 + m_TempImg[k][j + 3] * C3;
			m_LowImg[k][i + half] = m_TempImg[k][j] * C3 - m_TempImg[k][j + 1] * C2 + m_TempImg[k][j + 2] * C1 - m_TempImg[k][j + 3] * C0;
		}
		m_LowImg[k][i] = m_TempImg[k][n - 2] * C0 + m_TempImg[k][n - 1] * C1 + m_TempImg[k][0] * C2 + m_TempImg[k][1] * C3;
		m_LowImg[k][i + half] = m_TempImg[k][n - 2] * C3 - m_TempImg[k][n - 1] * C2 + m_TempImg[k][0] * C1 - m_TempImg[k][1] * C0;
	}

}

void Daubechies::transpose(int nn)
{
	int i, j;
	double temp;

	for (i = 0; i < nn; i++) {
		for (j = 0; j < i; j++) {
			temp = m_TempImg[i][j];
			m_TempImg[i][j] = m_TempImg[j][i];
			m_TempImg[j][i] = temp;
		}
	}
}

void Daubechies::change_int(cv::Mat inImage)
{
	double max = 0;
	double min = 0;

	int x, y;
	
	BYTE tmpVal;
	BYTE* inImg = (BYTE*)inImage.data;
	BYTE* ll_img = (BYTE*)LL_Img.data;
	BYTE* lh_img = (BYTE*)LH_Img.data;
	BYTE* hh_img = (BYTE*)HH_Img.data;
	BYTE* hl_img = (BYTE*)HL_Img.data;

	for (x = 0; x < nHeight /2; x++)
		for (y = 0; y < nWidth/2; y++) {
			if (max < fabs(m_TempImg[x][y]))
				max = fabs(m_TempImg[x][y]);
		}
	for (x = 0; x < nHeight /2; x++)
		for (y = 0; y < nWidth/2; y++) {
			tmpVal = (BYTE)((fabs(m_TempImg[x][y]) / max) * 255.0);
			inImg[x*inImage.cols + y] = tmpVal;
			ll_img[x*LL_Img.cols + y] = tmpVal;
		}

	max = 0;

	for (x = 0; x < nHeight /2; x++)
		for (y = nWidth/2; y < nWidth; y++) {
			if (max < fabs(m_TempImg[x][y]))
				max = fabs(m_TempImg[x][y]);
		}

	for (x = 0; x < nHeight /2; x++)
		for (y = nWidth/2; y < nWidth; y++) {
			tmpVal = (BYTE)((fabs(m_TempImg[x][y]) / max) * 255.0);
			inImg[x*inImage.cols + y] = tmpVal;
			hl_img[x*HL_Img.cols + (y-nWidth/2)] = tmpVal;
		}

	max = 0;

	for (x = nHeight /2; x < nHeight; x++)
		for (y = 0; y < nWidth; y++) {
			if (max < fabs(m_TempImg[x][y]))
				max = fabs(m_TempImg[x][y]);
		}

	for (x = nHeight /2; x < nHeight; x++)
		for (y = 0; y < nWidth / 2; y++) {
			tmpVal = (BYTE)((fabs(m_TempImg[x][y]) / max) * 255.0);
			inImg[x*inImage.cols + y] = tmpVal;
			lh_img[(x- nHeight /2)*LH_Img.cols + y] = tmpVal;
		}
	max = 0;

	for (x = nHeight / 2; x < nHeight; x++)
		for (y = nWidth / 2; y < nWidth; y++) {
			if (max < fabs(m_TempImg[x][y]))
				max = fabs(m_TempImg[x][y]);
		}

	for (x = nHeight /2; x < nHeight; x++)
		for (y = nWidth / 2; y < nWidth; y++) {
			tmpVal = (BYTE)((fabs(m_TempImg[x][y]) / max) * 255.0);
			inImg[x*inImage.cols + y] = tmpVal;
			hh_img[(x- nHeight /2)*HH_Img.cols + (y-nWidth/2)] = tmpVal;
		}
	cv::resize(LL_Img, LL_Img, cv::Size(nHeight, nWidth), 0, 0, CV_INTER_NN);
	cv::resize(LH_Img, LH_Img, cv::Size(nHeight, nWidth), 0, 0, CV_INTER_NN);
	cv::resize(HH_Img, HH_Img, cv::Size(nHeight, nWidth), 0, 0, CV_INTER_NN);
	cv::resize(HL_Img, HL_Img, cv::Size(nHeight, nWidth), 0, 0, CV_INTER_NN);

	
	//cv::imshow("LL_Img", LL_Img);
	//cv::imshow("LH_Img", LH_Img);
	//cv::imshow("HL_Img", HL_Img);
	//cv::imshow("HH_Img", HH_Img);
	//cv::waitKey(0);
}

void Daubechies::convert_mat_to_array(cv::Mat inImage)
{
	BYTE* inImg = (BYTE*)inImage.data;

	for (auto i = 0; i < nHeight; i++)
	{
		for (auto j = 0; j < nWidth; j++)
		{
			m_TempImg[i][j] = inImg[i*inImage.cols + j];
		}
	}
}