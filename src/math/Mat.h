#pragma once
#include "..\sys.h"
#include <iostream>
#include "Vec.h"

template<uint16 dimM, uint16 dimN>
class Mat
{
public:
	Mat() {
		assert(dimM && dimN);
	}
	Mat(double first, ...) {
		assert(dimM && dimN);
		va_list argp;
		va_start(argp, first);
		for (uint16 i=0; i < dimM; i++)
			for (uint16 j=0; j < dimN; j++)
			{
				if (i == 0 && j == 0) data[i][j]=first;
				else data[i][j]=va_arg(argp, double);
			}
		va_end(argp);
	}
	~Mat() { }

	Vec<dimN>& operator[] (uint16 n);
	const Vec<dimN>& operator[] (uint16 n) const;
	Mat<dimN,dimM> transpose ();
	void display ();
	Mat operator+ (const Mat &op);
	Mat& operator+= (const Mat &op);
	Mat operator- ();
	Mat& operator= (const Mat &op);
	Mat<dimM,1>& operator= (const Vec<dimM> &op); //TODO: CHECK this operation
	Mat operator- (const Mat &op);
	Mat operator* (const double op);
	void Identity ();
	void zero ();
	double det();
	Vec<dimM> eigenvalues();
	Mat<dimM,dimN> inv(double det);
	Mat<dimM-1,dimN-1> cross_cut (uint16 cuti, uint16 cutj);
	string toString ();
	//double* ptr ();
	//friend функции
	friend std::ostream& operator<< <dimM,dimN> (std::ostream& stream, const Mat<dimM,dimN> &obj);
	template <uint16 dimM, uint16 dimN, uint16 dimM2, uint16 dimN2> friend Mat<dimM,dimN2> operator* (const Mat<dimM,dimN> &op1, const Mat<dimM2,dimN2> &op2);
	template <uint16 dimM, uint16 dimN, uint16 dimM2> friend Vec<dimM> operator* (const Mat<dimM,dimN> &op1, const Vec<dimM2> &op2);

private:
	Vec<dimN> data[dimM];
};

//---------operator<<----------------------------------------------------------
template<uint16 dimM, uint16 dimN> std::ostream &operator<<(std::ostream &stream, const Mat<dimM,dimN> &obj) {
	for (uint16 i = 0; i < dimM; i++) {
		stream << "[" << obj.data[i] << "]" << std::endl;
	}
	return stream;
}

//----------operator[]---------------------------------------------------------
template<uint16 dimM, uint16 dimN>
inline Vec<dimN>& Mat<dimM,dimN>::operator[] (uint16 n) {
		assert(n < dimM);
		return data[n];
}
template<uint16 dimM, uint16 dimN>
inline const Vec<dimN>& Mat<dimM,dimN>::operator[] (uint16 n) const {
	assert(n < dimM);
	return data[n];
}
//----------Identity ()---------------------------------------------------------
template<uint16 dimM, uint16 dimN>
void Mat<dimM,dimN>::Identity () {
	assert (dimM==dimN);
	for (uint16 i=0; i < dimM; i++)
		for (uint16 j=0; j < dimN; j++)
			if (i==j)
				data[i][j]=1.0f;
			else
				data[i][j]=0.0f;
}
//----------zero()---------------------------------------------------------
template<uint16 dimM, uint16 dimN>
void Mat<dimM,dimN>::zero() {
	memset(data,0,sizeof(double)*dimM*dimN); //некрасиво
}
//-------------display()------------------------------------------------------
template<uint16 dimM, uint16 dimN>
void Mat<dimM,dimN>::display() {
	for (unsigned int i = 0; i < dimM; i++) {
		std::cout << "[";
		data[i].display();
		std::cout << "]" << std::endl;
	} 
}
//----------transpose()-------------------------------------------------------
template<uint16 dimM, uint16 dimN>
Mat<dimN,dimM> Mat<dimM,dimN>::transpose () {
	Mat<dimN,dimM> p;
	for (uint16 i=0; i < dimM; i++)
		for (uint16 j=0; j < dimN; j++)
			p[j][i]=this->data[i][j];
	return p;
}
//----------operator+(Mat)---------------------------------------------------------
template<uint16 dimM, uint16 dimN>
Mat<dimM,dimN> Mat<dimM,dimN>::operator+ (const Mat<dimM,dimN> &op) {
	Mat<dimM,dimN> p;
	for (uint16 i=0; i < dimM; i++)
		for (uint16 j=0; j < dimN; j++)
			p[i][j]=this->data[i][j]+op.data[i][j];
	return p;
}
//----------operator-(Mat)---------------------------------------------------------
//template<uint16 dimM, uint16 dimN>
//Mat<dimM,dimN> Mat<dimM,dimN>::operator- (const Mat<dimM,dimN> &op) {
//	Mat<dimM,dimN> p;
//	for (uint16 i=0; i < dimM; i++)
//		for (uint16 j=0; j < dimN; j++)
//			p[i][j]=this->data[i][j]-op.data[i][j];
//	return p;
//}
//-----------operator*(double)----------------------------------------------------------
template<uint16 dimM, uint16 dimN>
Mat<dimM,dimN> Mat<dimM,dimN>::operator*(const double op)
{
	Mat<dimM,dimN> p;
	for (uint16 i=0; i < dimM; i++)
		for (uint16 j=0; j < dimN; j++)
			p[i][j]=this->data[i][j]*op;
	return p;
}
//----------operator+=---------------------------------------------------------
template<uint16 dimM, uint16 dimN>
Mat<dimM,dimN>& Mat<dimM,dimN>::operator+= (const Mat<dimM,dimN> &op) {
	for (uint16 i=0; i < dimM; i++)
		for (uint16 j=0; j < dimN; j++)
			this->data[i][j]+=op.data[i][j];
	return *this;
}
//----------operator=(Mat)---------------------------------------------------------
template<uint16 dimM, uint16 dimN>
Mat<dimM,dimN>& Mat<dimM,dimN>::operator= (const Mat<dimM,dimN> &op)
{
	for (uint16 i = 0; i < dimM; i++)
		this->data[i]=op.data[i];
	return *this;
}
//----------operator=(Vec)---------------------------------------------------------
template<uint16 dimM, uint16 dimN>
Mat<dimM,1>& Mat<dimM,dimN>::operator= (const Vec<dimM> &op)
{
	assert(dimN == 1); //только для матрицы-столбца
	for (uint16 i = 0; i < dimM; i++)
		this->data[i][0] = op[i];
	return *this;
}
//----------operator*(Mat)---------------------------------------------------------
template<uint16 dimM, uint16 dimN, uint16 dimM2, uint16 dimN2> Mat<dimM,dimN2> operator*(const Mat<dimM,dimN> &op1, const Mat<dimM2,dimN2> &op2) {
	assert(dimN == dimM2);
	Mat<dimM, dimN2> p;
	double element;
	double *ptr1 = (double*) op1.data;
	double *ptr2 = (double*) op2.data;
	for (uint16 i=0; i < dimM; i++)
		for (uint16 j=0; j < dimN2; j++)
		{
			element=0.0f;
			for (uint16 l=0; l < dimN; l++) element+=ptr1[i*dimN+l]*ptr2[l*dimN2+j];
				//(op1.data[i])[l]*(op2.data[l])[j];
			p[i][j]=element;
		}
	return p;


	//assert(dimN == dimM2);
	//Mat<dimM, dimN2> p;
	//double element;
	//for (uint16 i=0; i < dimM; i++)
	//	for (uint16 j=0; j < dimN2; j++)
	//	{
	//		element=0.0f;
	//		for (uint16 l=0; l < dimN; l++) element+=(op1.data[i])[l]*(op2.data[l])[j];
	//		p[i][j]=element;
	//	}
	//return p;
}
//-----------operator*(Vec)---------------------------------------------------------
template <uint16 dimM, uint16 dimN, uint16 dimM2> Vec<dimM> operator* (const Mat<dimM,dimN> &op1, const Vec<dimM2> &op2) {
	assert(dimN==dimM2);
	Vec<dimM> p;
	double el;
	for (uint16 i=0; i < dimM; i++)
	{
		el = 0.0f;
		for (uint16 j=0; j < dimN; j++)
			el += op1.data[i][j]*op2[j];
		p[i] = el;
	}
	return p;
}
//--------------------------------------------------------------------------------
template<uint16 dimM, uint16 dimN>
string Mat<dimM,dimN>::toString()
{
	string p;
	for (unsigned int i = 0; i < dimM; i++) {
		p+= "[";
		p+=data[i].toString();
		p+= "]";
		if (i!=dimM-1) p+= "\n";
	} 
	return p;
}
//-------------------------------------------------------------------------------
template<> double Mat<1,1>::det();
template<> double Mat<2,2>::det();
template<> double Mat<3,3>::det();

template<uint16 dimM, uint16 dimN>
double Mat<dimM,dimN>::det()
{
	//разложение по строке i=0
	assert(dimM == dimN);
	double det = 0.0;
	double koef = 1;
	for (uint16 j=0; j < dimN; j++)
	{
		Mat<dimM-1,dimM-1> tmp;
		det += koef*data[0][j]*cross_cut(0,j).det();
		koef *= -1;
	}
	return det;
}

template<> Mat<1, 1> Mat<1,1>::inv(double det);
template<> Mat<2, 2> Mat<2,2>::inv(double det);
template<> Mat<3, 3> Mat<3,3>::inv(double det);

template<uint16 dimM, uint16 dimN>
Mat<dimM, dimN> Mat<dimM,dimN>::inv(double det)
{
	//через матрицу алг. дополнений
	assert(dimM == dimN);
	Mat<dimM,dimN> C;
	for (uint16 i=0; i < dimM; i++)
		for (uint16 j=0; j < dimN; j++)
			C[i][j] = npow(-1,i+j)*cross_cut(i,j).det();  //CHECK
	C = C.transpose()*(1/det);
	return C;
}

template<uint16 dimM, uint16 dimN>
Mat<dimM-1,dimN-1> Mat<dimM,dimN>::cross_cut (uint16 cuti, uint16 cutj)
{
	Mat<dimM-1, dimN-1> tmp;
	for (uint16 ii=0; ii < dimM-1; ii++)
		for (uint16 jj=0; jj < dimN-1; jj++)
			tmp[ii][jj] = data[(ii>=cuti)?ii+1:ii][(jj>=cutj)?jj+1:jj];
	return tmp;
}
template<uint16 dimM, uint16 dimN>
Vec<dimM>  Mat<dimM,dimN>::eigenvalues()
{
	assert(dimM == 3 && dimN == 3); //TODO: пока только для матриц 3 на 3
	double I1 = data[0][0] + data[1][1] + data[2][2];
	double I2 = data[0][0]*data[1][1] + data[1][1]*data[2][2] + data[0][0]*data[2][2] - 
		data[0][1]*data[0][1] - data[0][2]*data[0][2] - data[1][2]*data[1][2];
	double I3 = data[0][0]*data[1][1]*data[2][2] + 2*data[0][1]*data[1][2]*data[0][2] - 
		data[0][0]*data[1][2]*data[1][2]-data[1][1]*data[0][2]*data[0][2]-data[2][2]*data[0][1]*data[0][1];

	double a = -I1;
	double b = I2;
	double c = -I3;
	 // код из http://ru.wikipedia.org/wiki/%D0%A2%D1%80%D0%B8%D0%B3%D0%BE%D0%BD%D0%BE%D0%BC%D0%B5%D1%82%D1%80%D0%B8%D1%87%D0%B5%D1%81%D0%BA%D0%B0%D1%8F_%D1%84%D0%BE%D1%80%D0%BC%D1%83%D0%BB%D0%B0_%D0%92%D0%B8%D0%B5%D1%82%D0%B0
	// для трех вещественных корней
	// x*x*x + a * x * x + b * x  + c == 0
	double p = b - a * a / 3;
	double q = 2 * a * a * a / 27 - a * b / 3 + c;
	double A = sqrt(- 4 * p / 3); 
 
	double c3phi = - 4 * q / (A * A * A);
 
	double phi = acos(c3phi) / 3;  
 
	double root1 = A * cos(phi) - a / 3;
	double root2 = A * cos(phi + 2 * M_PI / 3) - a / 3;
	double root3 = A * cos(phi - 2 * M_PI / 3) - a / 3;

	//сортируем
	double s1,s2,s3;
	if (root1 > root2 && root1 > root3)
	{
		s1 = root1;
		if (root2 > root3)
		{
			s2 = root2;
			s3 = root3;
		}
		else
		{
			s2 = root3;
			s3 = root2;
		}
	}
	else if (root2 > root1 && root2 > root3)
	{
		s1 = root2;
		if (root1 > root3)
		{
			s2 = root1;
			s3 = root3;
		}
		else
		{
			s2 = root3;
			s3 = root1;
		}
	}
	else
	{
		s1 = root3;
		if (root1 > root2)
		{
			s2 = root1;
			s3 = root2;
		}
		else
		{
			s2 = root2;
			s3 = root1;
		}
	}
	return Vec<3>(s1,s2,s3);
}

//template<uint16 dimM, uint16 dimN>
//double* Mat<dimM,dimN>::ptr ()
//{
//	return data[0].ptr();
//}

















class dMat;
class dMat_interface
{
public:
	dMat_interface ():ptr(NULL),M(0),N(0),row(0) 
	{	}
	double& operator[] (uint16 col)
	{
		assert(ptr);
		assert (col < N);
		assert (row < M);
		return ptr[row*N + col];
	}
	friend class dMat;
private:
	uint16 M, N;
	uint32 row;
	double *ptr;
};

//dynamic matrix to pass arbitrary matrix to functions as arguments
class dMat
{
public:
	dMat(uint16 dim_m, uint16 dim_n) : dimM(0),dimN(0),data(NULL)
	{
		if (dim_m && dim_n)
			resize(dim_m, dim_n);
	}
	dMat(uint16 dim_m, uint16 dim_n, double first, ...) : dimM(0),dimN(0),data(NULL) {
		if (dim_m && dim_n)
		{
			resize(dim_m, dim_n);
			va_list argp;
			va_start(argp, first);
			data[0] = first;
			for (uint16 i=1; i < dimM*dimN; i++)
				data[i]=va_arg(argp, double);
			va_end(argp);
		}
	}
	dMat(dMat &from) : dimM(0),dimN(0),data(NULL)
	{
		operator=(from);
	}
	dMat_interface& operator[] (uint16 row) {
		dmat_int.row = row;
		return dmat_int;
	}
	void resize (uint16 dim_m, uint16 dim_n)
	{
		if (data) delete[] data;
		data = new double[dim_m*dim_n];
		dmat_int.ptr = data;
		dmat_int.M = dim_m;
		dmat_int.N = dim_n;
		dimM = dim_m;
		dimN = dim_n;
		zeros();
	}

	void fill (double first, ...) {
		va_list argp;
		va_start(argp, first);
		data[0] = first;
		for (uint16 i=1; i < dimM*dimN; i++)
			data[i]=va_arg(argp, double);
		va_end(argp);
	}
	~dMat()
	{
		if (data) delete[] data;
		data = NULL;
	}

	void zeros ()
	{
		memset(data,0,sizeof(double)*dimM*dimN);
	}
	template <uint16 M, uint16 N> Mat<M,N> toMat();
	template <uint16 M> Vec<M> toVec(uint16 col=0);
	template <uint16 M, uint16 N> void cpMat(Mat<M,N> &mat);
	template <uint16 M> void cpVec(Vec<M> &vec, uint16 col=0);
	dMat& operator= (const dMat &from)
	{
		resize(from.dimM, from.dimN);
		memcpy(this->data, from.data, sizeof(double)*dimM*dimN);
		return *this;
	}

	uint16 dM ()
	{
		return dimM;
	}
	uint16 dN ()
	{
		return dimN;
	}
	friend std::ostream& operator<< (std::ostream& stream, dMat &obj);
private:
	uint16 dimM;
	uint16 dimN;
	double *data;
	dMat_interface dmat_int;
};


//копирует dMat в правый верхний угол Mat<M,N>
template <uint16 M, uint16 N>
Mat<M,N> dMat::toMat()
{
	assert(M >= dimM && N >= dimN); //Mat >= dMat
	Mat<M,N> tmp;
	for (uint16 i = 0; i < M; i++)
		for (uint16 j = 0; j < N; j++)
			tmp[i][j] = (*this)[i][j];
	return tmp;
}

template <uint16 M> 
Vec<M> dMat::toVec(uint16 col)
{
	assert(M >= dimM && dimN > col); //Vec >= dMat по M 
	Vec<M> tmp;
	for (uint16 i = 0; i < M; i++)
		tmp[i] = (*this)[i][col];
	return tmp;
}

template <uint16 M, uint16 N> 
void dMat::cpMat(Mat<M,N> &mat)
{
	assert(dimM >= M && dimN >= N); //dMat >= Mat
	for (uint16 i = 0; i < M; i++)
		for (uint16 j = 0; j < N; j++)
			(*this)[i][j] = mat[i][j];
}

//функция копирует вектор в столбец col матрицы dMat
// col - от нуля
template <uint16 M> 
void dMat::cpVec(Vec<M> &vec, uint16 col)
{
	assert(dimM >= M && dimN > col);
	for (uint16 i = 0; i < M; i++)
		(*this)[i][col] = vec[i];
}
