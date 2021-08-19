#include "quadrature.h"

static const std::vector<std::vector<double>> QUADRATURE_POINTS = {
   {0},
   {-0.5773502691896257,0.5773502691896257},
   {0,-0.7745966692414834,0.7745966692414834},
   {0.6521451548625461,0.6521451548625461,0.3478548451374538,0.3478548451374538},
};

static const std::vector<std::vector<double>> QUADRATURE_WEIGHTS = {
	{2},
	{1,1},
	{0.8888888888888888,0.5555555555555556,0.5555555555555556},
	{-0.3399810435848563,0.3399810435848563,-0.8611363115940526,0.8611363115940526},
};

//Default Constructor
Quadrature::Quadrature() = default;

//Gaussian Quadrature for one-dimensional integration
double Quadrature::Integrate(const int& points, std::function<double(double)> func)
{
	int index = points - 1;
	double result = 0.0;
	for (int i = 0; i < points; i++)
	{
		result += QUADRATURE_WEIGHTS[index][i] * func(QUADRATURE_POINTS[index][i]);
	}
	return result;
}
//Gaussian Quadrature for two-dimensional integration
double Quadrature::Integrate(const int& points, std::function<double(double, double)> func)
{
	int index = points - 1;
	double result = 0.0;
	for (int i = 0; i < points; i++)
	{
		for (int j = 0; j < points; j++)
		{
			result += QUADRATURE_WEIGHTS[index][i] * QUADRATURE_WEIGHTS[index][j]
				* func(QUADRATURE_POINTS[index][i], QUADRATURE_POINTS[index][j]);
		}
	}
	return result;
}
//Gaussian Quadrature for three-dimensional integration
double Quadrature::Integrate(const int& points, std::function<double(double, double, double)> func)
{
	int index = points - 1;
	double result = 0.0;
	for (int i = 0; i < points; i++)
	{
		for (int j = 0; j < points; j++)
		{
			for (int k = 0; k < points; k++)
			{
				result += QUADRATURE_WEIGHTS[index][i] * QUADRATURE_WEIGHTS[index][j] * QUADRATURE_WEIGHTS[index][k]
					* func(QUADRATURE_POINTS[index][i], QUADRATURE_POINTS[index][j], QUADRATURE_POINTS[index][k]);
			}
		}
	}
	return result;
}
//Gaussian Quadrature for integration of matrices on 1-d fields
Matrix& Integrate(const int& points, std::function<Matrix& (double, Matrix)> func, const Matrix& mat)
{
	int index = points - 1;
	//Construct result matrix of the appropriate size
	Matrix result(mat.CountRows(), mat.CountCols());
	for (int i = 0; i < points; i++)
	{
			result = result + (func(QUADRATURE_POINTS[index][i], mat)
				* (QUADRATURE_WEIGHTS[index][i]));
	}
	return result;
}

//Gaussian Quadrature for integration of matrices on 2-d fields
Matrix& Integrate(const int& points, std::function<Matrix& (double, double, Matrix)> func, const Matrix& mat)
{
	int index = points - 1;
	//Construct result matrix of the appropriate size
	Matrix result(mat.CountRows(), mat.CountCols());
	for (int i = 0; i < points; i++)
	{
		for (int j = 0; j < points; j++)
		{
			result = result + (func(QUADRATURE_POINTS[index][i], QUADRATURE_POINTS[index][j], mat)
				* (QUADRATURE_WEIGHTS[index][i] * QUADRATURE_WEIGHTS[index][j]));
		}
	}
	return result;
}

//Gaussian Quadrature for integration of matrices on 3-d fields
Matrix& Quadrature::Integrate(const int& points, std::function<Matrix& (double, double, double, Matrix, std::shared_ptr<Element>, LinearElasticSolids*)> func, const Matrix& mat, std::shared_ptr<Element> el_ptr, LinearElasticSolids* model)
{
	int index = points - 1;
	//Construct result matrix of the appropriate size
	Matrix result(mat.CountRows(),mat.CountCols());
	for (int i = 0; i < points; i++)
	{
		for (int j = 0; j < points; j++)
		{
			for (int k = 0; k < points; k++)
			{
				result = result + (func(QUADRATURE_POINTS[index][i], QUADRATURE_POINTS[index][j], QUADRATURE_POINTS[index][k],mat,el_ptr,model)
					*(QUADRATURE_WEIGHTS[index][i] * QUADRATURE_WEIGHTS[index][j] * QUADRATURE_WEIGHTS[index][k]));
			}
		}
	}
	return result;
}