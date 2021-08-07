#include <vector>
#include <functional>

#include "quadrature.h"

//Pre-calculated quadrature points 
static const std::vector<std::vector<double>> QUADRATURE_POINTS = {
	{0},
	{-0.5773502691896257,0.5773502691896257},
	{0,-0.7745966692414834,0.7745966692414834},
	{0.6521451548625461,0.6521451548625461,0.3478548451374538,0.3478548451374538},
};

//Pre-calculated quuadrature weights
static const std::vector<std::vector<double>> QUADRATURE_WEIGHTS = {
	{2},
	{1,1},
	{0.8888888888888888,0.5555555555555556,0.5555555555555556},
	{-0.3399810435848563,0.3399810435848563,-0.8611363115940526,0.8611363115940526},
};


//Gaussian Quadrature for one-dimensional integration
double Quadrature::IntegrateOneD(const int& points, std::function<double(double)> Func)
{
	int index = points - 1;
	double result = 0.0;
	for (int i = 0; i < points; i++)
	{
		result += QUADRATURE_WEIGHTS[index][i] * Func(QUADRATURE_POINTS[index][i]);
	}
	return result;
}
//Gaussian Quadrature for two-dimensional integration
double Quadrature::IntegrateTwoD(const int& points, std::function<double(double, double)> Func)
{
	int index = points - 1;
	double result = 0.0;
	for (int i = 0; i < points; i++)
	{
		for (int j = 0; j < points; j++)
		{
			result += QUADRATURE_WEIGHTS[index][i] * QUADRATURE_WEIGHTS[index][j]
				* Func(QUADRATURE_POINTS[index][i], QUADRATURE_POINTS[index][j]);
		}
	}
	return result;
}
//Gaussian Quadrature for three-dimensional integration
double Quadrature::IntegrateThreeD(const int& points, std::function<double(double, double, double)> Func)
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
					* Func(QUADRATURE_POINTS[index][i], QUADRATURE_POINTS[index][j], QUADRATURE_POINTS[index][k]);
			}
		}
	}
	return result;
}