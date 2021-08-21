#pragma once
#include <vector>
#include <functional>

#include "matrix.h"
#include "element.h"
//Class to implement Gaussian Quadrature on an arbitrary function that takes and returns doubles or matrices of doubles
class Quadrature
{
private:

public:
	Quadrature();

	//Gaussian Quadrature for one-dimensional integration
	double Integrate(const int& points, std::function<double(double)> func);
	//Gaussian Quadrature for two-dimensional integration
	double Integrate(const int& points, std::function<double(double, double)> func);
	//Gaussian Quadrature for three-dimensional integration
	double Integrate(const int& points, std::function<double(double, double, double)> func);

	//Gaussian Quadrature for integration of matrices on 1-d fields
	Matrix& Integrate(const int& points, std::function<Matrix& (double, Matrix)> func, const Matrix& mat);
	//Gaussian Quadrature for integration of matrices on 2-d fields
	Matrix& Integrate(const int& points, std::function<Matrix& (double, double, Matrix)> func, const Matrix& mat);
	
};