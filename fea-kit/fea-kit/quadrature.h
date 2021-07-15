#pragma once


//Class to implement Gaussian Quadrature on an arbitrary function that takes and returns doubles
class Quadrature
{
public:
	//Gaussian Quadrature for one-dimensional integration
	double IntegrateOneD(const int& points, std::function<double(double)> Func);
	//Gaussian Quadrature for two-dimensional integration
	double IntegrateTwoD(const int& points, std::function<double(double, double)> Func);
	//Gaussian Quadrature for three-dimensional integration
	double IntegrateThreeD(const int& points, std::function<double(double, double, double)> Func);
};