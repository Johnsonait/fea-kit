#pragma once


//Class to implement Gaussian Quadrature on an arbitrary function that takes and returns doubles
class Quadrature
{
private:

public:
	Quadrature();


	//Gaussian Quadrature for one-dimensional integration
	double Integrate(const int& points, std::function<double(double)> Func);
	//Gaussian Quadrature for two-dimensional integration
	double Integrate(const int& points, std::function<double(double, double)> Func);
	//Gaussian Quadrature for three-dimensional integration
	double Integrate(const int& points, std::function<double(double, double, double)> Func);
};