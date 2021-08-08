#pragma once

#include <string>
#include <vector>

#include "linearsystem.h"

class Element
{
private:

protected:
	std::vector<std::vector<double>> jacobian;
	double jacobian_det;
	std::vector<std::vector<double>> shape_derivatives;
	std::vector<std::vector<double>> global_shape_derivatives; //Store global shape derivatives for each node
	std::vector<std::vector<double>> nodes; //Array storing node coordinates as X Y Z

public:

	virtual double ShapeFunction(double, double, double, uint32_t);
	virtual const std::vector<std::vector<double>>& GetGlobalShapeDerivatives(double, double, double);
	virtual const std::vector<std::vector<double>>& GetNodes();
};

