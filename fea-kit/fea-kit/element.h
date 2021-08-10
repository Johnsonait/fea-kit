#pragma once

#include <string>
#include <vector>

#include "linearsystem.h"

class Element
{
private:

	virtual double JacobianDet();

protected:
	std::vector<std::vector<double>> jacobian;
	double jacobian_det;
	std::vector<std::vector<double>> global_shape_derivatives; //Store global shape derivatives for each node
	std::vector<std::vector<double>> nodes;
	std::vector<uint32_t> global_nodes; //Store global node indices

public:

	Element();
	Element(const std::vector<double>&, const std::vector<uint32_t>&);

	//Useful
	virtual double ShapeFunction(const double&, const double&, const  double&, const uint32_t&) = 0;
	virtual double ShapeFunctionDerivatives(const double&, const double&, const  double&, const uint32_t&);


	//Accessors
	virtual const std::vector<std::vector<double>>& GetGlobalShapeDerivatives();
	const std::vector<std::vector<double>>& GetNodes();
	const double& GetJacobianDet();

	//Mutators
	void AddNode(const std::vector<double>& n);
	void SetGlobalElements(const std::vector<uint32_t>& el);
};

