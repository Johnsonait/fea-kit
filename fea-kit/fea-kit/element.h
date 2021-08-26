#pragma once

#include <string>
#include <vector>

#include "linearsystem.h"
#include "matrix.h"

class Element
{
private:

	virtual double JacobianDet();

protected:
	std::vector<std::vector<double>> jacobian; //Store the 3x3 jacobian matrix
	double jacobian_det; //Store element jacobian value (can potentially change in space)
	std::vector<std::vector<double>> global_shape_derivatives; //Store global shape derivatives for each node
	std::vector<std::vector<double>> nodes;
	std::vector<std::vector<uint32_t>> bounds;
	std::vector<uint32_t> global_nodes; //Store global node indices

public:

	Element();
	Element(const std::vector<double>&, const std::vector<uint32_t>&);

	//Useful
	virtual double ShapeFunction(const double&, const double&, const  double&, const uint32_t&) = 0;
	virtual double ShapeFunctionDerivatives(const double&, const double&, const  double&, const uint32_t&, const uint32_t&);
	virtual void CalcGlobalShapeDerivatives(const double& zeta, const double& eta, const  double& mu, const uint32_t& m);
	virtual Matrix& Integrate(const int& points, std::function<Matrix& (double, double, double, std::shared_ptr<Element>, LinearElasticSolids*)> func, const Matrix& mat, std::shared_ptr<Element>, LinearElasticSolids*);


	//Accessors
	const std::vector<std::vector<double>>& GetGlobalShapeDerivatives();
	const std::vector<std::vector<double>>& GetNodes();
	const std::vector<uint32_t>& GetGlobalIDs();
	const double& GetJacobianDet();
	const std::vector<std::vector<uint32_t>>& GetBounds();

	virtual const double& GetJacobianDet(double,double,double);
	//Mutators
	void AddNode(const std::vector<double>& n);
	void SetGlobalElements(const std::vector<uint32_t>& el);
};

