#include "element.h"

//Constructors
Element::Element() : jacobian({}), jacobian_det(0) {}

Element::Element(const std::vector<double>&n, const std::vector<uint32_t>& el) : jacobian({}), jacobian_det(0)
{
	AddNode(n);
	SetGlobalElements(el);
}

//Useful
double Element::ShapeFunction(const double& zeta,const double& eta,const  double& mu, const uint32_t& m)
{
	return 0;
}

double Element::JacobianDet()
{
	return 0;
}

double Element::ShapeFunctionDerivatives(const double& zeta, const double& eta, const  double& mu, const uint32_t& m)
{
	return 0;
}

void CalcGlobalShapeDerivatives(const double& zeta, const double& eta, const  double& mu, const uint32_t& m)
{

}

//Accessors
const std::vector<std::vector<double>>& Element::GetGlobalShapeDerivatives()
{
	return global_shape_derivatives;
}

const std::vector<std::vector<double>>& Element::GetNodes()
{
	return nodes;
}

const double& Element::GetJacobianDet()
{
	return jacobian_det;
}

//Mutators
void Element::AddNode(const std::vector<double>& n)
{
	nodes.push_back(n);
}

void Element::SetGlobalElements(const std::vector<uint32_t>& el)
{
	global_nodes = el;
}
