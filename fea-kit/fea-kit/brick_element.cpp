#include "brick_element.h"

static const std::vector<std::vector<double>> BRICK_COEFFS = {
	{-1,-1,-1},
	{1,-1,-1},
	{1,1,-1},
	{-1,1,-1},
	{-1,-1,1},
	{1,-1,1},
	{1,1,1},
	{-1,1,1}
};

BrickElement::BrickElement() = default;

BrickElement::BrickElement(const std::vector<std::vector<double>>& body_nodes, std::vector<uint32_t>& nodal_ids)
{
	//Note that body_nodes is a reference to ALL nodes in the body
	for (int i = 0; i < nodal_ids.size(); i++)
	{
		nodes.push_back(body_nodes[nodal_ids[i]-1]);
		global_shape_derivatives.push_back({ {}, {},{} });
	}
	global_nodes = nodal_ids;
	//Store the global nodes that create the six surfaces of the brick
	std::vector<std::vector<uint32_t>> temp_bounds = {
		{nodal_ids[0],nodal_ids[1],nodal_ids[5],nodal_ids[4]},
		{nodal_ids[1],nodal_ids[2],nodal_ids[6],nodal_ids[5]},
		{nodal_ids[2],nodal_ids[3],nodal_ids[7],nodal_ids[6]},
		{nodal_ids[3],nodal_ids[0],nodal_ids[4],nodal_ids[7]},
		{nodal_ids[0],nodal_ids[3],nodal_ids[2],nodal_ids[1]},
		{nodal_ids[4],nodal_ids[5],nodal_ids[6],nodal_ids[7]}
	};
	bounds = temp_bounds;
}

BrickElement::~BrickElement() {}


//Returns 3x3 jacobian matrix at a given (local) points
std::vector<std::vector<double>>& BrickElement::Jacobian(const double& xsi, const double& eta, const double& zeta)
{
	std::vector<std::vector<double>> ret;
	//Each row in jacobian is const x,y,or z
	for (size_t row = 0; row < nodes[0].size(); ++row)
	{
		ret.push_back({});
		//Each col in jacobian differentiates xsi, eta, or zeta
		for (size_t col = 0; col < nodes[0].size(); ++col)
		{
			double sum = 0;
			for (size_t a = 0; a < global_nodes.size(); ++a)
			{
				sum += ShapeFunctionDerivatives(xsi,eta,zeta,a,col+1) * nodes[a][row];
			}
			ret[row].push_back(sum);
		}
	}
	return ret;
}

void BrickElement::Jacobian(const double& xsi, const double& eta, const double& zeta, std::vector<std::vector<double>>& ret)
{
	for (size_t row = 0; row < nodes[0].size(); ++row)
	{
		ret.push_back({});
		//Each col in jacobian differentiates xsi, eta, or zeta
		for (size_t col = 0; col < nodes[0].size(); ++col)
		{
			double sum = 0;
			for (size_t a = 0; a < global_nodes.size(); ++a)
			{
				sum += ShapeFunctionDerivatives(xsi, eta, zeta, a, col + 1) * nodes[a][row];
			}
			ret[row].push_back(sum);
		}
	}
}

double BrickElement::ShapeFunction(const double& xsi, const double& eta, const double& zeta, const uint32_t& m)
{
	return ((1+BRICK_COEFFS[m][0]*xsi)*(1+BRICK_COEFFS[m][1]*eta)*(1+BRICK_COEFFS[m][2]*zeta)) / 8;
}
//Returns local shape function value at a point for a node m in a given direction
double BrickElement::ShapeFunctionDerivatives(const double& xsi, const double& eta, const  double& zeta, const uint32_t& m, const uint32_t& direction)
{
	switch (direction)
	{
	case 1:
		return (BRICK_COEFFS[m][0]) * (1 + BRICK_COEFFS[m][1] * eta) * (1 + BRICK_COEFFS[m][2] * zeta) / 8;
	case 2:
		return (BRICK_COEFFS[m][1]) * (1 + BRICK_COEFFS[m][0] * xsi) * (1 + BRICK_COEFFS[m][2] * zeta) / 8;
	case 3:
		return (BRICK_COEFFS[m][2]) * (1 + BRICK_COEFFS[m][0] * xsi) * (1 + BRICK_COEFFS[m][1] * eta) / 8;
	default:
		break;
	}
}

//This updates the global_shape_derivatives with the values at a given (local coord) point
void BrickElement::CalcGlobalShapeDerivatives(const double& xsi, const double& eta, const double& zeta)
{
	std::vector<std::vector<double>> J = {};
	Jacobian(xsi, eta, zeta, J);
	//Store jacobian cofactor matrix
	double cof[3][3] = {
		{(J[1][1] * J[2][2]) - (J[1][2] * J[2][1]),(J[1][2] * J[2][0]) - (J[1][0] * J[2][2]),(J[1][0] * J[2][1]) - (J[1][1] * J[2][0])},
		{(J[2][1] * J[0][2]) - (J[2][2] * J[0][1]),(J[2][2] * J[0][0]) - (J[2][0] * J[0][2]),(J[2][0] * J[0][1]) - (J[2][1] * J[0][0])},
		{(J[0][1] * J[1][2]) - (J[0][2] * J[1][1]),(J[0][2] * J[1][0]) - (J[0][0] * J[1][2]),(J[0][0] * J[1][1]) - (J[0][1] * J[1][0])}
	};
	jacobian_det = (J[0][0] * cof[0][0]) + (J[0][1] * cof[0][1]) + (J[0][2] * cof[0][2]);
	for (size_t a = 0; a < global_nodes.size(); a++)
	{
		global_shape_derivatives[a][0] = ((ShapeFunctionDerivatives(xsi, eta, zeta, a, 1) * cof[0][0]) + (ShapeFunctionDerivatives(xsi, eta, zeta, a, 2) * cof[0][1]) + (ShapeFunctionDerivatives(xsi, eta, zeta, a, 3) * cof[0][2])) / jacobian_det;
		global_shape_derivatives[a][1] = ((ShapeFunctionDerivatives(xsi, eta, zeta, a, 1) * cof[1][0]) + (ShapeFunctionDerivatives(xsi, eta, zeta, a, 2) * cof[1][1]) + (ShapeFunctionDerivatives(xsi, eta, zeta, a, 3) * cof[1][2])) / jacobian_det;
		global_shape_derivatives[a][2] = ((ShapeFunctionDerivatives(xsi, eta, zeta, a, 1) * cof[2][0]) + (ShapeFunctionDerivatives(xsi, eta, zeta, a, 2) * cof[2][1]) + (ShapeFunctionDerivatives(xsi, eta, zeta, a, 3) * cof[2][2])) / jacobian_det;
	}
}

const std::vector<std::vector<double>>& BrickElement::GetJacobian(const double& xsi, const double& eta, const double& zeta)
{
	return Jacobian(xsi, eta, zeta);
}

const double& BrickElement::GetJacobianDet(const double& xsi, const double& eta, const double& zeta)
{
	std::vector<std::vector<double>> J = {};
	Jacobian(xsi, eta, zeta,J);
	//Store jacobian cofactor matrix
	double cof[3][3] = {
		{(J[1][1] * J[2][2]) - (J[1][2] * J[2][1]),(J[1][2] * J[2][0]) - (J[1][0] * J[2][2]),(J[1][0] * J[2][1]) - (J[1][1] * J[2][0])},
		{(J[2][1] * J[0][2]) - (J[2][2] * J[0][1]),(J[2][2] * J[0][0]) - (J[2][0] * J[0][2]),(J[2][0] * J[0][1]) - (J[2][1] * J[0][0])},
		{(J[0][1] * J[1][2]) - (J[0][2] * J[1][1]),(J[0][2] * J[1][0]) - (J[0][0] * J[1][2]),(J[0][0] * J[1][1]) - (J[0][1] * J[1][0])}
	};

	jacobian_det = (J[0][0] * cof[0][0]) + (J[0][1] * cof[0][1]) + (J[0][2] * cof[0][2]);

	return jacobian_det;
}

