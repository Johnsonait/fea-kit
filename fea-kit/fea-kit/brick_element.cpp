#include "brick_element.h"

static const std::vector<std::vector<double>> BRICK_QUADRATURE_POINTS = {
   {0},
   {-0.5773502691896257,0.5773502691896257},
   {0,-0.7745966692414834,0.7745966692414834},
   {0.6521451548625461,0.6521451548625461,0.3478548451374538,0.3478548451374538},
};

static const std::vector<std::vector<double>> BRICK_QUADRATURE_WEIGHTS = {
	{2},
	{1,1},
	{0.8888888888888888,0.5555555555555556,0.5555555555555556},
	{-0.3399810435848563,0.3399810435848563,-0.8611363115940526,0.8611363115940526},
};

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
	for (int i = 0; i < body_nodes.size(); i++)
	{
		nodes.push_back(body_nodes[nodal_ids[i]-1]);
		global_shape_derivatives.push_back({});
	}
	std::vector<std::vector<uint32_t>> temp_bounds = {
		{nodal_ids[0],nodal_ids[1],nodal_ids[3]},
		{nodal_ids[0],nodal_ids[1],nodal_ids[2]},
		{nodal_ids[0],nodal_ids[2],nodal_ids[3] },
		{nodal_ids[1],nodal_ids[2],nodal_ids[3]}
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
	std::vector<std::vector<double>> J = GetJacobian(xsi,eta,zeta);
	//Store jacobian cofactor matrix
	double cof[3][3] = {
		{(J[1][1] * J[2][2]) - (J[1][2] * J[2][1]),(J[1][2] * J[2][0]) - (J[1][0] * J[2][2]),(J[1][0] * J[2][1]) - (J[1][1] * J[2][0])},
		{(J[2][1] * J[0][2]) - (J[2][2] * J[0][1]),(J[2][2] * J[0][0]) - (J[2][0] * J[0][2]),(J[2][0] * J[0][1]) - (J[2][1] * J[0][0])},
		{(J[0][1] * J[1][2]) - (J[0][2] * J[1][1]),(J[0][2] * J[1][0]) - (J[0][0] * J[1][2]),(J[0][0] * J[1][1]) - (J[0][1] * J[1][0])}
	};
	double j_det = (J[0][0] * cof[0][0]) + (J[0][1] * cof[0][1]) + (J[0][2] * cof[0][2]);
	for (size_t a = 0; a < global_nodes.size(); a++)
	{
		global_shape_derivatives[a][0] = ((ShapeFunctionDerivatives(xsi, eta, zeta, a, 1) * cof[0][0]) * (ShapeFunctionDerivatives(xsi, eta, zeta, a, 2) * cof[0][1]) * (ShapeFunctionDerivatives(xsi, eta, zeta, a, 3) * cof[0][2])) / j_det;
		global_shape_derivatives[a][1] = ((ShapeFunctionDerivatives(xsi, eta, zeta, a, 1) * cof[1][0]) * (ShapeFunctionDerivatives(xsi, eta, zeta, a, 2) * cof[1][1]) * (ShapeFunctionDerivatives(xsi, eta, zeta, a, 3) * cof[1][2])) / j_det;
		global_shape_derivatives[a][2] = ((ShapeFunctionDerivatives(xsi, eta, zeta, a, 1) * cof[2][0]) * (ShapeFunctionDerivatives(xsi, eta, zeta, a, 2) * cof[2][1]) * (ShapeFunctionDerivatives(xsi, eta, zeta, a, 3) * cof[2][2])) / j_det;
	}
}

const std::vector<std::vector<double>>& BrickElement::GetJacobian(const double& xsi, const double& eta, const double& zeta)
{
	jacobian = Jacobian(xsi,eta,zeta);
	return jacobian;
}

Matrix& BrickElement::Integrate(const int& points, std::function<Matrix& (double, double, double, std::shared_ptr<Element>, LinearElasticSolids*)> func, const Matrix& mat, std::shared_ptr<Element>el_ptr, LinearElasticSolids* model)
{
	int index = points - 1;
	//Construct result matrix of the appropriate size
	Matrix result(mat.CountRows(), mat.CountCols());
	for (int i = 0; i < points; i++)
	{
		for (int j = 0; j < points; j++)
		{
			for (int k = 0; k < points; k++)
			{
				result = result + (func(BRICK_QUADRATURE_POINTS[index][i], BRICK_QUADRATURE_POINTS[index][j], BRICK_QUADRATURE_POINTS[index][k], el_ptr, model)
					* (BRICK_QUADRATURE_WEIGHTS[index][i] * BRICK_QUADRATURE_WEIGHTS[index][j] * BRICK_QUADRATURE_WEIGHTS[index][k]));
			}
		}
	}
	return result;
}