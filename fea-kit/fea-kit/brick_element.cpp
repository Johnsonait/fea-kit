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

	jacobian_det = JacobianDet();
}

BrickElement::~BrickElement() {}

const double& BrickElement::GetJacobianDet(double, double, double)
{
	return jacobian_det;
}
double BrickElement::ShapeFunction(const double& zeta, const double& eta, const double& mu, const uint32_t& m)
{
	TODO;
}
//Returns local shape function value at a point for a node m in a given direction
double BrickElement::ShapeFunctionDerivatives(const double& zeta, const double& eta, const  double& mu, const uint32_t& m, const uint32_t& direction)
{
	TODO;
}

//This updates the global_shape_derivatives with the values at a given (local coord) point
void BrickElement::CalcGlobalShapeDerivatives(const double& zeta, const double& eta, const double& mu, const size_t& node_num)
{
	//Each row of global_shape_derivatives has 3 values: dN/dx dN/dy dN/dz
	for (size_t m = 0; m < nodes.size(); m++)
	{
		std::vector<std::vector<double>> mat = { {},{},{} }; //Matrix for jacobian system solution
		std::vector<std::vector<double>> b = { {} };//Vector "solution" for jacobian system 

		// For all coordinates associated with nodes (3)
		// Seems iffy, double check this
		for (int i = 0; i < nodes[0].size(); i++) //Rows of jacobian matrix 
		{
			//Update b vector with chosen node shape derivatives
			b[0].push_back(ShapeFunctionDerivatives(0, 0, 0, m, i));

			for (int j = 0; j < nodes[0].size(); j++) //Columns of jacobian matrix
			{
				double sum = 0;
				for (int n = 0; n < nodes.size(); n++) //Number of nodes (x = x1*N1 + ... +xn*Nn)
				{
					sum += nodes[n][i] * ShapeFunctionDerivatives(0, 0, 0, n, j);
				}
				mat[i].push_back(sum);
			}
		}
		LinearSystem temp_system(mat, b);
		//Update globe_shape_derivatives with solution jacobian equation
		//Can create B matrix after getting global shape derivatives
		temp_system.Solve(global_shape_derivatives[m]); //Update global shape derivates
	}
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