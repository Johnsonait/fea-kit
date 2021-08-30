#include "linear_elastic_solids.h"


static const double GRAVITY = 9.81;

Matrix& StiffnessIntegrand(double eta, double zeta, double mu, std::shared_ptr<Element> el_ptr, LinearElasticSolids* model)
{
	const Matrix B = model->ConstructBMatrix(eta, zeta, mu, el_ptr);
	const Matrix B_T = B.GetTranspose();
	const Matrix C = model->GetElasticMatrix();

	//Final integrand value
	Matrix ret = B_T * C * B * el_ptr->GetJacobianDet();
	return ret;
}

//Finite element form of integrand used for updating local force vector to include body forces
Matrix& BodyForceIntegrand(double eta, double zeta, double mu, std::shared_ptr<Element> el_ptr, LinearElasticSolids* model)
{
	Matrix N = model->ConstructShapeMatrix(eta,zeta,mu,el_ptr);
	std::vector<std::vector<double>> temp = { {0,0,-GRAVITY} };
	Matrix body_force(temp);

	Matrix ret = (N*body_force)*el_ptr->GetJacobianDet();
	return ret;
}

Matrix& SurfaceForceIntegrand(double eta, double zeta, double mu, std::shared_ptr<Element> el_ptr, LinearElasticSolids* model)
{

	Matrix N = model->ConstructShapeMatrix(eta, zeta, mu, el_ptr);
	
}


void LinearElasticSolids::Lame() //Pronounced Lam-eh
{
	lambda = (E * poisson) / ((1 + poisson) * (1 - (2 * poisson)));
	G = E / (2 * (1 + poisson));
}

//Populates isotropic elasticity matrix with Lame parameters
void LinearElasticSolids::ConstructElasticMatrix()
{
	std::vector<std::vector<double>> temp = {
		{((2 * G) + lambda),lambda,lambda,0,0,0},
		{lambda,((2 * G) + lambda),lambda,0,0,0},
		{lambda,lambda,((2 * G) + lambda),0,0,0},
		{0,0,0,G,0,0},
		{0,0,0,0,G,0},
		{0,0,0,0,0,G} };
	Matrix mat(temp);
	elastic_matrix = mat;
}

void LinearElasticSolids::InitMatrices(std::shared_ptr<std::vector<std::vector<double>>> ptr,const uint32_t& rows, const uint32_t& cols)
{
	(*ptr) = {};
	for (uint32_t i = 0; i < rows; ++i)
	{
		ptr->push_back({});
		for (uint32_t j = 0; j < cols; ++j)
		{
			(*ptr)[i].push_back(0);
		}
	}
}

//Constructs 6x12 elemental B matrix
//Requires global shape function derivatives
Matrix LinearElasticSolids::ConstructBMatrix(const double& xsi,const double& eta, const double& zeta, std::shared_ptr<Element> el)
{
	std::vector<std::vector<double>> nodes = el->GetNodes();
	el->CalcGlobalShapeDerivatives(xsi, eta, zeta);
	std::vector<std::vector<double>> G = el->GetGlobalShapeDerivatives();

	std::vector<std::vector<double>> Mat; //Temporary matrix to be returned as matrix

	Mat.resize(6); //Set rows of Mat to 6, there are 6 rows in B matrix and dof*m_nodes number of columns

	for (int m = 0; m < nodes.size(); ++m)
	{
		std::vector<std::vector<double>> sub_matrix = {
				{G[m][0],0,0},
				{0,G[m][1],0},
				{0,0,G[m][2]},
				{0,G[m][2],G[m][1]},
				{G[m][2],0,G[m][0]},
				{G[m][1],G[m][0],0}
		};
		for (int i = 0; i < 6; i++) //6 rows
		{
			for (int j = 0; j < 3; j++) //For col of sub-matrix j<dof
			{
				Mat[i].push_back(sub_matrix[i][j]); //Add rows to full B matrix
			}
		}
	}
	Matrix ret(Mat);
	return ret;
}

//Generates a 3x3*m matrix storing useful version of global shape derivatives in a matrix
//Used to multiply properties by shape function during formulation of element force vector
Matrix LinearElasticSolids::ConstructShapeMatrix(const double& zeta, const double& eta, const double& mu, std::shared_ptr<Element> el)
{
	std::vector<std::vector<double>> nodes = el->GetNodes();
	el->CalcGlobalShapeDerivatives(zeta, eta, mu);
	std::vector<std::vector<double>> G = el->GetGlobalShapeDerivatives();

	std::vector<std::vector<double>> Mat; //Temporary matrix to be returned as matrix

	for (size_t m = 0; m < nodes.size();++m)
	{
		std::vector<std::vector<double>> sub_matrix = {
			{el->ShapeFunction(eta,zeta,mu,m),0,0},
			{0,el->ShapeFunction(eta,zeta,mu,m),0},
			{0,0,el->ShapeFunction(eta,zeta,mu,m)}
		};

		Mat.push_back(sub_matrix[0]);
		Mat.push_back(sub_matrix[1]);
		Mat.push_back(sub_matrix[2]);
	}
	Matrix ret(Mat);
	return ret;
}

//This function is used to update the provided local_k matrix using quadrature 
void LinearElasticSolids::CalculateLocalK(Matrix& local_k, std::shared_ptr<Element> el_ptr)
{
	local_k = el_ptr->Integrate(2, StiffnessIntegrand, local_k ,el_ptr, this);
}

//This function updates the global stiffness matrix from the element stiffness matrix
void LinearElasticSolids::AssembleStiffness(Matrix& local_k,const std::vector<uint32_t>& node_ids)
{
	for (size_t i = 0; i < local_k.CountRows(); ++i)
	{
		for (size_t j = 0; j < local_k.CountCols(); j++)
		{
			(*global_k)[node_ids[i] - 1][node_ids[j] - 1] += local_k[i][j];
		}
	}
}

void LinearElasticSolids::EnforceSurfaceBounds(Matrix& local_k, Matrix& local_f, std::shared_ptr<Element> el_ptr)
{
	for (size_t boundary = 0; boundary<body_ptr->GetBoundaryCount(); ++boundary)
	{
		std::string type = "";
		std::vector<uint32_t> boundary_nodes = {};
		std::vector<double> boundary_vector = {};
		body_ptr->GetBoundaryInfo(boundary_nodes,type,boundary_vector,boundary);

		if (type == "*TRACTION")
		{
			std::vector<uint32_t> temp_track;
			for (uint32_t node : boundary_nodes)
			{
				for (auto e_node : el_ptr->GetGlobalIDs())
				{
					if (node == e_node)
					{
						temp_track.push_back(node);
						break;
					}
				}
			}
			//We now have every element node associated with the boundary stored in temp_track.
			//Now we need to associate the nodes with element-defined bounds
			//GetBounds should return a vector<uint32_t> that stores the node IDs for each bound in the element (eg: a Tetrahedral element
			//has 4 faces and so four bounds containing 3 global node IDs each)
			std::vector<std::vector<uint32_t>> affected_bound = {};
			for (auto bound : el_ptr->GetBounds())
			{
				bool is_in = true;
				for (auto node : bound)
				{
					if (!binary_search(temp_track.begin(), temp_track.end(), node))
					{
						is_in = false;
					}
				}
				if (!is_in)
				{
					affected_bound.push_back(bound);
				}
			}
			//Now we have the bounds of the element that are affected by the traction vector 
			//(Which is stored in boundary_vector)
			//The next step is to integrate the bounds surfaces in affected_bound with the vector and shape matrices, then update the local_f
			//For now we will just use the bound area to find the total force and average that across the bound nodes

		}
	}
}

//This function is used to update the provided local_f matrix using quadrature 
//It needs to consider the volume integral over the body and the surface integral over the tractioned surfaces
void LinearElasticSolids::CalculateLocalForce(Matrix& local_f, std::shared_ptr<Element> el_ptr)
{
	CalculateBodyForce(local_f, el_ptr);
	CalculateSurfaceForce(local_f, el_ptr);
}

//Integrate the body-force values to be used in global force-vector
void LinearElasticSolids::CalculateBodyForce(Matrix& local_f, std::shared_ptr<Element> el_ptr)
{
	local_f = local_f + el_ptr->Integrate(2, BodyForceIntegrand, local_f, el_ptr, this);
}

//Integrate the surface tractions to be used in global force-vector
void LinearElasticSolids::CalculateSurfaceForce(Matrix& local_f, std::shared_ptr<Element> el_ptr)
{
	local_f = local_f + el_ptr->Integrate(2,SurfaceForceIntegrand, local_f,el_ptr, this);
}

LinearElasticSolids::LinearElasticSolids()//Default constructors
{
	E = 0;
	G = 0;
	lambda = 0;
	poisson = 0;
	body_ptr = nullptr;
	read_ptr = nullptr;

}

LinearElasticSolids::LinearElasticSolids(Reader& reader, Body& body)
{
	body_ptr = &body; 
	read_ptr = &reader;
	global_sol = std::make_shared<std::vector<std::vector<double>>>();//Prepare heap-allocated solution vector
	global_f = std::make_shared<std::vector<std::vector<double>>>();//Prepare heap-allocated f vector
	global_k = std::make_shared<std::vector<std::vector<double>>>();//Prepare heap-allocated K matrix

	uint32_t sz = (body_ptr->GetNodeNum()) * 3;
	InitMatrices(global_sol,sz,1); 
	InitMatrices(global_f, sz, 1);
	InitMatrices(global_k, sz, sz);//Create square matrix of zeros to be added to during solution;

	E = body_ptr->GetStiffness();
	poisson = body_ptr->GetPoisson();
	Lame();
	ConstructElasticMatrix();
}

Body& LinearElasticSolids::GetBody()
{
	return *body_ptr;
}

void LinearElasticSolids::EnforceDisplacements(std::shared_ptr<std::vector<std::vector<double>>> k, std::shared_ptr<std::vector<std::vector<double>>> f)
{
	for (size_t boundary = 0; boundary < body_ptr->GetBoundaryCount(); ++boundary)
	{
		std::string type = "";
		std::vector<uint32_t> boundary_nodes = {};
		std::vector<double> boundary_vector = {};
		body_ptr->GetBoundaryInfo(boundary_nodes, type, boundary_vector, boundary);
		if (type == "*DISPLACEMENT")
		{
			for (uint32_t node : boundary_nodes)
			{
				//Clear the given node's row and set the node stiffness to 1
				for (uint32_t col = 0; col < (*k)[node-1].size(); ++col)
				{
					(*k)[node - 1][col] = 0;
				}
				(*k)[node-1][node-1] = 1;
				//Set the Fx, Fy, Fz values to the given displacement so they show up in the solution
				for (size_t n = 0; n<3;++n)
				{
					(*f)[node - 1 + n][0] = boundary_vector[n];
				}
			}
		}
	}
}
//This is where the overall problem is solved
//Each element in the problem is used to first construct the element-wise stiffness matrix which is then assembled into the global stiffness matrix
//The global problem is then solved and the results stored (the result being the displacement values at each node)
void LinearElasticSolids::Solve()
{
	const std::vector<std::vector<double>>& nodes = body_ptr->GetNodes();
	const std::vector<std::vector<uint32_t>>& elems = body_ptr->GetElements();

	for (auto e = elems.begin(); e != elems.end(); ++e)
	{
		//Construct a local K and f matrix sized for the element (dof*#nodes)
		Matrix local_k(3*e->size(),3*e->size());
		Matrix local_f(3*e->size(),1);

		std::vector<uint32_t> local_elems = *e; //Store the set of element id's
		
		//Four nodes means a linear tetrahedral element (3D elements)
		if (e->size() == 4)
		{
			//Create a heap-allocated tetrahedral element & pointer to it which will be passed to 
			std::shared_ptr<TetrahedralElement> tet_ptr = std::make_shared<TetrahedralElement>(nodes, local_elems);

			//Update local stiffness and force vectors based on element and model information
			CalculateLocalK(local_k,tet_ptr);
			EnforceSurfaceBounds(local_k,local_f,tet_ptr);
			//CalculateLocalForce(local_f,tet_ptr);
		}
		if (e->size() == 8)
		{
			//Create a heap-allocated brick element & pointer to be passed to the calulation functions
			std::shared_ptr<BrickElement> brick_ptr = std::make_shared<BrickElement>(nodes,local_elems);

			CalculateLocalK(local_k, brick_ptr);
			EnforceSurfaceBounds(local_k, local_f, brick_ptr);
		}

		//Insert element stiffness and force into global system
		AssembleStiffness(local_k,*e);
		AssembleForce(local_f,*e);
	}
	EnforceDisplacements(global_k,global_f);
	//Prepare to solve the whole system
	LinearSystem system(global_k,global_f);
	//Update the  global solution vector using global force and stiffness matrices
	system.Solve(global_sol);
}

Matrix& LinearElasticSolids::GetElasticMatrix()
{
	return this->elastic_matrix;
}