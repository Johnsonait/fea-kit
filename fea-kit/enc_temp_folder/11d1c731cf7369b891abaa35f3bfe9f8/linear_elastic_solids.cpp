#include "linear_elastic_solids.h"

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

static const double GRAVITY = -9.81;

Matrix StiffnessIntegrand(double xsi, double eta, double zeta, std::shared_ptr<Element> el_ptr, LinearElasticSolids* model)
{
	const Matrix B = model->ConstructBMatrix(xsi, eta, zeta, el_ptr);
	const Matrix B_T = B.GetTranspose();
	const Matrix C = model->GetElasticMatrix();

	//Final integrand value
	Matrix ret = (B_T * (C * B)) * el_ptr->GetJacobianDet(xsi,eta,zeta);
	return ret;
}

//Finite element form of integrand used for updating local force vector to include body forces
Matrix& BodyForceIntegrand(double xsi, double eta, double zeta, std::shared_ptr<Element> el_ptr, LinearElasticSolids* model)
{
	Matrix N = model->ConstructShapeMatrix(xsi,eta,zeta,el_ptr);
	std::vector<std::vector<double>> temp = { {0},{0},{GRAVITY} };
	Matrix body_force(temp);

	Matrix ret = (N*body_force)*el_ptr->GetJacobianDet(xsi,eta,zeta);
	return ret;
}

Matrix& SurfaceForceIntegrand(double xsi, double eta, std::shared_ptr<Element> el_ptr, LinearElasticSolids* model,std::vector<std::vector<double>> traction)
{

	Matrix N = model->ConstructShapeMatrix(xsi, eta,0, el_ptr);
	const std::vector<std::vector<double>> J = el_ptr->GetJacobian(xsi, eta, 0);

	Matrix tract(traction);
	
	Matrix ret = N * (tract.Transpose()) * std::sqrt((std::pow(((J[1][0] * J[2][1]) - (J[1][1] * J[2][0])), 2)) + (std::pow(((J[0][0] * J[2][1]) - (J[0][1] * J[2][0])), 2)) + (std::pow(((J[0][0] * J[1][1]) - (J[0][1] * J[1][0])), 2)));
	return ret;
}

void LinearElasticSolids::Log(const std::string& str)
{
	std::cout << "\n" << str << std::endl;
}

void LinearElasticSolids::Lame() //Pronounced Lam-eh
{
	lambda = (E * poisson) / ((1 + poisson) * (1 - (2 * poisson)));
	G = E / (2 * (1 + poisson));
}
double LinearElasticSolids::EquivalentStrain(const std::vector<double>& strain)
{
	return (2.0 / 3.0) * std::sqrt(((3.0 / 2.0) * (std::pow(strain[0], 2) + std::pow(strain[1], 2) + std::pow(strain[2], 2))) + ((3.0 / 4.0) * (std::pow(strain[3], 2) + std::pow(strain[4], 2) + std::pow(strain[5], 2))));
}
double LinearElasticSolids::EquivalentStress(const std::vector<double>& stress)
{
	return std::sqrt(0.5 * (std::pow((stress[0] - stress[1]), 2) + std::pow((stress[1] - stress[2]), 2), std::pow((stress[2] - stress[0]), 2) + (6 * (std::pow(stress[3], 2) + std::pow(stress[4], 2) + std::pow(stress[5], 2)))));
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
	local_k = Integrate(3, StiffnessIntegrand,local_k,el_ptr,this);
}

//This function updates the global stiffness matrix from the element stiffness matrix
void LinearElasticSolids::AssembleStiffness(Matrix& local_k,const std::vector<uint32_t>& node_ids)
{
	for (size_t i = 0; i < node_ids.size(); ++i)
	{
		for (size_t j = 0; j < node_ids.size(); ++j)
		{
			for (size_t n = 0; n < 3; ++n)
			{
				(*global_k)[(3*node_ids[i]-2) - 1 + n][(3*node_ids[j]-2) - 1 + n] += local_k[i + n][j + n];
			}
		}
	}
}
void LinearElasticSolids::AssembleForce(Matrix& local_f, const std::vector<uint32_t>& node_ids)
{
	for (size_t i = 0; i < node_ids.size(); ++i)
	{
		for (size_t n = 0; n < 3; ++n)
		{
			(*global_f)[3*node_ids[i]-2 - 1 + n][0] += local_f[i + n][0];
		}
	}
}

void LinearElasticSolids::EnforceSurfaceBounds(Matrix& local_k, std::shared_ptr<Element> el_ptr)
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
				if (is_in)
				{
					affected_bound.push_back(bound);
				}
			}
			//Now we have the bounds of the element that are affected by the traction vector 
			//(Which is stored in boundary_vector)
			//The next step is to integrate the bounds surfaces in affected_bound with the vector and shape matrices, then update the local_f
			//For now forces will be calculated assuming the traction is constant across the surface and averaging the force across all nodes
			for (auto bound : affected_bound)
			{
				//Bound stores id's of nodes in the associated bound. Can use ids to get node coordinates. Can use node coordinates (and cross product) to estimate area of surface
				//Area of surface can then be multiplied by traction value to get average force Fx, Fy, Fz. Number of nodes in bound can be used to find average force for each node
				//Average force is then added to appropriate local_f

				std::vector<double> node1 = (body_ptr->GetNodes())[bound[0]-1];
				std::vector<double> node2 = (body_ptr->GetNodes())[bound[1]-1];
				std::vector<double> node3 = (body_ptr->GetNodes())[bound[2]-1];
				std::vector<double> node4 = (body_ptr->GetNodes())[bound[3]-1];

				uint32_t x = 0, y = 1, z = 2;

				double A1 = 0.5 * std::sqrt((pow((((node2[y] - node1[y]) * (node4[z] - node1[z])) - ((node4[y] - node1[y]) * (node2[z] - node1[z]))), 2)) + (pow((((node2[x] - node1[x]) * (node4[z] - node1[z])) - ((node4[x] - node1[x]) * (node2[z] - node1[z]))), 2)) + (pow((((node2[x] - node1[x]) * (node4[y] - node1[y])) - ((node4[x] - node1[x]) * (node2[y] - node1[y]))), 2)));
				double A2 = 0.5 * std::sqrt((pow((((node4[y] - node3[y]) * (node2[z] - node3[z])) - ((node2[y] - node3[y]) * (node4[z] - node3[z]))), 2)) + (pow((((node4[x] - node3[x]) * (node2[z] - node3[z])) - ((node2[x] - node3[x]) * (node4[z] - node3[z]))), 2)) + (pow((((node4[x] - node3[x]) * (node2[y] - node3[y])) - ((node2[x] - node3[x]) * (node4[y] - node3[y]))), 2)));
				
				for (auto node : bound)
				{
					uint32_t index = 3 * node - 2;
					for (size_t n = 0; n < boundary_vector.size(); ++n)
					{
						(*global_f)[index - 1 + n][0] += boundary_vector[n] * (A1 + A2)/ bound.size();
					}
				}
			}
		}
	}
}
void LinearElasticSolids::GenerateSymmetry(std::shared_ptr<std::vector<std::vector<double>>> k)
{
	for (size_t i = 0; i < k->size(); i++)
	{
		for (size_t j = 0; j < (*k)[i].size(); j++)
		{
			(*k)[j][i] = (*k)[i][j];
		}
	}
}

//This function is used to update the provided local_f matrix using quadrature 
//It needs to consider the volume integral over the body and the surface integral over the tractioned surfaces
void LinearElasticSolids::CalculateLocalForce(Matrix& local_f, std::shared_ptr<Element> el_ptr,std::vector<std::vector<double>>& traction)
{
	CalculateBodyForce(local_f, el_ptr,this);
	CalculateSurfaceForce(local_f, el_ptr,this,traction);
}

//Integrate the body-force values to be used in global force-vector
void LinearElasticSolids::CalculateBodyForce(Matrix& local_f, std::shared_ptr<Element> el_ptr,LinearElasticSolids* model)
{
	local_f = local_f + Integrate(2, BodyForceIntegrand,local_f, el_ptr,model);
}

//Integrate the surface tractions to be used in global force-vector
void LinearElasticSolids::CalculateSurfaceForce(Matrix& local_f, std::shared_ptr<Element> el_ptr,LinearElasticSolids* model,std::vector<std::vector<double>>& traction)
{
	local_f = local_f + IntegrateSurf(2,SurfaceForceIntegrand,local_f, el_ptr,model,traction);
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

//Changes global stiffness matrix entries and force vector
void LinearElasticSolids::EnforceDisplacements(std::shared_ptr<std::vector<std::vector<double>>> k, std::shared_ptr<std::vector<std::vector<double>>> f)
{
	//Go through all boundaries associated with body and apply displacements if there are any
	for (size_t boundary = 0; boundary < body_ptr->GetBoundaryCount(); ++boundary)
	{
		std::string type = "";
		std::vector<uint32_t> boundary_nodes = {};
		std::vector<double> boundary_vector = {};
		body_ptr->GetBoundaryInfo(boundary_nodes, type, boundary_vector, boundary);
		if (type == "*DISPLACEMENT")
		{
			//Go through global ids of nodes in boundary_node and enforce displacement
			for (uint32_t node : boundary_nodes)
			{
				uint32_t k_index = 3 * node - 3;
				//What we do here applies to each dof (there are 3)
				for (size_t dof = 0; dof < 3; dof++)
				{
					//clear the row of all values
					//for (size_t col = 0; col < (*k)[node-1+dof].size(); col++)
					//{
					//	(*k)[k_index+ dof][col] = 0;
					//}
					//Now set the appropriate entry of stiffness to 1
					(*k)[k_index + dof][k_index + dof] = 10e20;
					(*f)[k_index+dof][0] = boundary_vector[dof];
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
		
		if (e->size() == 8)
		{
			//Create a heap-allocated brick element & pointer to be passed to the calulation functions
			std::shared_ptr<BrickElement> brick_ptr = std::make_shared<BrickElement>(nodes,local_elems);

			CalculateLocalK(local_k, brick_ptr);
			EnforceSurfaceBounds(local_k, brick_ptr);
		}
		//Insert element stiffnes into global system
		AssembleStiffness(local_k,*e);
		//AssembleForce(local_f,*e);
	}
	Log("Enforcing displacement conditions");
	EnforceDisplacements(global_k,global_f);
	//Prepare to solve the whole system
	LinearSystem system(global_k,global_f);
	//Update the  global solution vector using global force and stiffness matrices
	Log("Starting system solution...");
	system.Solve(global_sol);
	Log("System solution complete!");
	//*
	std::cout << std::left;
	for (size_t i = 0; i < global_k->size(); i++)
	{
		for (size_t j = 0; j < global_k->size(); j++)
		{
			double out = (*global_k)[i][j];
			if (out != 0)
			{
				out = 1;
				std::cout << std::setw(3) << std::setprecision(2) << "\033[31m";
				std::cout << out << " ";
			}
			else
			{
				
				std::cout << std::setw(3) << std::setprecision(2)<<"\033[0m" << out << " ";
			}
		}
		std::cout << std::endl;
	}
	std::cout << "\033[0m";
	//*/
	for (size_t i = 0; i < global_f->size(); i++)
	{
		std::cout << (*global_sol)[i][0] << std::endl;
	}
	//Solution finished! Time to post-process
}

//method generates stress and strain data from global_sol and put data into body_ptr to be writen to results file
void LinearElasticSolids::PostProcess()
{
	Log("Starting post-processing...");
	//Starting out, we have the displacements U V W stored for each node in the global_sol pointer, lets add them to the body_ptr via AddDisplacement
	for (size_t node = 0; node < global_sol->size(); node+=3)
	{
		std::vector<double> temp;
		for (size_t n = 0; n < 3; n++)
		{
			temp.push_back((*global_sol)[node+n][0]);
		}
		body_ptr->AddDisplacement(temp);
	}
	//Strains can be computed using the B matrix approach (e = B*u)
	//The body_ptr contains the strains vector we want to update, we can go through each element (like we do in Solve()) and additively construct the strains from it
	//Get references to the body nodes and elements for convienence
	const std::vector<std::vector<double>>& nodes = body_ptr->GetNodes();
	const std::vector<std::vector<uint32_t>>& elems = body_ptr->GetElements();
	const std::vector<std::vector<double>> BRICK_POINTS = {
		{-1,-1,-1},
		{1,-1,-1},
		{1,1,-1},
		{-1,1,-1},
		{-1,-1,1},
		{1,-1,1},
		{1,1,1},
		{-1,1,1}
	};

	std::vector<std::vector<double>> global_strain;
	std::vector<std::vector<double>> global_stress;
	std::vector<double> global_temp;
	//Prepare array for size of total nodes in the body
	for (size_t n = 0; n < body_ptr->GetNodeNum(); n++)
	{
		std::vector<double> temp(6,0);
		global_strain.push_back(temp);
		global_stress.push_back(temp);
		global_temp.push_back(0);
	}
	//Iterate through each element and add each elements contribution to the results
	for (auto e = elems.begin(); e != elems.end(); ++e)
	{
		std::vector<uint32_t> local_elems = *e; //Store the set of element id's
		//If element size is 8, use brick_element formulations
		if (e->size() == 8)
		{
			//Create a heap-allocated brick element & pointer to be passed to the calulation functions
			std::shared_ptr<BrickElement> brick_ptr = std::make_shared<BrickElement>(nodes, local_elems);
			//Go through each node of the element and calculate the strain at that node via B*u
			//Remember we need to extract the proper "local" u from the global solution vector
			std::vector<std::vector<double>> temp_u;
			//Iterate through all the nodes in the element
			for (size_t n = 0; n < local_elems.size(); n++)
			{
				//3 dof to take into account
				for (size_t i = 0; i < 3; ++i)
				{
					temp_u.push_back({});
					double val = (*global_sol)[(3*local_elems[n]-2)-1+i][0];
					temp_u[3*n + i].push_back(val);
				}
			}
			Matrix local_u(temp_u);
			//Iterate throgu all nodes in the element
			for (size_t n = 0; n < local_elems.size(); ++n)
			{
				double xsi = BRICK_POINTS[n][0], eta = BRICK_POINTS[n][1], zeta = BRICK_POINTS[n][2];
				Matrix B = ConstructBMatrix(xsi,eta,zeta, brick_ptr);
				Matrix strain_mat = B * local_u;
				Matrix stress_mat = elastic_matrix * strain_mat;
				//We now have the proper values of strain and stress in strain_mat and stress_mat at the node n
				//Time to add it to the global stress vector!
				for (size_t i = 0; i < strain_mat.CountRows(); ++i)
				{
					global_strain[(local_elems[n]) - 2 + i][0] += strain_mat[i][0];
					global_stress[(local_elems[n]) - 2 + i][0] += stress_mat[i][0];
				}
			}

		}
		//end of if for element types
		//Now add the calculated results to the body_ptr
		for (size_t i = 0; i < body_ptr->GetNodeNum(); i++)
		{
			body_ptr->AddStrain(global_strain[i]);
			body_ptr->AddStress(global_stress[i]);
			body_ptr->AddEquivalentStrain(EquivalentStrain(global_strain[i]));
			body_ptr->AddEquivalentStress(EquivalentStress(global_stress[i]));
			body_ptr->AddTemperature(global_temp[i]);
		}
	}
	Log("Data processing complete!");
	
}

Matrix& LinearElasticSolids::GetElasticMatrix()
{
	return this->elastic_matrix;
}

//3D integrals
Matrix LinearElasticSolids::Integrate(
	const int& points, 
	std::function<Matrix (double, double, double, std::shared_ptr<Element>,LinearElasticSolids*)> func,
	Matrix& mat, 
	std::shared_ptr<Element> el_ptr, 
	LinearElasticSolids* model)
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
				result = result + (func(BRICK_QUADRATURE_POINTS[index][i], BRICK_QUADRATURE_POINTS[index][j], BRICK_QUADRATURE_POINTS[index][k], el_ptr,model)
					* (BRICK_QUADRATURE_WEIGHTS[index][i] * BRICK_QUADRATURE_WEIGHTS[index][j] * BRICK_QUADRATURE_WEIGHTS[index][k]));
			}
		}
	}
	return result;
}

//2D integrals
Matrix& LinearElasticSolids::IntegrateSurf(const int& points, std::function<Matrix& (double, double, std::shared_ptr<Element>,LinearElasticSolids*,std::vector<std::vector<double>>&)> func, Matrix& mat, std::shared_ptr<Element> el_ptr,LinearElasticSolids* model,std::vector<std::vector<double>>& traction)
{
	int index = points - 1;
	//Construct result matrix of the appropriate size
	Matrix result(mat.CountRows(), mat.CountCols());
	for (int i = 0; i < points; i++)
	{
		for (int j = 0; j < points; j++)
		{
			result = result + (func(BRICK_QUADRATURE_POINTS[index][i], BRICK_QUADRATURE_POINTS[index][j],el_ptr,model,traction))
				* (BRICK_QUADRATURE_WEIGHTS[index][i] * BRICK_QUADRATURE_WEIGHTS[index][j]);
		}
	}
	return result;
}