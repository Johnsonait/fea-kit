#include "linear_elastic_solids.h"

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
Matrix LinearElasticSolids::ConstructBMatrix(const double& zeta,const double& eta, const double& mu, std::shared_ptr<Element> el)
{
	std::vector<std::vector<double>> nodes = el->GetNodes();
	el->CalcGlobalShapeDerivatives(zeta, eta, mu, nodes.size());
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

//This function is used to update the provided local_k matrix using quadrature integrator
void LinearElasticSolids::CalculateLocalK(Matrix& local_k, std::shared_ptr<Element> el_ptr)
{
	Quadrature integrator;
	local_k = integrator.Integrate(2, StiffnessIntegrand, local_k, el_ptr, this);
}

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
	InitMatrices(global_sol, sz, sz);//Create square matrix of zeros to be added to during solution;

	E = body_ptr->GetStiffness();
	poisson = body_ptr->GetPoisson();
	Lame();
	ConstructElasticMatrix();
}

void LinearElasticSolids::Solve()
{
	const std::vector<std::vector<double>>& nodes = body_ptr->GetNodes();
	const std::vector<std::vector<uint32_t>>& elems = body_ptr->GetElements();

	for (auto e = elems.begin(); e != elems.end(); ++e)
	{
		Matrix local_k(e->size(),e->size());
		Matrix local_f(e->size(),1);

		std::vector<uint32_t> local_elems = *e; //Store the set of element id's
		std::vector<std::vector<double>> local_nodes = {};
		for (auto local_e = local_elems.begin(); local_e != local_elems.end(); ++local_e)
		{
			local_nodes.push_back(nodes[(*local_e)-1]);//Note that element id's start at 1 so we need to subtract one to access the proper index
		}
		
		if (e->size() == 4)
		{
			//Create a heap-allocated tetrahedral element & pointer to it which will be passed to 
			std::shared_ptr<TetrahedralElement> tet_ptr = std::make_shared<TetrahedralElement>(local_nodes, local_elems);
			//TODO
			//Update local stiffness and force vectors based on element and model information
			CalculateLocalK(local_k,tet_ptr);
			CalculateLocalForce(local_f,tet_ptr);
		}
		//TODO
		//
		AssembleStiffness(local_k,*e);
		AssembleForce(local_f,*e);
	}
}

Matrix& LinearElasticSolids::GetElasticMatrix()
{
	return this->elastic_matrix;
}

Matrix& StiffnessIntegrand(double eta, double zeta, double mu, Matrix mat,std::shared_ptr<Element> el_ptr,LinearElasticSolids* model)
{
	const Matrix B = model->ConstructBMatrix(eta,zeta,mu, el_ptr);
	const Matrix B_T = B.GetTranspose();
	const Matrix C = model->GetElasticMatrix();
	Matrix ret = B_T * C * B * el_ptr->GetJacobianDet();
	return ret;
}