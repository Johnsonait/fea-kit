#pragma once
#include <vector>
#include <iostream>
#include <memory>

#include "model.h"
#include "matrix.h"
#include "element.h"
#include "tetrahedral_element.h"
#include "body.h"
#include "reader.h"
#include "quadrature.h"


//Class to solve problems in linear elasticity 
class LinearElasticSolids : public Model
{
private:
	double E, poisson, lambda, G; //Lame parameters.
	Body* body_ptr;
	Reader* read_ptr;
	std::shared_ptr<std::vector<std::vector<double>>> global_k;
	std::shared_ptr<std::vector<std::vector<double>>> global_f;
	std::shared_ptr<std::vector<std::vector<double>>> global_sol;
	Matrix elastic_matrix;

	void Lame();
	//Populates isotropic elasticity matrix with Lame parameters
	void ConstructElasticMatrix();

	void InitMatrices(std::shared_ptr<std::vector<std::vector<double>>>,const uint32_t&,const uint32_t&);

	void CalculateLocalK(Matrix&, std::shared_ptr<Element> );
	void CalculateLocalForce(Matrix&,std::shared_ptr<Element>);

	void AssembleStiffness(Matrix&,const std::vector<uint32_t>&);
	void AssembleForce(Matrix&, const std::vector<uint32_t>&);

public:
	LinearElasticSolids();

	LinearElasticSolids(Reader& reader, Body& body);

	void Solve(); 
	Matrix ConstructBMatrix(const double&, const double&, const double&, std::shared_ptr<Element>);
	Matrix& GetElasticMatrix();
};