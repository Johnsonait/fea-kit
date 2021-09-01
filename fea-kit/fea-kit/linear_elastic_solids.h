#pragma once
#include <vector>
#include <iostream>
#include <memory>
#include <math.h>
#include <functional>

#include "model.h"
#include "matrix.h"
#include "element.h"
#include "brick_element.h"
#include "body.h"
#include "reader.h"
//#include "quadrature.h"

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
	void CalculateLocalForce(Matrix&,std::shared_ptr<Element>,std::vector<std::vector<double>>&);
	void CalculateBodyForce(Matrix& local_f, std::shared_ptr<Element> el_ptr,LinearElasticSolids*);
	void CalculateSurfaceForce(Matrix& local_f, std::shared_ptr<Element> el_ptr,LinearElasticSolids*,std::vector<std::vector<double>>&);

	void AssembleStiffness(Matrix&,const std::vector<uint32_t>&);
	void AssembleForce(Matrix&, const std::vector<uint32_t>&);
	void EnforceSurfaceBounds(Matrix& local_k, std::shared_ptr<Element> el_ptr);
	void EnforceDisplacements(std::shared_ptr<std::vector<std::vector<double>>> k, std::shared_ptr<std::vector<std::vector<double>>> f);

	Matrix Integrate(const int& points, std::function<Matrix (double, double, double, std::shared_ptr<Element>,LinearElasticSolids*)> func, Matrix& mat, std::shared_ptr<Element> el_ptr,LinearElasticSolids*);
	Matrix& IntegrateSurf(const int& points, std::function<Matrix& (double, double,std::shared_ptr<Element>, LinearElasticSolids*,std::vector<std::vector<double>>&)> func, Matrix& mat, std::shared_ptr<Element> el_ptr,LinearElasticSolids*,std::vector<std::vector<double>>&);

	void Log(const std::string&);

public:
	LinearElasticSolids();

	LinearElasticSolids(Reader& reader, Body& body);

	void Solve(); 
	static Matrix ConstructBMatrix(const double&, const double&, const double&, std::shared_ptr<Element>);
	static Matrix ConstructShapeMatrix(const double& zeta, const double& eta, const double& mu, std::shared_ptr<Element> el);
	Matrix& GetElasticMatrix();

	Body& GetBody();

};
