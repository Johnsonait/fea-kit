#pragma once
#include <vector>
#include <iostream>

#include "model.h"
#include "matrix.h"
#include "element.h"
#include "tetrahedral_element.h"
#include "body.h"
#include "reader.h"


//Class to solve problems in linear elasticity 
class LinearElasticSolids : public Model
{
private:
	double E, poisson, lambda, G; //Lame parameters.

	std::vector<std::vector<double>> global_k;

	Matrix elastic_matrix;

	void Lame();
	//Populates isotropic elasticity matrix with Lame parameters

	void ConstructElasticMatrix();

	//Constructs 6x12 elemental B matrix
	//Requires global shape function derivatives
	Matrix ConstructBMatrix(const double&,const double&,const double&,Element* el);

public:
	LinearElasticSolids();

	LinearElasticSolids(Reader& reader, Body& body);

	void Solve();
};