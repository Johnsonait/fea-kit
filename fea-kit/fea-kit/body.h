#pragma once
#include <vector>
#include <string>
#include <algorithm>

#include "element.h"

class Body
{
private:
	std::vector<std::vector<double>> nodes;
	std::vector<std::vector<uint32_t>> element_id;
	std::vector<std::vector<double>> displacement;
	std::vector<std::vector<double>> strain;
	std::vector<std::vector<double>> stress;
	std::vector<double> temperature;
	std::vector<std::vector<uint32_t>> boundary_nodes; //Store nodes (id) for bounds
	std::vector<std::string> boundary_types; //Store bound types for each bound (eg disp...)
	std::vector<std::vector<double>> boundary_values; //Can be scalar, vector
	
	//Material properties
	double conductivity;
	double elastic_modulus;
	double poisson_ratio;	
	double density;

public:
	//Constructors
	Body();

	//Accessors
	const std::vector<std::vector<double>>& GetNodes();
	const std::vector<std::vector<uint32_t>>& GetElements();
	const std::vector<std::vector<double>>& GetDisplacement();
	const std::vector<std::vector<double>>& GetStrain();
	const std::vector<std::vector<double>>& GetStress();
	const std::vector<double>& GetTemperature();
	void GetBoundaryInfo(std::vector<uint32_t>& b_n, std::string& s, std::vector<double>& b_v, const uint32_t& index);
	uint32_t GetBoundaryCount();
	const double& GetStiffness();
	const double& GetPoisson();
	const double& GetConductivity();
	const double& GetDensity();
	uint32_t GetNodeNum();
	uint32_t GetDOF();
	uint32_t GetElementCount();

	void SearchBoundaryInfo(std::vector<uint32_t>& b_n, std::vector<std::string>& s, std::vector<std::vector<double>>& b_v, std::shared_ptr<Element> el_ptr);

	//Mutators
	void AddNode(const std::vector<double>& n);
	void AddElement(const std::vector<uint32_t>& e);
	void AddBoundary(const std::vector<uint32_t>& b_n, const std::string& s, const std::vector<double>& b_v);

	void SetMaterialProp(const std::string& type, const double& val);

};