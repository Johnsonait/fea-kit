#include "body.h"

//Constructors
Body::Body()
{
	nodes = {};
	element_id = {};
	displacement = {};
	strain = {};
	stress = {};
	temperature = {};
	boundary_nodes = {};
	boundary_types = {};
	boundary_values = {};
}

//Accessors
const std::vector<std::vector<double>>& Body::GetNodes() { return nodes; }
const std::vector<std::vector<uint32_t>>& Body::GetElements() { return element_id; }
const std::vector<std::vector<double>>& Body::GetDisplacement() { return displacement; }
const std::vector<std::vector<double>>& Body::GetStrain() { return strain; }
const std::vector<std::vector<double>>& Body::GetStress() { return stress; }
const std::vector<double>& Body::GetTemperature() { return temperature; }
void Body::GetBoundaryInfo(std::vector<uint32_t>& b_n, std::string& s, std::vector<double>& b_v, const uint32_t& index)
{
	b_n = boundary_nodes[index];
	s = boundary_types[index];
	b_v = boundary_values[index];
}
const double& Body::GetStiffness() { return elastic_modulus; }
const double& Body::GetConductivity() { return conductivity; }
const double& Body::GetPoisson() { return poisson_ratio; }

uint32_t Body::GetNodeNum()
{
	return nodes.size();
}
uint32_t Body::GetDOF()
{
	return nodes[0].size();
}
uint32_t Body::GetElementCount()
{
	return element_id.size();
}

//Mutators
void Body::AddNode(const std::vector<double>& n)
{ 
	nodes.push_back(n);
}
void Body::AddElement(const std::vector<uint32_t>& e) 
{
	element_id.push_back(e);
}
void Body::AddBoundary(const std::vector<uint32_t>& b_n, const std::string& s, const std::vector<double>& b_v)
{
	boundary_nodes.push_back(b_n);
	boundary_types.push_back(s);
	boundary_values.push_back(b_v);
}

void Body::SetMaterialProp(const std::string& type, const double& val)
{
	if (type == "*CONDUCTIVITY")
	{
		conductivity = val;
	}
	else if (type == "*ELASTIC_MODULUS")
	{
		elastic_modulus = val;
	}
	else if (type == "*POISSON_RATIO")
	{
		poisson_ratio = val;
	}
}	
