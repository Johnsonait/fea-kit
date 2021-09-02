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
const std::vector<double>& Body::GetEquivalentStrain() { return equiv_strain; }
const std::vector<std::vector<double>>& Body::GetStress() { return stress; }
const std::vector<double>& Body::GetEquivalentStress() { return equiv_stress; }
const std::vector<double>& Body::GetTemperature() { return temperature; }

//Function that takes in references to a set of boundary nodes, a string storing boundary type, a string storing the values associated with the boundary
//and an index. The index identifies which boundary is the one you want to know about
void Body::GetBoundaryInfo(std::vector<uint32_t>& b_n, std::string& s, std::vector<double>& b_v, const uint32_t& index)
{
	b_n = boundary_nodes[index];
	s = boundary_types[index];
	b_v = boundary_values[index];
}

void Body::SearchBoundaryInfo(std::vector<uint32_t>& b_n, std::vector<std::string>& s, std::vector<std::vector<double>>& b_v, std::shared_ptr<Element> el_ptr)
{
	std::vector<uint32_t> global_ids = el_ptr->GetGlobalIDs();

	for (size_t n = 0; n < global_ids.size(); ++n)
	{
		for (size_t bound = 0; bound < boundary_nodes.size(); ++bound)
		{
			//Search the current boundary to see if the current node is contained in it
			//If true, update the provided vectors so that the caller knows which nodes have which boundaries
			if (std::binary_search(boundary_nodes[bound].begin(),boundary_nodes[bound].end(),global_ids[n]))
			{
				b_n.push_back(global_ids[n]);
				s.push_back(boundary_types[bound]);
				b_v.push_back(boundary_values[bound]);
			}
		}
	}
}

uint32_t Body::GetBoundaryCount()
{
	return boundary_nodes.size();
}

const double& Body::GetStiffness() { return elastic_modulus; }
const double& Body::GetConductivity() { return conductivity; }
const double& Body::GetPoisson() { return poisson_ratio; }
const double& Body::GetDensity() { return density; }


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
	else if (type == "*DENSITY")
	{
		density = val;
	}
}	

void Body::AddDisplacement(const std::vector<double>& vector)
{
	displacement.push_back(vector);
}

void Body::AddStrain(const std::vector<double>& vector)
{
	strain.push_back(vector);
}

void Body::AddEquivalentStrain(const double& value)
{
	equiv_strain.push_back(value);
}

void Body::AddStress(const std::vector<double>& vector)
{
	stress.push_back(vector);
}

void Body::AddEquivalentStress(const double& value)
{
	equiv_stress.push_back(value);
}

void Body::AddTemperature(const double& value)
{
	temperature.push_back(value);
}