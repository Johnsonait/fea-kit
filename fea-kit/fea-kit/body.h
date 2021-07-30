#pragma once

class Body
{
private:
	std::vector<std::vector<double>> nodes;
	std::vector<std::vector<uint32_t>> elements;
	std::vector<std::vector<double>> displacement;
	std::vector<std::vector<double>> strain;
	std::vector<std::vector<double>> stress;
	std::vector<double> temperature;
	std::vector<std::vector<uint32_t>> boundary_nodes; //Store nodes (id) for bounds
	std::vector<std::string> boundary_types; //Store bound types for each bound (eg disp...)
	std::vector<std::vector<double>> boundary_values; //Can be scalar, vector

public:
	//Constructors
	Body()
	{
		nodes = {};
		elements = {};
		displacement = {};
		strain = {};
		stress = {};
		temperature = {};
		boundary_nodes = {};
		boundary_types = {};
		boundary_values = {};
	}

	//Accessors
	const std::vector<std::vector<double>>& GetNodes();
	const std::vector<std::vector<uint32_t>>& GetElements();
	const std::vector<std::vector<double>>& GetDisplacement();
	const std::vector<std::vector<double>>& GetStrain();
	const std::vector<std::vector<double>>& GetStress();
	const std::vector<double>& GetTemperature();
	void GetBoundaryInfo(std::vector<uint32_t>& b_n, std::string& s, std::vector<double>& b_v, const uint32_t& index);

	//Mutators
	void AddNode(const std::vector<double>& n);
	void AddElement(const std::vector<uint32_t>& e);
	void AddBoundary(const std::vector<uint32_t>& b_n, const std::string& s, const std::vector<double>& b_v);
};