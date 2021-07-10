#pragma once

class Body
{
public:
	std::vector<std::vector<double>> nodes;
	std::vector<std::vector<int>> traction_boundary_nodes;
	std::vector<std::vector<int>> displacement_boundary_nodes;

	std::vector<std::vector<double>> displacement;
	std::vector<std::vector<double>> stress;
	std::vector<std::vector<double>> elastic_matrix;
	std::vector<double> temperature;

	void GetStrain();

	void GetStress();

	void OutputData();
};