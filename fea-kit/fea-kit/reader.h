#pragma once

#include <string>
#include <vector>
#include <fstream>

#include "body.h"
//Class used to read and parse model data, parameters, and boundary conditions
class Reader
{
private:
	std::string instruction;

	std::string& CheckComma(std::string& str);

	void ParseInstruction(std::fstream& f, const std::string& instruction, Body& body);

	void AddScalarBound(const std::string& instruction, std::fstream& f, Body& body);

	void AddVectorBound(const std::string& inst, std::fstream& f, Body& body);

public:
	Reader();
	Reader(const std::string& file_name, Body& body);

};