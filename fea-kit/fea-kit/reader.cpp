#include "reader.h"

std::string& Reader::CheckComma(std::string& str)
{
	if (str.back() == ',')
	{
		str.pop_back();
	}
	return str;
}

void  Reader::ParseInstruction(std::fstream& f, const std::string& instruction, Body& body)
{
	if (instruction == "*NODE")
	{
		std::string input = "";
		while (f >> input && input != "*")
		{
			std::string id = input, x, y, z;
			f >> x;
			f >> y;
			f >> z;
			CheckComma(x);
			CheckComma(y);
			CheckComma(z);
			std::vector<double> n = { std::stod(x),std::stod(y),std::stod(z) };
			body.AddNode(n);
			input = "";
		}
		return;
	}
	else if (instruction == "*ELEMENT_SOLID")
	{
		std::string input = "";
		while (f >> input && input != "*")
		{
			if (input.front() == '$')
			{
				continue;
			}
			std::string id1 = input, id2, n1, n2, n3, n4, n5, n6, n7, n8;
			f >> id2 >> n1 >> n2 >> n3 >> n4 >> n5 >> n6 >> n7 >> n8;
			if (n5 == n4)
			{
				std::vector<uint32_t> el = {
					static_cast<uint32_t>(std::stoul(CheckComma(n1))),
					static_cast<uint32_t>(std::stoul(CheckComma(n2))),
					static_cast<uint32_t>(std::stoul(CheckComma(n3))),
					static_cast<uint32_t>(std::stoul(CheckComma(n4)))
				};
				body.AddElement(el);
			}
			else
			{
				std::vector<uint32_t> el = {
					static_cast<uint32_t>(std::stoul(CheckComma(n1))),
					static_cast<uint32_t>(std::stoul(CheckComma(n2))),
					static_cast<uint32_t>(std::stoul(CheckComma(n3))),
					static_cast<uint32_t>(std::stoul(CheckComma(n4))),
					static_cast<uint32_t>(std::stoul(CheckComma(n5))),
					static_cast<uint32_t>(std::stoul(CheckComma(n6))),
					static_cast<uint32_t>(std::stoul(CheckComma(n7))),
					static_cast<uint32_t>(std::stoul(CheckComma(n8)))
				};
				body.AddElement(el);
			}
		}
		return;
	}
	else if (instruction == "*SET_LINEAR_MATERIAL")
	{
		std::string sub_instruction = "";
		std::string val = "";

		while (f >> sub_instruction && sub_instruction != "*")
		{
			f >> val;
			body.SetMaterialProp(sub_instruction, std::stod(CheckComma(val)));
		}
		return;
	}
	else if (instruction == "*DISPLACEMENT" || instruction == "*TRACTION")
	{
		AddVectorBound(instruction, f, body);

		return;
	}
	else if (instruction == "*PRESSURE" || instruction == "*TEMPERATURE")
	{
		AddScalarBound(instruction, f, body);
		return;
	}
	else if (instruction == "*END")
	{
		f.close();
		return;
	}
}

void Reader::AddScalarBound(const std::string& instruction, std::fstream& f, Body& body)
{
	std::string input = "";
	std::string val;

	f >> val;
	std::vector<double> scal = {
		std::stod(CheckComma(val))
	};

	std::vector<uint32_t> node_list = {};

	while (f >> input && input != "*")
	{
		node_list.push_back(static_cast<uint32_t>(std::stoul(CheckComma(input))));
	}

	body.AddBoundary(node_list, "*DISPLACEMENT", scal);

	return;
}

void Reader::AddVectorBound(const std::string& inst, std::fstream& f, Body& body)
{
	std::string input = "";
	std::string u, v, w;

	f >> u >> v >> w;
	std::vector<double> vec = {
		std::stod(CheckComma(u)),
		std::stod(CheckComma(v)),
		std::stod(CheckComma(w))
	};

	std::vector<uint32_t> node_list = {};

	while (f >> input && input != "*")
	{
		node_list.push_back(static_cast<uint32_t>(std::stoul(CheckComma(input))));
	}

	body.AddBoundary(node_list, "*DISPLACEMENT", vec);

	return;
}


Reader::Reader() : instruction("")
{}

Reader::Reader(const std::string& file_name, Body& body) : instruction("")
{
	std::fstream file;
	file.open(file_name);
	while (file.is_open())
	{
		std::getline(file, instruction);
		ParseInstruction(file, instruction, body);
	}

}
