#include "writer.h"

Writer::Writer()
{
	file_name = "";
	body_ptr = nullptr;
}

Writer::Writer(const std::string& file, Body& body)
{
	file_name = file;
	body_ptr = &body;
}

void Writer::Write()
{
	//Create file to write into
	std::ofstream outfile(file_name);
	//We are writing into the file such that each node has a line dedicated to each node. On this node's line is every property associated with the node 
	//like position (X Y Z) displacement (U V W) strain (e_xx e_yy e_zz e_yz e_xz e_xy) equivalent strain (Eeq) stress (Sxx Syy Szz Syz Sxz Sxy) equivalent stress (Seq) and temperature (T)
	
	//outfile << "(X Y Z) (U V W) (e_xx e_yy e_zz e_yz e_xz e_xy) E_eq (Sxx Syy Szz Syz Sxz Sxy) Seq T " << "\n";
	for (size_t node = 0; node < body_ptr->GetNodes().size(); ++node)
	{
		//Write in positions
		for (size_t i = 0; i < body_ptr->GetNodes()[node].size();++i)
		{
			outfile << body_ptr->GetNodes()[node][i] << ", ";
		}
		//Write in displacements
		for (size_t i = 0; i < body_ptr->GetDisplacement()[node].size(); i++)
		{
			outfile << body_ptr->GetDisplacement()[node][i] << ", ";
		}
		//Write in strains
		for (size_t i = 0; i < body_ptr->GetStrain()[node].size(); i++)
		{
			outfile << body_ptr->GetStrain()[node][i] << ", ";
		}
		//Write in equivalent strain
		outfile << body_ptr->GetEquivalentStrain()[node] << ", ";
		//Write in stresses
		for (size_t i = 0; i < body_ptr->GetStress()[node].size(); i++)
		{
			outfile << body_ptr->GetStress()[node][i] << ", ";
		}
		//Write in equivalent stress
		outfile << body_ptr->GetEquivalentStress()[node] << ", ";
		//Write in the temperature
		outfile << body_ptr->GetTemperature()[node] << " ";
		//Write in endline, all nodal data has been inputed
		outfile << "\n";
	}
	//Once all nodes have been stored, close the file
	outfile.close();
}