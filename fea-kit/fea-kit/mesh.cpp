#include "mesh.h"




struct mesh::bounds
{
    int FREE_SURFACE = 0;
    int FIXED_SURFACE = 1;
    int TRACTION_SURFACE = 2;
};

int boundary_types_number = 3;

mesh::mesh(std::string fileName) {
    std::vector<std::vector<double>> nodal;
    std::vector<std::vector<int>> elemental;
    std::vector<std::vector<int>> shal_list;
    std::vector<std::vector<int>> nodal_list;

    int boundary_types_number = 3;

    for (int i = 0; i < boundary_types_number; i++) { shal_list.push_back({}); nodal_list.push_back({}); }

    //Prepare to read mesh text file
    std::ifstream inFile;
    std::string x;
    inFile.open(fileName);
    if (!inFile)
    {
        std::cerr << "Unable to open file!" << std::endl;
        exit(1);
    }
    //"While values are being read"
    //Reading and storing mesh and instructions
    std::string instruction = "";
    while (inFile >> x)
    {
        int nodal_track = 0;
        int elemental_track = 0;

        if (x.front() == '*')
        {
            instruction = x;
        }
        //Read in node data
        if (instruction == "*NODE")
        {
            std::string n, N_x, N_y, N_z;
            if (x == instruction)
            {
                inFile >> n;
            }
            inFile >> N_x;
            inFile >> N_y;
            inFile >> N_z;

            nodal.push_back({ std::stod(N_x),std::stod(N_y),std::stod(N_z) });
        }

        //Read in beam (line) elements
        if (instruction == "*ELEMENT_BEAM")
        {
            readInNodalStructures(3, x, inFile, elemental);
        }
        // Read in 2D shell element declarations
        if (instruction == "*ELEMENT_SHELL")
        {
            readInNodalStructures(5, x, inFile, elemental);
        }
        //Read in solid element declarations
        if (instruction == "*ELEMENT_SOLID")
        {
            readInNodalStructures(9, x, inFile, elemental);
        }

        //Read and set boundary condition surfaces
        if (instruction == "*SET_SHELL_LIST")
        {
            bounds bound_track;
            
            int what_boundary = 0;
            std::string prelude, name, saver;
            inFile >> prelude;
            inFile >> name;

            if (name == "free_surface") { what_boundary = bound_track.FREE_SURFACE; }
            if (name == "fixed_surface") { what_boundary = bound_track.FIXED_SURFACE; }
            if (name == "traction_surface") { what_boundary = bound_track.TRACTION_SURFACE; }
            inFile >> saver;

            std::string rd;
            inFile >> rd;
            while (rd.front() != '*' && rd.front() != '$')
            {
                shal_list[what_boundary].push_back(std::stoi(rd));
                inFile >> rd;
            }
        }

        //Read and set nodes that have applied boundary conditions
        if (instruction == "*SET_NODE_LIST")
        {
            bounds bound_track;
            int what_boundary{};
            std::string prelude, name, saver;
            inFile >> prelude;
            inFile >> name;

            if (name == "free_surface") { what_boundary = bound_track.FREE_SURFACE; }
            if (name == "fixed_surface") { what_boundary = bound_track.FIXED_SURFACE; }
            if (name == "traction_surface") {what_boundary = bound_track.TRACTION_SURFACE; }
            inFile >> saver;

            std::string rd;
            inFile >> rd;
            while (rd.front() != '*' && rd.front() != '$')
            {
                nodal_list[what_boundary].push_back(std::stoi(rd));
                inFile >> rd;
            }
        }
    }
    inFile.close();

    //Set private variables to data that has been read
    nodes = nodal;
    elements = elemental;
    shell_list = shal_list;
}

void mesh::readInNodalStructures(int num, std::string x, std::ifstream& inFile, std::vector<std::vector<int>>& readArray)
{
    std::string prelude, n;
    std::vector<int> out;

    if (x.front() == '*')
    {
        inFile >> prelude;
        inFile >> n;
    }
    else if (x.front() == '$')
    {
        inFile >> n;
    }
    else
    {
        n = x;
    }
    for (int i = 0; i < num; i++)
    {
        std::string v;
        inFile >> v;
        out.push_back(std::stoi(v));
    }
    readArray.push_back(out);
}

void mesh::printNodes()
{
    for (int i = 1; i < nodes.size(); i++) {
        std::cout << "Node: " << i << std::endl;
        std::cout << "x: " << nodes[i][0] << " y: " << nodes[i][1] << " z: " << nodes[i][2] << std::endl;
    }
}

void mesh::printElements()
{
    for (int i = 0; i < elements.size(); i++) {
        std::cout << "Element:" << i << " ";
        for (int j = 0; j < elements[i].size(); j++) {
            std::cout << elements[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void mesh::printShellList()
{
    for (int i = 0; i < boundary_types_number; i++)
    {
        std::cout << "Boundary: " << i << std::endl;
        for (int j = 0; j < shell_list[i].size(); j++)
        {
            std::cout << shell_list[i][j] << " ";
        }
        std::cout << "" << std::endl;
    }

}