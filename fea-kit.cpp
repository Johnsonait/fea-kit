#include <iostream>
#include <vector>
#include <fstream>
#include <string>

//using namespace std;

class mesh{
    private:
    std::vector<std::vector<double>> nodes;
    std::vector<std::vector<int>> elements;
    std::vector<std::vector<int>> boundaries;

    public:
    
    enum bounds
    {
        FREE_SURFACE,
        FIXED_SURFACE,
        TRACTION_SURFACE
    };

    mesh(std::string fileName){
        std::vector<std::vector<double>> nodal;
        std::vector<std::vector<int>> elemental;
        std::vector<std::vector<int>> boundal;
        
        //Prepare to read mesh text file
        std::ifstream inFile;
        std::string x;
        inFile.open(fileName);
        if(!inFile){
            std::cerr<<"Unable to open file!"<<std::endl;
            exit(1);
        }
        //"While values are being read"
        //Reading and storing mesh and instructions
        std::string instruction = "";
        while(inFile >> x)
        {
            int nodal_track = 0;
            int elemental_track = 0;

            if(x.front() =='*')
            {
                instruction = x;
            }
            //Read in node data
            if(instruction == "*NODE")
            {
                std::string n,N_x,N_y,N_z;
                if(x == instruction)
                {
                    inFile >> n;
                }
                inFile >> N_x;
                inFile >> N_y;
                inFile >> N_z;        

                nodal.push_back({std::stod(N_x),std::stod(N_y),std::stod(N_z)});
            }
            //Read in beam (line) elements
            if(instruction == "*ELEMENT_BEAM")
            {
                std::string prelude, n, e_Type, node_1, node_2;
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
                inFile >> e_Type;
                inFile >> node_1;
                inFile >> node_2;
                elemental.push_back({std::stoi(e_Type),std::stoi(node_1),std::stoi(node_2)});
            }
            // Read in 2D shell element declarations
            if(instruction == "*ELEMENT_SHELL")
            {
                std::string prelude, n, e_Type, node_1, node_2, node_3, node_4;
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
                for (int i = 0; i < 5; i++)
                {
                    std::string v;
                    inFile >> v;
                    out.push_back(std::stoi(v));
                }
                elemental.push_back(out);              
            }
            //Read in solid element declarations
            if (instruction == "*ELEMENT_SOLID")
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
                for (int i = 0; i < 9; i++)
                {
                    std::string v;
                    inFile >> v;
                    out.push_back(std::stoi(v));
                }
                elemental.push_back(out);                     
            }
            //Read and set boundary condition surfaces
            if(instruction == "*SET_SHELL_LIST")
            {
                bounds bound_track;
                std::string prelude,name;
                inFile >> prelude;
                inFile >> name;
                
                if(name == "free_surface"){bound_track = FREE_SURFACE;}
                if(name == "fixed_surface"){bound_track = FIXED_SURFACE;}
                if(name == "traction_surface"){bound_track = TRACTION_SURFACE;}

            }
        }
        inFile.close();

        //Set private variables to data that has been read
        nodes = nodal;
        elements = elemental;
    }
    void printNodes(){
        for(int i =1;i<nodes.size();i++){
            std::cout<<"Node: "<<i<<std::endl;
            std::cout<<"x: "<<nodes[i][0]<<" y: "<<nodes[i][1]<<" z: "<<nodes[i][2]<<std::endl;
        }
    }
    void printElements(){
        for(int i=0; i<elements.size();i++){
            std::cout<<"Element:"<<i<<" ";
            for(int j = 0;j<elements[i].size();j++){
                std::cout<<elements[i][j]<< " ";
            }
            std::cout<<std::endl;
        }
    }
};

int main(){
    mesh MESH("mesh.txt");

    //MESH.printNodes();
    //MESH.printElements();

    return 0;
}