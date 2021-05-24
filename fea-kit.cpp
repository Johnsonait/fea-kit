#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>


class mesh{
    private:
    std::vector<std::vector<double>> nodes;
    std::vector<std::vector<int>> elements;
    std::vector<std::vector<int>> shell_list;
    std::vector<std::vector<int>> node_list;

    int boundary_types_number = 3;

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
        std::vector<std::vector<int>> shal_list;
        std::vector<std::vector<int>> nodal_list;

        for (int i = 0; i <= boundary_types_number; i++){shal_list.push_back({}); nodal_list.push_back({});}
        
        //Prepare to read mesh text file
        std::ifstream inFile;
        std::string x;
        inFile.open(fileName);
        if(!inFile)
        {
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
                readInNodalStructures(3, x , inFile, elemental);
            }
            // Read in 2D shell element declarations
            if(instruction == "*ELEMENT_SHELL")
            {
                readInNodalStructures(5, x, inFile, elemental);
            }
            //Read in solid element declarations
            if (instruction == "*ELEMENT_SOLID")
            {
                readInNodalStructures(9, x, inFile, elemental);
            }

            //Read and set boundary condition surfaces
            if(instruction == "*SET_SHELL_LIST")
            {
                bounds bound_track;
                std::string prelude,name,saver;
                inFile >> prelude;
                inFile >> name;
                
                if(name == "free_surface"){bound_track = FREE_SURFACE;}
                if(name == "fixed_surface"){bound_track = FIXED_SURFACE;}
                if(name == "traction_surface"){bound_track = TRACTION_SURFACE;}
                inFile >> saver;
                
                std::string rd;
                inFile >> rd;
                while(rd.front() != '*' && rd.front() != '$')
                {
                    shal_list[bound_track].push_back(std::stoi(rd));
                    inFile >> rd;
                }                
            }

            //Read and set nodes that have applied boundary conditions
            if (instruction == "*SET_NODE_LIST")
            {
                bounds bound_track;
                std::string prelude,name,saver;
                inFile >> prelude;
                inFile >> name;
                
                if(name == "free_surface"){bound_track = FREE_SURFACE;}
                if(name == "fixed_surface"){bound_track = FIXED_SURFACE;}
                if(name == "traction_surface"){bound_track = TRACTION_SURFACE;}   
                inFile >> saver;

                std::string rd;
                inFile >> rd;
                while (rd.front() != '*' && rd.front() != '$')
                {
                    nodal_list[bound_track].push_back(std::stoi(rd));
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

    void readInNodalStructures(int num, std::string x, std::ifstream &inFile, std::vector<std::vector<int>> &readArray)
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

    void printNodes()
    {
        for(int i =1;i<nodes.size();i++){
            std::cout<<"Node: "<<i<<std::endl;
            std::cout<<"x: "<<nodes[i][0]<<" y: "<<nodes[i][1]<<" z: "<<nodes[i][2]<<std::endl;
        }
    }

    void printElements()
    {
        for(int i=0; i<elements.size();i++){
            std::cout<<"Element:"<<i<<" ";
            for(int j = 0;j<elements[i].size();j++){
                std::cout<<elements[i][j]<< " ";
            }
            std::cout<<std::endl;
        }
    }

    void printShellList()
    {
        for (int i = 0; i < boundary_types_number; i++)
        {
            std::cout << "Boundary: " << i << std::endl;
            for (int j = 0; j < shell_list[i].size(); j++)
            {
                std::cout<< shell_list[i][j] << " ";
            }
            std::cout<<""<<std::endl;
        }
        
    }
};

class linearSystem{
    private:

    std::vector<std::vector<double>> A;
    std::vector<double> b;

    //Return length of vector, useful for calculating residuals
    double vectorNorm(std::vector<double> &x)
    {
        double sum = 0;
        for (int i = 0; i < x.size(); i++)
        {
            sum += std::pow(x[i],2);
        }
        return std::pow(sum,0.5);
    }
    
    //Return absolute value of a double, useful for calculating residuals
    double abs(double x)
    {
        if (x < 0)
        {
            return -1*x;
        }
        else
        {
            return x;
        }
    }

    public:

    //Class constructor takes the square system matrix and constraint vector Ax = b --> Matrix*x = Vec
    linearSystem(std::vector<std::vector<double>> Matrix, std::vector<double> Vec)
    {
        A = Matrix;
        b = Vec;
    }
    //Solve system of equations, solution outputted to x that is passed in by reference
    //Uses Gauss-siedel method
    void solve(std::vector<double> &x)
    {
        x.resize(b.size(),0.1);

        double residual = 1.0;
        int count = 0;

        while (residual > 0.000001 && count < 10000)
        {
            double normStore = vectorNorm(x);
            for (int i = 0; i < b.size(); i++)
            {
                x[i] = b[i];
                for (int j  = 0; j < b.size(); j++)
                {
                    if (i!=j)
                    {
                        x[i] -= x[j]*A[i][j];
                    }                
                }
                x[i] *= 1/A[i][i];            
            }
            residual = abs(normStore - vectorNorm(x));
            count += 1;
        }
    }
    //Print solution vector (x) to console
    void printSol(std::vector<double> &x)
    {
        for (int i = 0; i < x.size(); i++)
        {
            std::cout<< x[i] <<std::endl;
        }
    }
};

int main()
{
    mesh MESH("mesh.txt");

    //MESH.printNodes();
    //MESH.printElements();
    //MESH.printShellList();
    std::vector<std::vector<double>> Mat = {{1,1,2},{0,1,3},{0,0,1}};
    std::vector<double> Con = {1,2,1};
    std::vector<double> Sol;

    linearSystem testSystem(Mat,Con);
    testSystem.solve(Sol);
    testSystem.printSol(Sol);
    return 0;
}