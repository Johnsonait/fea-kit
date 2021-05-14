#include <iostream>
#include <vector>
#include <fstream>
#include <string>

//using namespace std;

class mesh{
    private:
    std::vector<std::vector<double>> nodes;
    std::vector<std::vector<double>> elements;

    public:
    mesh(std::string fileName){
        std::vector<std::vector<double>> nodal;
        std::vector<std::vector<double>> elemental;
        
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
        while(inFile >> x){
            int nodal_track = 0;
            int elemental_track = 0;

            if(x.front() =='*'){
                instruction = x;
            }
            if(instruction == "*NODE"){
                std::string n,N_x,N_y,N_z;
                inFile >> n;
                inFile >> N_x;
                inFile >> N_y;
                inFile >> N_z;            
                nodal.push_back({std::stod(N_x),std::stod(N_y),std::stod(N_z)});

                nodal_track = nodal_track + 1;
            }
            if(instruction == "*ELEMENT_SHELL"){
                elemental.push_back({});                
            }
        }
        inFile.close();

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
            std::cout<<"Element:"<<i<<std::endl;
            for(int j = 0;j<elements[i].size();j++){
                std::cout<<elements[i][j]<<std::endl;
            }
        }
    }
};


int main(){
    mesh MESH("mesh.txt");

    MESH.printNodes();

    return 0;
}