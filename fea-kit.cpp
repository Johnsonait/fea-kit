#include <iostream>
#include <vector>

class mesh{
    private:
    std::vector<std::vector<double>> nodes;
    std::vector<std::vector<double>> elements;

    public:
    mesh(std::vector<std::vector<double>> nod, std::vector<std::vector<double>> ele){
        nodes = nod;
        elements = ele;
    }
    void printNodes(){
        for(int i =1;i<nodes.size();i++){
            std::cout<<"Node: "<<i<<std::endl;
            for(int j = 0;j<nodes[i].size();j++){
                std::cout<<nodes[i][j]<<std::endl;
            }
        }
    }
};

int main(){
    std::vector<std::vector<double>> nodal = {{1,2,3},{0,0,0},{1,1,1},{2,2,2}};
    std::vector<std::vector<double>> elemental = {{1,2},{2,3}};
    
    mesh MESH(nodal,elemental);

    MESH.printNodes();

    return 0;
}