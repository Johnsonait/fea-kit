import vtk
import csv
from numpy import loadtxt

global_nodes = []
global_elements = []
global_u = []
global_strain = []
equiv_strain = []
global_stress = []
equiv_stress = []
global_temp = []

#Read in data from file
def Read(results,model):
    #Open file and read each line in to lines
    #Note that the first line is fluff, the rest are dedicated to data for each node
    with open(results) as results_file:
        lines = loadtxt(results, comments="#", delimiter=",", unpack=False)
        for line in lines:
            #Read line data as expected
            x = [line[0], line[1], line[2]] 
            u = [line[3], line[4], line[5]] 
            e = [line[6], line[7], line[8],line[9], line[10], line[11]]
            e_eq = line[12]
            s = [line[13], line[14], line[15],line[16], line[17], line[18]]
            s_eq = line[19]
            T = line[20]

            #Append data to global vars
            global_nodes.append(x)
            global_u.append(u)
            global_strain.append(e)
            equiv_strain.append(e_eq)
            global_stress.append(s)
            equiv_stress.append(s_eq)
            global_temp.append(T)
    #Most global data has been read, only remaining data is the element data to store global ids for each element
    #This is stored as part of the model (not result) data and so must be read from the model file
    with open(model) as model_file:
        line = ["0"]
        while len(line) > 0:
            line = model_file.readline().split()
            if len(line) > 0 and line[0] == "*ELEMENT_SOLID":
                line = model_file.readline().split()
                while line != ['*'] and len(line)>0:
                    temp = []
                    for i in range(2,10):
                        temp.append(int(line[i].replace(',','')))
                    line = model_file.readline().split()
                    
                    global_elements.append(temp)
    
def PrintData():
    print("Nodes: ")
    for node in global_nodes:
        print(node)
    print("Elements: ")
    for element in global_elements:
        print(element)
    print("Displacements: ")
    for u in global_u:
        print(u)
    print("Strains: ")
    for e in global_strain:
        print(e)
    print("Equivalent Strains: ")
    for e_eq in equiv_strain:
        print(e_eq)
    print("Stresses: ")
    for s in global_stress:
        print(s)
    print("Equivalent Stresses: ")
    for s_eq in equiv_stress:
        print(s_eq)
    print("Temperatures: ")
    for t in global_temp:
        print(t)
    

def main():
    #First read data from "results.txt" into global variables
    Read("results.txt","model.txt")
    #PrintData()

    field = []
    for row in range(0,len(global_u)):
        field.append(global_u[row][2])
    print(field)
    #This is where VTK starts to work its magic
    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()

    mesh = vtk.vtkPolyData()


    for i in range(0,len(global_nodes)):
        points.InsertNextPoint(global_nodes[i])
    
    hexahedron = vtk.vtkHexahedron()
    
    uGrid = vtk.vtkUnstructuredGrid()
    uGrid.SetPoints(points)

    #Set up each element from global_elements and add it to cells and grid
    for el in range(0,len(global_elements)):
        for i in range(0,8):
            hexahedron.GetPointIds().SetId(i,global_elements[el][i]-1) #SetId(i,j) i is local node id, j is global 
        cells.InsertNextCell(hexahedron)
        uGrid.InsertNextCell(hexahedron.GetCellType(), hexahedron.GetPointIds())


    mesh.SetPoints(points)
    mesh.SetPolys(cells)

    point_data = vtk.vtkDoubleArray()
    point_data.SetNumberOfComponents(1)
    
    for i in range(0,len(field)):
        point_data.InsertNextTuple([field[i]])
    
    #mesh.GetPointData().SetScalars(point_data)
    uGrid.GetPointData().SetScalars(point_data)

    #Create lookup table for color mapping scalar data
    lu_table = vtk.vtkLookupTable()
    lu_table.SetNumberOfColors(256)
    lu_table.SetHueRange(0.0,0.667)
    lu_table.Build()
    
    #mapper = vtk.vtkPolyDataMapper()
    #mapper.SetInputData(mesh)
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(uGrid)
    mapper.SetLookupTable(lu_table)
    mapper.SetScalarRange(min(field),max(field))
    
    # a colorbar to display the colormap
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable( mapper.GetLookupTable() )
    scalar_bar.SetTitle("Point scalar value")
    scalar_bar.SetOrientationToHorizontal()
    scalar_bar.GetLabelTextProperty().SetColor(1,1,1) 
    scalar_bar.GetTitleTextProperty().SetColor(1,1,1)

    # position it in window
    coord = scalar_bar.GetPositionCoordinate()
    coord.SetCoordinateSystemToNormalizedViewport()
    coord.SetValue(0.1,0.05)
    scalar_bar.SetWidth(.8)
    scalar_bar.SetHeight(.15)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetEdgeVisibility(True) 

    window = vtk.vtkRenderWindow()
    window.SetSize(500,500)

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)

    renderer = vtk.vtkRenderer()
    window.AddRenderer(renderer)

    renderer.AddActor(actor)
    renderer.AddActor(scalar_bar)
    renderer.SetBackground(0,0,0)

    window.Render()
    interactor.Start()

main()
