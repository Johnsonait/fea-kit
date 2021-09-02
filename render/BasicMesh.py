import vtk

points = vtk.vtkPoints()
cells = vtk.vtkCellArray()

mesh_1 = vtk.vtkPolyData()
mesh_2 = vtk.vtkPolyData()

tri = vtk.vtkTriangle()

points.InsertNextPoint((0.0,0.0,0.0))
points.InsertNextPoint((0.0,1.0,0.0))
points.InsertNextPoint((1.0,0.0,0.0))
points.InsertNextPoint((1.0,1.0,1.0))

for i in range(0,3):
    tri.GetPointIds().SetId(i,i) #SetId(i,j) i is local node id, j is global

cells.InsertNextCell(tri)

tri.GetPointIds().SetId(0,3)

cells.InsertNextCell(tri)

mesh_1.SetPoints(points)
mesh_1.SetPolys(cells)

mesh_2.SetPoints(points)
mesh_2.SetPolys(cells)

#Set up cell data (1 cell)
cell_data = vtk.vtkDoubleArray()
cell_data.SetNumberOfComponents(1) #1 component is scalar, 3 is vector, 9 is tensor etc
cell_data.InsertNextTuple([0.5]) #Provide the actual data associated with cell(s)
cell_data.InsertNextTuple([0.25])
mesh_1.GetCellData().SetScalars(cell_data) #Combine the cell data with the mesh itself


#Set up point data (3 points)
point_data = vtk.vtkDoubleArray()
point_data.SetNumberOfComponents(1)
point_data.InsertNextTuple([0.0])
point_data.InsertNextTuple([1.0])
point_data.InsertNextTuple([0.5])
point_data.InsertNextTuple([0.25])
mesh_2.GetPointData().SetScalars(point_data)


mapper = vtk.vtkPolyDataMapper()
mapper.SetInputData(mesh_2)

actor = vtk.vtkActor()
actor.SetMapper(mapper)

window = vtk.vtkRenderWindow()
window.SetSize(500,500)

interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(window)

renderer = vtk.vtkRenderer()
window.AddRenderer(renderer)

renderer.AddActor(actor)
renderer.SetBackground(0.1,0.1,0.4)

window.Render()
interactor.Start()
