import vtk
import os
import sys 

pdbID, vtkform = sys.argv[1].split('.')
inputfile = pdbID + '.vtk'
output_file = pdbID + "_Surface_Potential.vtk"


iso_scalar_name = "ls"


reader = vtk.vtkDataSetReader()
reader.SetFileName(inputfile)
reader.ReadAllScalarsOn()
reader.Update()


iso_surface = vtk.vtkContourFilter()
iso_surface.SetInputConnection(reader.GetOutputPort())
iso_surface.SetValue(0, 0.0)

iso_surface.SetInputArrayToProcess(
    0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, iso_scalar_name
)
iso_surface.Update()



writer = vtk.vtkDataSetWriter()
writer.SetFileName(output_file)
writer.SetFileTypeToBinary()
writer.SetInputConnection(iso_surface.GetOutputPort())
writer.Write()


