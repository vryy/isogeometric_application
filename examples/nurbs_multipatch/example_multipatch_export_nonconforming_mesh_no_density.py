##################################################################
##################################################################
import sys
import os
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.StructuralApplication import *
kernel = Kernel()   #defining kernel

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()

import example_multipatch_insert_knots

######################################################################

mpatch = example_multipatch_insert_knots.CreateMultiPatch()
print(mpatch)

iga_mesh = NonConformingMultipatchLagrangeMesh2D(mpatch)
iga_mesh.SetBaseElementName("KinematicLinear")
iga_mesh.SetLastNodeId(1)
iga_mesh.SetLastElemId(1)
iga_mesh.SetLastPropId(1)
iga_mesh.SetUniformDivision(10)

model_part = ModelPart("iga mesh")
iga_mesh.WriteModelPart(model_part)

print(model_part)

#######WRITE TO GID
write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
write_elements = WriteConditionsFlag.WriteConditions
#write_elements = WriteConditionsFlag.WriteElementsOnly
post_mode = GiDPostMode.GiD_PostAscii
multi_file_flag = MultiFileFlag.MultipleFiles
gid_io = StructuralGidIO( "iga_mesh", post_mode, multi_file_flag, write_deformed_flag, write_elements )
gid_io.InitializeMesh( 0.0 )
mesh = model_part.GetMesh()
#self.gid_io.WriteNodeMesh( mesh )
gid_io.WriteMesh( mesh )
print("mesh written...")
gid_io.FinalizeMesh()

