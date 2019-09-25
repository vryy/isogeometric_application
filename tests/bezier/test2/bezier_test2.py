import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##setting up paths
kratos_libs_path = kratos_root_path+'libs' ##kratos_root/libs
kratos_applications_path = kratos_root_path+'applications' ##kratos_root/applications
##################################################################
##################################################################
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.DiscontinuitiesApplication import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

model_part = ModelPart("isogeometric_model")
model_part.AddNodalSolutionStepVariable(DISPLACEMENT)

model_part_io = IsogeometricModelPartIO("bezier_test2")
model_part_io.ReadModelPart(model_part)

model_part.Properties[1].SetValue(DENSITY, 0 )
model_part.Properties[1].SetValue(YOUNG_MODULUS, 1e+07 )
model_part.Properties[1].SetValue(POISSON_RATIO, 0.48 )
model_part.Properties[1].SetValue(THICKNESS, 1 )
model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, PlaneStressSD() )

print model_part

scheme = ResidualBasedIncrementalUpdateStaticScheme()
scheme.InitializeElements(model_part)

print "Test1: probing shape function values"
test_utils = NURBSTestUtils()
test_utils.Test2(model_part)



