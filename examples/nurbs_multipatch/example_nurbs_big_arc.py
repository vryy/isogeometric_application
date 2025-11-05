##################################################################
import sys
import os
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel
##################################################################

import geometry_factory

mpatch_export3 = MultiNURBSPatchMatlabExporter()

center = [0.0, 0.0, 0.0]
axis = 'x'
r = 10.0
start_angle = -45.0 # 15.0
end_angle = 95.0 # 45.0
arc_ptr = geometry_factory.CreateSmallArc(center, axis, r, start_angle, end_angle)

mpatch = MultiPatch1D()
mpatch.AddPatch(arc_ptr)
mpatch.Enumerate()

mpatch_export3.Export(mpatch, "arc.m")
