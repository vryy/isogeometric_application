##################################################################
##################################################################
import sys
import os
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

import geometry_factory
mpatch_util = MultiPatchUtility()
mpatch_export = MultiNURBSPatchMatlabExporter()

arc_ptr = geometry_factory.CreateHalfCircle([0.0, 0.0, 0.0], 'z', 2.1, start_angle=0.0)
arc = arc_ptr.GetReference()
print(arc)
mpatch_export.Export(arc, "half_cicle.m")

## use this command to check in matlab
# for x: arc;nrbplot(nurbs,10);xlabel('x');ylabel('y');zlabel('z');view([1 0 0]);
# for y: arc;nrbplot(nurbs,10);xlabel('x');ylabel('y');zlabel('z');view([0 1 0]);
# for z: arc;nrbplot(nurbs,10);xlabel('x');ylabel('y');zlabel('z');view([0 0 1]);
