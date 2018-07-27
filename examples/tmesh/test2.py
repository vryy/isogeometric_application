##################################################################
# test T-splines mesh creation. The example is taken from Fig. 23,
# Isogeometric Analysis using T-splines
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

## test the tsplines bezier extraction functionality
# utils = BezierUtils()
# utils.test_tsplines_1() # test the compute_extended_knot_vector


## test the functionality of the Tmesh
Tmesh = TsMesh2D()
util = TSplineUtils()

util.ReadFromFile(Tmesh, 'test2.tmesh')
# print(Tmesh)

Tmesh.BuildExtendedTmesh()
print "Tmesh is analysis suitable: ", Tmesh.IsAnalysisSuitable()

util.ExportMatlab(Tmesh, "tmesh2_topo.m", "topology")
util.ExportMatlab(Tmesh, "tmesh2_knots.m", "knots")

Tmesh.BuildAnchors('test2.coordinates')

Tmesh.BuildCells()

util.ExportMDPA(Tmesh, 'tmesh2.mdpa', 1, 1)

