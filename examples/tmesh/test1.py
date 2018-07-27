##################################################################
# test T-splines mesh creation. The example is taken from Fig. 23,
# Isogeometric Analysis using T-splines
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

Tmesh = TsMesh2D()
util = TSplineUtils()

util.ReadFromFile(Tmesh, 'test1.tmesh')
# print(Tmesh)

# 
# u_knots = Vector()
# v_knots = Vector()
# 
# u = 8
# v = 6
# print("knots at " + str(u) + "," + str(v))
# Tmesh.FindKnots2(u, v, u_knots, v_knots)
# print(u_knots)
# print(v_knots)

Tmesh.BuildExtendedTmesh()
# Tmesh.BuildExtendedTmesh()
print "Tmesh is analysis suitable: ", Tmesh.IsAnalysisSuitable()

util.ExportMatlab(Tmesh, "tmesh1_topo.m", "topology")
util.ExportMatlab(Tmesh, "tmesh1_knots.m", "knots")

