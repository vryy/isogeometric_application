##################################################################
# test T-splines mesh creation. The example is taken from Fig. 23,
# Isogeometric Analysis using T-splines
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

## create the BSplines multipatch
sys.path.append("../nurbs_multipatch")
import nurbs_ring_patch_from_two_lines
mpatch = nurbs_ring_patch_from_two_lines.CreateMultiPatch()
multipatch_refine_util = MultiPatchRefinementUtility()
multipatch_refine_util.DegreeElevate(mpatch[1], [1, 2]) # elevate degree to cubic

## convert to T-Mesh
Tmesh = TsMesh2D()
util = TSplineUtils()
util.CreateFromBSplines(Tmesh, mpatch[1].GetReference().FESpace())

Tmesh.BuildExtendedTmesh()
print "Tmesh is analysis suitable: ", Tmesh.IsAnalysisSuitable()

util.ExportMatlab(Tmesh, "tmesh_ring_topo.m", "topology")
util.ExportMatlab(Tmesh, "tmesh_ring_knots.m", "knots")

