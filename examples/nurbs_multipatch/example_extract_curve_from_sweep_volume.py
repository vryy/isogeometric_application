##################################################################
##################################################################
import sys
import os
write_gid = False
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.IsogeometricBRepApplication import *
if write_gid:
    from KratosMultiphysics.LayerApplication import *
    from KratosMultiphysics.ExternalSolversApplication import * # for model_iga_include
    from KratosMultiphysics.StructuralApplication import * # for model_iga_include
    from KratosMultiphysics.IsogeometricStructuralApplication import * # for model_iga_include
kernel = Kernel()   #defining kernel

import geometry_factory

if write_gid:
    import model_iga_include

mpatch_export3 = MultiNURBSPatchMatlabExporter()
bsplines_patch_util = BSplinesPatchUtility()

def CreatePatch():
    ## generate trajectory curve
    cpoints = []
    cpoints.append([0.0, 0.0, 0.0])
    cpoints.append([0.0, 0.0, 1.0])
    cpoints.append([1.0, 1.0, 2.0])

    order = 2
    curve_ptr = geometry_factory.CreateCurve(cpoints, order)
    curve = curve_ptr.GetReference()
    curve.Prefix = "curve"
    print(curve)
    # mpatch_export3.Export(curve, "curve.m")

    ## generate the local Frenet frame along the trajectory
    npoints = 10
    trans_list = geometry_factory.GenerateLocalFrenetFrame(curve, npoints)
    # for trans in trans_list:
    #     print(trans)
    # geometry_factory.ExportLocalFrenetFrameToMatlab(trans_list, "frame.m", 2e-1)

    ## generate a list of cut section
    rin = 0.5
    rout = 1.0
    start_angle = 45
    end_angle = 90
    axis = "z"
    ring_patches = []
    for i in range(0, npoints):
        ring_ptr = geometry_factory.CreateSmallRing([0.0, 0.0, 0.0], axis, rin, rout, start_angle, end_angle)
        ring = ring_ptr.GetReference()
        ring.Id = i+1
        ring.Prefix = "ring"
        ring.ApplyTransformation(trans_list[i])
        # mpatch_export3.Export(ring, ring.Name() + "_def.m")
        ring_patches.append(ring_ptr)

    ## create the sweep volume
    order_w = 2
    vol_ptr = bsplines_patch_util.CreateLoftPatchFromList3D(ring_patches, order_w)
    return vol_ptr

def main():
    vol_ptr = CreatePatch()
    vol = vol_ptr.GetReference()
    vol.Id = 1
    vol.LayerIndex = 1
    vol.Prefix = "volume"
    # mpatch_export3.Export(vol, "sweep_volume.m")

    if write_gid:
        params_post = {}
        params_post['name'] = "sweep_volume"
        mpatch = MultiPatch3D()
        mpatch.AddPatch(vol)
        model_iga_include.PostMultiPatch(mpatch, 3, 0.0, params_post)

    spatch_ptr = vol.ConstructSlicedPatch(0, 0.5)
    spatch = spatch_ptr.GetReference()

    lpatch_ptr = spatch.ConstructSlicedPatch(0, 0.5)
    lpatch = lpatch_ptr.GetReference()
    print(lpatch)

    pcurve = PatchCurve(lpatch_ptr)
    print(pcurve.GetValue(0.5))
    print(pcurve.GetDerivative(0, 0.5))
    print(pcurve.GetSecondDerivative(0, 0, 0.5))

    pcurve[CURVE_SEARCH_TOLERANCE] = 1e-8
    pcurve[CURVE_NUMBER_OF_SAMPLING] = 20
    print(pcurve.ComputeLength(0.0, 1.0))

    # create a sample uniform distribution
    tvec = pcurve.ComputeUniformDivision(0.0, 1.0, 10)
    print(tvec)
    for t in tvec:
        P = pcurve.GetValue(t)
        print("%f,%f,%f" % (P[0], P[1], P[2]))

if __name__ == "__main__":
    main()
