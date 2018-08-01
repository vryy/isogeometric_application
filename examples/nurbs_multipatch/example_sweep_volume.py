##################################################################
######################## include.py   ############################
##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Institute for Structural Mechanics, RUB   #####
##### all rights reserved                                    #####
##################################################################
##################################################################
##################################################################
##################################################################
import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

import geometry_factory

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
    mpatch_export3.Export(curve, "curve.m")

    ## generate the local Frenet frame along the trajectory
    npoints = 10
    trans_list = geometry_factory.GenerateLocalFrenetFrame(curve, npoints)
    for trans in trans_list:
        print(trans)
    geometry_factory.ExportLocalFrenetFrameToMatlab(trans_list, "frame.m", 2e-1)

    ## generate a list of cut section
    rin = 0.5
    rout = 1.0
    start_angle = 45
    end_angle = 90
    axis = "z"
    cnt = 1
    ring_patches = []
    for i in range(0, npoints):
        ring_ptr = geometry_factory.CreateSmallRing([0.0, 0.0, 0.0], axis, rin, rout, start_angle, end_angle)
        ring = ring_ptr.GetReference()
        ring.Id = cnt
        cnt = cnt+1
        ring.Prefix = "ring"
        ring.ApplyTransformation(trans_list[i])
        # mpatch_export3.Export(ring, ring.Name() + "_def.m")
        ring_patches.append(ring_ptr)

    ## create the sweep volume
    order_w = 2
    vol_ptr = bsplines_patch_util.CreateConnectedPatch(ring_patches, order_w)
    return vol_ptr

def main():
    vol_ptr = CreatePatch()
    vol = vol_ptr.GetReference()
    vol.Id = 1
    vol.Prefix = "volume"
    mpatch_export3.Export(vol, "sweep_volume.m")

if __name__ == "__main__":
    main()

