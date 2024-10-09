##################################################################
##################################################################
import sys
import os
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

mpatch_export = MultiNURBSPatchMatlabExporter()
mpatch_import = MultiNURBSPatchGeoImporter2D()

mp = mpatch_import.Import("../../tests/geo_Lshaped_mp.txt")
print(mp)
#patch = mp[0].GetReference()
#print(patch)

mpatch_export.Export(mp, "geo_Lshaped_mp.m")

