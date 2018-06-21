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

mpatch_export = MultiNURBSPatchMatlabExporter()
mpatch_import = MultiNURBSPatchGeoImporter2D()

mp = mpatch_import.Import("../../tests/geo_Lshaped_mp.txt")
print(mp)
#patch = mp[0].GetReference()
#print(patch)

mpatch_export.Export(mp, "geo_Lshaped_mp.m")

