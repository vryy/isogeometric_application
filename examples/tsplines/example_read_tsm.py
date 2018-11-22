##################################################################
# test the domain manager in 2D
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

mpatch_import = TSplinePatchTSMImporter()

def CreatePatch(fn):
    patch_ptr = mpatch_import.ImportSingle(fn)

    patch = patch_ptr.GetReference()
    patch.FESpace().Enumerate()
    #print(patch)

    return patch_ptr

def main():
    if len(sys.argv) > 1:
        fn = sys.argv[1]
    else:
        fn = "/home/hbui/workspace2/T-SPLINE/rhino/simple.tsm"

    patch_ptr = CreatePatch(fn)
    patch = patch_ptr.GetReference()
    echo_level = 1
    patch.FESpace().Enumerate()
    patch.FESpace().UpdateCells()
    print(patch)
    print("Import T-Splines from " + fn + " successfully")

if __name__ == "__main__":
    main()

