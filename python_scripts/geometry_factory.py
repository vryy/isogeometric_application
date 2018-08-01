import math
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *

###
### This module is a factory to generate typical geometries for isogeometric analysis, e.g. circle, l-shape, ...
###

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()
bsplines_patch_util = BSplinesPatchUtility()

### Compute cross product
def cross(c, a, b):
    c[0] = a[1]*b[2] - a[2]*b[1]
    c[1] = a[2]*b[0] - a[0]*b[2]
    c[2] = a[0]*b[1] - a[1]*b[0]
    return c

### Compute dot product
def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

### Normalize a vector
def normalize(a):
    norma = math.sqrt(a[0]**2 + a[1]**2 + a[2]**2)
    a[0] = a[0] / norma
    a[1] = a[1] / norma
    a[2] = a[2] / norma
    return a

### Create a line from start_point to end_point with knot vector [0 0 0 ... 1 1 1]
### On output the pointer to the patch will be returned
def CreateLine(start_point, end_point, order = 1):
    Id = 0
    fes = nurbs_fespace_library.CreatePrimitiveFESpace(order)
    ctrl_grid = grid_lib.CreateLinearControlPointGrid(start_point[0], start_point[1], start_point[2], fes.Number(0), end_point[0], end_point[1], end_point[2])
    patch_ptr = multipatch_util.CreatePatchPointer(Id, fes)
    patch = patch_ptr.GetReference()
    patch.CreateControlPointGridFunction(ctrl_grid)
    return patch_ptr

### Create a curve from the point list, given as [ [x0, y0, z0], ... ]
### All the weight is assumed 1
def CreateCurve(points, order):
    Id = 0
    number = len(points)
    fes = nurbs_fespace_library.CreateUniformFESpace(number, order)
    ctrl_grid = StructuredControlPointGrid1D(number)
    for i in range(0, number):
        ctrl_grid.SetValue(i, ControlPoint(points[i][0], points[i][1], points[i][2], 1.0))
    curve_ptr = multipatch_util.CreatePatchPointer(Id, fes)
    curve = curve_ptr.GetReference()
    curve.CreateControlPointGridFunction(ctrl_grid)
    return curve_ptr

### Create an arc at center on the surface perpendicular to the given axis. By default, the quadratic arc is generated. The knot vector will be [0 0 0 1 1 1]
### On output the pointer to the patch will be returned. Small arc means that the open angle is less than 90 degrees.
def CreateSmallArc(center, axis, radius, start_angle, end_angle):
    ## firstly create an arc in xy plane at (0, 0)
    Id = 0
    fes = nurbs_fespace_library.CreatePrimitiveFESpace(2)
    ctrl_grid = grid_lib.CreateLinearControlPointGrid(0.0, 0.0, 0.0, fes.Number(0), radius, 0.0, 0.0)

    sweep = end_angle - start_angle
#    if abs(sweep > 90):
#        raise ValueError('the open angle must be in [-90, 90] degrees, sweep =', sweep)

    dsweep = 0.5*sweep/180.0*math.pi
    wm = math.cos(dsweep)
    x = radius*wm
    y = radius*math.sin(dsweep)
    xm = x + y*math.tan(dsweep)

    if axis == 'z':
        trans = RotationZ(start_angle + 0.5*sweep)
    elif axis == 'y':
        trans = RotationZ(start_angle + 0.5*sweep)
        trans.AppendTransformation(RotationX(90.0))
    elif axis == 'x':
        trans = RotationZ(start_angle + 0.5*sweep + 90.0)
        trans.AppendTransformation(RotationY(90.0))
    trans.AppendTransformation(Translation(center[0], center[1], center[2]))

    pt1 = ctrl_grid[0]
    pt1.WX = x
    pt1.WY = -y
    pt1.WZ = 0.0
    pt1.W = 1.0
    pt1.ApplyTransformation(trans)
    ctrl_grid[0] = pt1

    pt2 = ctrl_grid[1]
    pt2.WX = wm*xm
    pt2.WY = 0.0
    pt2.WZ = 0.0
    pt2.W = wm
    pt2.ApplyTransformation(trans)
    ctrl_grid[1] = pt2

    pt3 = ctrl_grid[2]
    pt3.WX = x
    pt3.WY = y
    pt3.WZ = 0.0
    pt3.W = 1.0
    pt3.ApplyTransformation(trans)
    ctrl_grid[2] = pt3

    patch_ptr = multipatch_util.CreatePatchPointer(Id, fes)
    patch = patch_ptr.GetReference()
    patch.CreateControlPointGridFunction(ctrl_grid)
    return patch_ptr

### Create a 2D ring at center on the surface perpendicular to the axis. By default, the quadratic arc is generated. The knot vector will be [0 0 0 1 1 1]
### On output the pointer to the patch will be returned. Small ring means that the open angle is less than 90 degrees.
def CreateSmallRing(center, axis, rin, rout, start_angle, end_angle):
    ## create inner arc
    iarc_ptr = CreateSmallArc(center, axis, rin, start_angle, end_angle)
    iarc = iarc_ptr.GetReference()

    ## create outer arc
    oarc_ptr = CreateSmallArc(center, axis, rout, start_angle, end_angle)
    oarc = oarc_ptr.GetReference()

    ## create ring
    ring_patch_ptr = bsplines_patch_util.CreateConnectedPatch(iarc, oarc)
    return ring_patch_ptr

### Create a list of Frenet frame along a curve. The Frenet frame is stored as a transformation matrix.
### zvec is a reference vector to compute B at the first sampling point. It shall not be parallel with the tangent vector of the first sampling point.
def GenerateLocalFrenetFrame(curve, num_sampling_points, zvec = [1.0, 0.0, 0.0]):
    trans_list = []
    B = Array3()
    ctrl_pnt_grid_func = curve.GridFunction(CONTROL_POINT_COORDINATES)
    for i in range(0, num_sampling_points):
        xi = float(i) / (num_sampling_points-1)
        pnt = [xi, 0.0, 0.0]
        P = ctrl_pnt_grid_func.GetValue(pnt)
        T = ctrl_pnt_grid_func.GetDerivative(pnt)
        T = normalize(T[0])

        if i == 0:
            cross(B, zvec, T)
            B = normalize(B)
        else:
            B = B - dot(B, T)*T
            B = normalize(B)

        trans = Transformation(B, T, P)
        trans_list.append(trans)

    return trans_list

def ExportLocalFrenetFrameToMatlab(trans_list, fn, s = 1.0):
    ifile = open(fn, "w")
    cnt = 1
    ifile.write("s = " + str(s) + ";\n")
    for trans in trans_list:
        P = trans.P()
        B = trans.V1()
        N = trans.V2()
        T = trans.V3()
        ifile.write("C{" + str(cnt) + "} = [" + str(P[0]) + " " + str(P[1]) + " " + str(P[2]) + "];\n")
        ifile.write("T{" + str(cnt) + "} = [" + str(T[0]) + " " + str(T[1]) + " " + str(T[2]) + "];\n")
        ifile.write("B{" + str(cnt) + "} = [" + str(B[0]) + " " + str(B[1]) + " " + str(B[2]) + "];\n")
        ifile.write("N{" + str(cnt) + "} = [" + str(N[0]) + " " + str(N[1]) + " " + str(N[2]) + "];\n")
        ifile.write("hold on; plot_frame(C{" + str(cnt) + "}, B{" + str(cnt) + "}, N{" + str(cnt) + "}, T{" +str(cnt) + "}, s);\n")
        cnt = cnt + 1
    ifile.close()

