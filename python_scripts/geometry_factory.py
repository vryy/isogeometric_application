import math
from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.IsogeometricApplication import *

###
### This module is a factory to generate typical geometries for isogeometric analysis, e.g. circle, l-shape, ...
###

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()
multipatch_refine_util = MultiPatchRefinementUtility()
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

### Compute Gaussian function
def gaussian(mu, sigma, x):
    return math.exp(-0.5*((x-mu)/sigma)**2)/sigma/math.sqrt(2.0*math.pi)

### Compute inverse Gaussian function
def inv_gaussian1(mu, sigma, g):
    return -sigma*math.sqrt(-2.0*math.log(sigma*math.sqrt(2*math.pi)*g)) + mu

### Compute inverse Gaussian function
def inv_gaussian2(mu, sigma, g):
    return sigma*math.sqrt(-2.0*math.log(sigma*math.sqrt(2*math.pi)*g)) + mu

# ### Generate distributed Gaussian array in span (min, max). It is useful to generate a knot vector with Gaussian distribution for testing
# def GenerateGaussianArray(half_n, min_k, max_k, sigma):
#     mu = 0.5*(min_k + max_k)
#     max_g = gaussian(mu, sigma, mu)
#     min_g = gaussian(mu, sigma, 0.0)
#     print("min_g:", min_g)
#     print("max_g:", max_g)
#     print("mu:", inv_gaussian1(mu, sigma, max_g))
#     print("mu:", inv_gaussian2(mu, sigma, max_g))

#     k_list = []

#     for i in range(0, half_n+1):
#         t = float(i+1)/(half_n+1)
#         g = t*(max_g-min_g) + min_g
#         # print("g:", g)
#         k = inv_gaussian1(mu, sigma, g)
#         # print("k:", k)
#         k_list.append(k)

#     for i in range(0, half_n):
#         t = float(half_n-i)/(half_n+1)
#         g = t*(max_g-min_g) + min_g
#         k = inv_gaussian2(mu, sigma, g)
#         k_list.append(k)

#     return k_list

# ### Generate distributed Gaussian array in span (min, max). It is useful to generate a knot vector with Gaussian distribution for testing
# def GenerateGaussianArray(n, min_k, max_k):
#     mu = 0.0
#     sigma = 1.0
#     k_list = []
#     g_list = []
#     min_g = 0.0
#     max_g = gaussian(mu, sigma, mu)
#     sum_g = 0.0
#     for i in range(0, n):
#         t = 6.0*float(i)/(n-1) - 3.0;
#         g = gaussian(mu, sigma, t)
#         g_list.append(g)
#         sum_g = sum_g + g
#     t = 0.0
#     for g in g_list:
#         t = t + g
#         k_list.append(t/sum_g)
#     return k_list

### Generate distributed Gaussian array in span (min, max). It is useful to generate a knot vector with Gaussian distribution for testing
def GenerateGaussianArray(num_span, max_elem_in_span, sigma, min_k, max_k):
    k_list = []
    # make a span in [-3, 3]
    for i in range(0, num_span):
        min_t = -3.0 + float(i)/num_span*6.0
        max_t = -3.0 + float(i+1)/num_span*6.0
        # print("min_t:", min_t)
        # print("max_t:", max_t)
        # get a samping value in [min_t, max_t]
        t = 0.5*(min_t + max_t)
        # print("t:", t)
        g = gaussian(0.0, sigma, t*sigma)
        # print("g:", g)
        n = int(max_elem_in_span*g)
        # print("n:", n)
        # generate n numbers from min_t to max_t
        for j in range(0, n):
            t = min_t + float(j+0.5)/n*(max_t-min_t)
            t_scale = (t+3.0)/6.0;
            k_list.append(t_scale*(max_k-min_k)+min_k)

    return k_list

### Create a line from start_point to end_point with knot vector [0 0 0 ... 1 1 1]
### On output the pointer to the patch will be returned
def CreateLine(start_point, end_point, order = 1, sample_patch = Patch1DSelector.RealPatch):
    Id = 0
    fes = nurbs_fespace_library.CreatePrimitiveFESpace(order)
    ctrl_grid = grid_lib.CreateLinearControlPointGrid(start_point[0], start_point[1], start_point[2], fes.Number(0), end_point[0], end_point[1], end_point[2])
    patch_ptr = sample_patch.Create(Id, fes)
    patch = patch_ptr.GetReference()
    patch.CreateControlPointGridFunction(ctrl_grid)
    return patch_ptr

### Create a curve from the control point list, given as [ [x0, y0, z0], ... ]
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
def CreateSmallArc(center, axis, radius, start_angle, end_angle, sample_patch=Patch1DSelector.RealPatch):
    ## firstly create an arc in xy plane at (0, 0)
    Id = 0
    fes = nurbs_fespace_library.CreatePrimitiveFESpace(2)
    ctrl_grid = grid_lib.CreateLinearControlPointGrid(0.0, 0.0, 0.0, fes.Number(0), radius, 0.0, 0.0)

    sweep = end_angle - start_angle
    if (sweep > 180.0):
        raise Exception("Small arc cannot be larger than 180 degres, the input angle is %f" % (sweep))

    if (sweep == 180.0):
        print("WARNING!!!the control weight of the curve will be zero when angle=180. One must do h-refinement later to fix it.")

    # Note: the sweep angle=180.0 still works, but one must do h-refinement to avoid zero weight
    # Alternatively, use CreateHalfCircle (initial order will be 3)

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

    patch_ptr = sample_patch.Create(Id, fes)
    patch = patch_ptr.GetReference()
    patch.CreateControlPointGridFunction(ctrl_grid)
    return patch_ptr

### Create a single curve for NURBS half circle. Reference: https://www.geometrictools.com/Documentation/NURBSCircleSphere.pdf
### Note that the default degree is 3
def CreateHalfCircle(center, axis, radius, start_angle=0.0):
    Id = 0
    fes = nurbs_fespace_library.CreatePrimitiveFESpace(3)
    ctrl_grid = grid_lib.CreateLinearControlPointGrid(0.0, 0.0, 0.0, fes.Number(0), radius, 0.0, 0.0)

    if axis == 'z':
        trans = RotationZ(start_angle)
    elif axis == 'y':
        trans = RotationZ(start_angle)
        trans.AppendTransformation(RotationX(90.0))
    elif axis == 'x':
        trans = RotationZ(start_angle + 90.0)
        trans.AppendTransformation(RotationY(90.0))
    trans.AppendTransformation(Translation(center[0], center[1], center[2]))

    pt1 = ctrl_grid[0]
    pt1.WX = radius
    pt1.WY = 0.0
    pt1.WZ = 0.0
    pt1.W = 1.0
    pt1.ApplyTransformation(trans)
    ctrl_grid[0] = pt1

    pt2 = ctrl_grid[1]
    pt2.WX = 1.0/3*radius
    pt2.WY = 2.0/3*radius
    pt2.WZ = 0.0
    pt2.W = 1.0/3
    pt2.ApplyTransformation(trans)
    ctrl_grid[1] = pt2

    pt3 = ctrl_grid[2]
    pt3.WX = -1.0/3*radius
    pt3.WY = 2.0/3*radius
    pt3.WZ = 0.0
    pt3.W = 1.0/3
    pt3.ApplyTransformation(trans)
    ctrl_grid[2] = pt3

    pt4 = ctrl_grid[3]
    pt4.WX = -radius
    pt4.WY = 0.0
    pt4.WZ = 0.0
    pt4.W = 1.0
    pt4.ApplyTransformation(trans)
    ctrl_grid[3] = pt4

    patch_ptr = multipatch_util.CreatePatchPointer(Id, fes)
    patch = patch_ptr.GetReference()
    patch.CreateControlPointGridFunction(ctrl_grid)
    return patch_ptr

### Create a 2D ring at center on the surface perpendicular to the axis. By default, the quadratic arc is generated. The knot vector will be [0 0 0 1 1 1]
### On output the pointer to the patch will be returned. Small ring means that the open angle is typically less than 180 degrees.
def CreateSmallRing(center, axis, rin, rout, start_angle, end_angle):
    ## create inner arc
    iarc_ptr = CreateSmallArc(center, axis, rin, start_angle, end_angle)
    iarc = iarc_ptr.GetReference()

    ## create outer arc
    oarc_ptr = CreateSmallArc(center, axis, rout, start_angle, end_angle)
    oarc = oarc_ptr.GetReference()

    ## create ring
    ring_patch_ptr = bsplines_patch_util.CreateLoftPatch(iarc, oarc)
    return ring_patch_ptr

### Create the 2D rectangle aligned with Cartesian axes
def CreateRectangle(start_point, end_point):
    line1 = CreateLine(start_point, [end_point[0], start_point[1], start_point[2]])
    line2 = CreateLine([start_point[0], end_point[1], start_point[2]], [end_point[0], end_point[1], start_point[2]])
    face_ptr = bsplines_patch_util.CreateLoftPatch(line1, line2)
    return face_ptr

### Create the 2D parallelogram
###   P4---P3
###   |    |
###   P1---P2
def CreateParallelogram(P1, P2, P3, P4):
    line1 = CreateLine(P1, P2)
    line2 = CreateLine(P4, P3)
    face_ptr = bsplines_patch_util.CreateLoftPatch(line1, line2)
    return face_ptr

### Create the 3D slab aligned with Cartesian axes
def CreateSlab(start_point, end_point):
    line1 = CreateLine(start_point, [end_point[0], start_point[1], start_point[2]])
    line2 = CreateLine([start_point[0], end_point[1], start_point[2]], [end_point[0], end_point[1], start_point[2]])
    face1_ptr = bsplines_patch_util.CreateLoftPatch(line1, line2)
    face1 = face1_ptr.GetReference()

    line3 = CreateLine([start_point[0], start_point[1], end_point[2]], [end_point[0], start_point[1], end_point[2]])
    line4 = CreateLine([start_point[0], end_point[1], end_point[2]], [end_point[0], end_point[1], end_point[2]])
    face2_ptr = bsplines_patch_util.CreateLoftPatch(line3, line4)
    face2 = face2_ptr.GetReference()

    volume_patch_ptr = bsplines_patch_util.CreateLoftPatch(face1, face2)
    return volume_patch_ptr

### Create the 3D parallelpiped
def CreateParallelepiped(start_point, dir1, dir2, dir3):
    # l1: start -> start + d1
    line1 = CreateLine(start_point, [start_point[0] + dir1[0], start_point[1] + dir1[1], start_point[2] + dir1[2]])
    # l2: start + d2 -> start + d1 + d2
    line2 = CreateLine([start_point[0] + dir2[0], start_point[1] + dir2[1], start_point[2] + dir2[2]], [start_point[0] + dir1[0] + dir2[0], start_point[1] + dir1[1] + dir2[1], start_point[2] + dir1[2] + dir2[2]])
    face1_ptr = bsplines_patch_util.CreateLoftPatch(line1, line2)
    face1 = face1_ptr.GetReference()

    # l3: start + d3 -> start + d1 + d3
    line3 = CreateLine([start_point[0] + dir3[0], start_point[1] + dir3[1], start_point[2] + dir3[2]], [start_point[0] + dir1[0] + dir3[0], start_point[1] + dir1[1] + dir3[1], start_point[2] + dir1[2] + dir3[2]])
    # l4: start + d2 + d3 -> start + d1 + d2 + d3
    line4 = CreateLine([start_point[0] + dir2[0] + dir3[0], start_point[1] + dir2[1] + dir3[1], start_point[2] + dir2[2] + dir3[2]], [start_point[0] + dir1[0] + dir2[0] + dir3[0], start_point[1] + dir1[1] + dir2[1] + dir3[1], start_point[2] + dir1[2] + dir2[2] + dir3[2]])
    face2_ptr = bsplines_patch_util.CreateLoftPatch(line3, line4)
    face2 = face2_ptr.GetReference()

    volume_patch_ptr = bsplines_patch_util.CreateLoftPatch(face1, face2)
    return volume_patch_ptr

### Create a half circle with 4 patches configuration
def CreateHalfCircle4(center, axis, radius, rotation_angle, params={}):

    if 'make_interface' in params:
        make_interface = params['make_interface']
    else:
        make_interface = True

    if 'square_control' in params:
        square_control = params['square_control']
    else:
        square_control = 1.0/3

    ### create arcs
    arc1_ptr = CreateSmallArc(center, axis, radius, 0.0, 45.0)
    arc1 = arc1_ptr.GetReference()

    arc2_ptr = CreateSmallArc(center, axis, radius, 45.0, 135.0)
    arc2 = arc2_ptr.GetReference()

    arc3_ptr = CreateSmallArc(center, axis, radius, 135.0, 180.0)
    arc3 = arc3_ptr.GetReference()

    square_size = square_control*radius

    ### create lines
    if axis == 'x':
        p1 = [center[0], center[1] + square_size, center[2]]
        p2 = [center[0], center[1] + square_size, center[2] + square_size]
        p3 = [center[0], center[1] - square_size, center[2]]
        p4 = [center[0], center[1] - square_size, center[2] + square_size]
    elif axis == 'y':
        p1 = [center[0], center[1], center[2] + square_size]
        p2 = [center[0] + square_size, center[1], center[2] + square_size]
        p3 = [center[0], center[1], center[2] - square_size]
        p4 = [center[0] + square_size, center[1], center[2] - square_size]
    elif axis == 'z':
        p1 = [center[0] + square_size, center[1], center[2]]
        p2 = [center[0] + square_size, center[1] + square_size, center[2]]
        p3 = [center[0] - square_size, center[1], center[2]]
        p4 = [center[0] - square_size, center[1] + square_size, center[2]]

    u_order = arc1.Order(0)

    line1_ptr = CreateLine(p1, p2, u_order)
    line1 = line1_ptr.GetReference()

    line2_ptr = CreateLine(p2, p4, u_order)
    line2 = line2_ptr.GetReference()

    line3_ptr = CreateLine(p4, p3, u_order)
    line3 = line3_ptr.GetReference()

    line4_ptr = CreateLine(p1, p3, u_order)
    line4 = line4_ptr.GetReference()

    patch1_ptr = bsplines_patch_util.CreateLoftPatch(arc1, line1)
    patch2_ptr = bsplines_patch_util.CreateLoftPatch(arc2, line2)
    patch3_ptr = bsplines_patch_util.CreateLoftPatch(arc3, line3)
    patch4_ptr = bsplines_patch_util.CreateLoftPatchFromList2D([line2, line4], 1)
    multipatch_refine_util.DegreeElevate(patch4_ptr, [0, u_order-1])

    patch1 = patch1_ptr.GetReference()
    patch1.Id = 1
    patch2 = patch2_ptr.GetReference()
    patch2.Id = 2
    patch3 = patch3_ptr.GetReference()
    patch3.Id = 3
    patch4 = patch4_ptr.GetReference()
    patch4.Id = 4

    print("patch1:" + str(patch1))
    print("patch4:" + str(patch4))

    if rotation_angle != 0.0:
        trans = Transformation()
        trans.AppendTransformation(Translation(-center[0], -center[1], -center[2]))
        if axis == 'z':
            trans.AppendTransformation(RotationZ(rotation_angle))
        elif axis == 'y':
            trans.AppendTransformation(RotationY(rotation_angle))
        elif axis == 'x':
            trans.AppendTransformation(RotationX(rotation_angle))
        trans.AppendTransformation(Translation(center[0], center[1], center[2]))
        patch1.ApplyTransformation(trans)
        patch2.ApplyTransformation(trans)
        patch3.ApplyTransformation(trans)
        patch4.ApplyTransformation(trans)

    if make_interface:
        bsplines_patch_util.MakeInterface(patch1, BoundarySide2D.U1, patch2, BoundarySide2D.U0, BoundaryDirection.Forward)
        bsplines_patch_util.MakeInterface(patch2, BoundarySide2D.U1, patch3, BoundarySide2D.U0, BoundaryDirection.Forward)
        bsplines_patch_util.MakeInterface(patch1, BoundarySide2D.V1, patch4, BoundarySide2D.U0, BoundaryDirection.Reversed)
        bsplines_patch_util.MakeInterface(patch2, BoundarySide2D.V1, patch4, BoundarySide2D.V0, BoundaryDirection.Forward)
        bsplines_patch_util.MakeInterface(patch3, BoundarySide2D.V1, patch4, BoundarySide2D.U1, BoundaryDirection.Forward)

    return [patch1_ptr, patch2_ptr, patch3_ptr, patch4_ptr]

### Create a list of Frenet frame for specific points on a curve. The Frenet frame is stored as a transformation matrix.
### zvec is a reference vector to compute B at the first sampling point. It shall not be parallel with the tangent vector of the first sampling point.
def GenerateLocalFrenetFrameAt(curve, xi_list, zvec = [1.0, 0.0, 0.0]):
    trans_list = []
    B = Array3()
    ctrl_pnt_grid_func = curve.GridFunction(CONTROL_POINT_COORDINATES)
    # print(ctrl_pnt_grid_func)
    cnt = 0
    for xi in xi_list:
        pnt = [xi, 0.0, 0.0]
        P = ctrl_pnt_grid_func.GetValue(pnt)
        T = ctrl_pnt_grid_func.GetDerivative(pnt)
        T = normalize(T[0])

        if cnt == 0:
            cross(B, zvec, T)
            B = normalize(B)
        else:
            B = B - dot(B, T)*T
            B = normalize(B)

        trans = Transformation(B, T, P)
        trans_list.append(trans)

        cnt = cnt + 1

    return trans_list

### Create a list of Frenet frame along a curve. The Frenet frame is stored as a transformation matrix.
### zvec is a reference vector to compute B at the first sampling point. It shall not be parallel with the tangent vector of the first sampling point.
def GenerateLocalFrenetFrame(curve, num_sampling_points, zvec = [1.0, 0.0, 0.0]):
    xi_list = []
    for i in range(0, num_sampling_points):
        xi_list.append(float(i) / (num_sampling_points-1))
    return GenerateLocalFrenetFrameAt(curve, xi_list, zvec = zvec)

def ExportLocalFrenetFrameToMatlab(trans_list, fn, s = 1.0):
    ifile = open(fn, "w")
    cnt = 1
    ifile.write("s = " + str(s) + ";\n")
    ifile.write("C = {}; B = {}; T = {}; N = {};\n")
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

## Create multipatch of two slabs in 3D along the x-direction. The dimension of each cube along x-, y-, and z-direction is L, w, h correspondingly.
## On output returns the multipatch and the corresponding parameter direction of L, w and h side of the first patch
def CreateTwoSlabs(setting=1, configuration = 1, dist = 0.1, L = 1.0, w = 1.0, h = 1.0):
    if setting == 1:
        return CreateTwoSlabs1(configuration = configuration, dist = dist, L = L, w = w, h = h)
    elif setting == 2:
        return CreateTwoSlabs2(configuration = configuration, dist = dist, L = L, w = w, h = h)
    elif setting == 3:
        return CreateTwoSlabs3(configuration = configuration, dist = dist, L = L, w = w, h = h)
    elif setting == 4:
        return CreateTwoSlabs4(configuration = configuration, dist = dist, L = L, w = w, h = h)
    elif setting == 5:
        return CreateTwoSlabs5(configuration = configuration, dist = dist, L = L, w = w, h = h)
    elif setting == 6:
        return CreateTwoSlabs6(configuration = configuration, dist = dist, L = L, w = w, h = h)
    else:
        raise Exception("Unknown setting %d", setting)

## Create multipatch of two slabs in 3D along the x-direction. The dimension of each cube along x-, y-, and z-direction is L, w, h correspondingly.
## in this setting, the two surface patch matches on v->v and w->w, that means uv_or_vu == True
def CreateTwoSlabs1(configuration = 1, dist = 0.1, L = 1.0, w = 1.0, h = 1.0):

    # create patch 1
    p1 = [0.0, 0.0, 0.0]
    p2 = [L, w, h]
    patch1_ptr = CreateSlab(p1, p2)
    patch1 = patch1_ptr.GetReference()
    patch1.Id = 1
    patch1.LayerIndex = 1

    # create patch 2
    p3 = [L + dist, 0.0, 0.0]
    p4 = [2.0*L + dist, w, h]
    patch2_ptr = CreateSlab(p3, p4)
    patch2 = patch2_ptr.GetReference()
    patch2.Id = 2
    patch2.LayerIndex = 1

    if configuration == 1:

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.U1, patch2, BoundarySide3D.U0, True, BoundaryDirection.Forward, BoundaryDirection.Forward)

    elif configuration == 2:

        bsplines_patch_util.Reverse(patch1, 0)
        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.U0, patch2, BoundarySide3D.U0, True, BoundaryDirection.Forward, BoundaryDirection.Reversed)

    elif configuration == 3:

        bsplines_patch_util.Reverse(patch1, 0)
        bsplines_patch_util.Reverse(patch1, 1)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.U0, patch2, BoundarySide3D.U0, True, BoundaryDirection.Reversed, BoundaryDirection.Forward)

    elif configuration == 4:

        bsplines_patch_util.Reverse(patch1, 1)
        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.U1, patch2, BoundarySide3D.U0, True, BoundaryDirection.Reversed, BoundaryDirection.Reversed)

    else:
        raise Exception("Unknown configuration %d" % (configuration))

    return mpatch, [0, 1, 2]

## Create multipatch of two slabs in 3D along the x-direction. The dimension of each cube along x-, y-, and z-direction is L, w, h correspondingly.
## in this setting, the two surface patch matches on v->w and w->v, that means uv_or_vu == False
def CreateTwoSlabs2(configuration = 1, dist = 0.1, L = 1.0, w = 1.0, h = 1.0):

    # create patch 1
    p1 = [0.0, 0.0, 0.0]
    d1 = [L, 0.0, 0.0]
    d2 = [0.0, 0.0, h]
    d3 = [0.0, w, 0.0]
    patch1_ptr = CreateParallelepiped(p1, d1, d2, d3)
    patch1 = patch1_ptr.GetReference()
    patch1.Id = 1
    patch1.LayerIndex = 1

    # create patch 2
    p3 = [L + dist, 0.0, 0.0]
    p4 = [2.0*L + dist, w, h]
    patch2_ptr = CreateSlab(p3, p4)
    patch2 = patch2_ptr.GetReference()
    patch2.Id = 2
    patch2.LayerIndex = 1

    if configuration == 1:

        bsplines_patch_util.Reverse(patch1, 0)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.U0, patch2, BoundarySide3D.U0, False, BoundaryDirection.Forward, BoundaryDirection.Forward)

    elif configuration == 2:

        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.U1, patch2, BoundarySide3D.U0, False, BoundaryDirection.Forward, BoundaryDirection.Reversed)

    elif configuration == 3:

        bsplines_patch_util.Reverse(patch1, 1)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.U1, patch2, BoundarySide3D.U0, False, BoundaryDirection.Reversed, BoundaryDirection.Forward)

    elif configuration == 4:

        bsplines_patch_util.Reverse(patch1, 0)
        bsplines_patch_util.Reverse(patch1, 1)
        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.U0, patch2, BoundarySide3D.U0, False, BoundaryDirection.Reversed, BoundaryDirection.Reversed)

    else:
        raise Exception("Unknown configuration %d" % (configuration))

    return mpatch, [0, 2, 1]

## Create multipatch of two slabs in 3D along the x-direction. The dimension of each cube along x-, y-, and z-direction is L, w, h correspondingly.
## in this setting, the two surface patch matches on u->v and w->w, that means uv_or_vu == True
def CreateTwoSlabs3(configuration = 1, dist = 0.1, L = 1.0, w = 1.0, h = 1.0):

    # create patch 1
    p1 = [0.0, 0.0, 0.0]
    d1 = [0.0, w, 0.0]
    d2 = [L, 0.0, 0.0]
    d3 = [0.0, 0.0, h]
    patch1_ptr = CreateParallelepiped(p1, d1, d2, d3)
    patch1 = patch1_ptr.GetReference()
    patch1.Id = 1
    patch1.LayerIndex = 1

    # create patch 2
    p3 = [L + dist, 0.0, 0.0]
    p4 = [2.0*L + dist, w, h]
    patch2_ptr = CreateSlab(p3, p4)
    patch2 = patch2_ptr.GetReference()
    patch2.Id = 2
    patch2.LayerIndex = 1

    if configuration == 1:

        bsplines_patch_util.Reverse(patch1, 1)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.V0, patch2, BoundarySide3D.U0, True, BoundaryDirection.Forward, BoundaryDirection.Forward)

    elif configuration == 2:

        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.V1, patch2, BoundarySide3D.U0, True, BoundaryDirection.Forward, BoundaryDirection.Reversed)

    elif configuration == 3:

        bsplines_patch_util.Reverse(patch1, 0)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.V1, patch2, BoundarySide3D.U0, True, BoundaryDirection.Reversed, BoundaryDirection.Forward)

    elif configuration == 4:

        bsplines_patch_util.Reverse(patch1, 0)
        bsplines_patch_util.Reverse(patch1, 1)
        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.V0, patch2, BoundarySide3D.U0, True, BoundaryDirection.Reversed, BoundaryDirection.Reversed)

    else:
        raise Exception("Unknown configuration %d" % (configuration))

    return mpatch, [1, 0, 2]

## Create multipatch of two slabs in 3D along the x-direction. The dimension of each cube along x-, y-, and z-direction is L, w, h correspondingly.
## in this setting, the two surface patch matches on w->v and u->w, that means uv_or_vu == False
def CreateTwoSlabs4(configuration = 1, dist = 0.1, L = 1.0, w = 1.0, h = 1.0):

    # create patch 1
    p1 = [0.0, 0.0, 0.0]
    d1 = [0.0, 0.0, h]
    d2 = [L, 0.0, 0.0]
    d3 = [0.0, w, 0.0]
    patch1_ptr = CreateParallelepiped(p1, d1, d2, d3)
    patch1 = patch1_ptr.GetReference()
    patch1.Id = 1
    patch1.LayerIndex = 1

    # create patch 2
    p3 = [L + dist, 0.0, 0.0]
    p4 = [2.0*L + dist, w, h]
    patch2_ptr = CreateSlab(p3, p4)
    patch2 = patch2_ptr.GetReference()
    patch2.Id = 2
    patch2.LayerIndex = 1

    if configuration == 1:

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.V1, patch2, BoundarySide3D.U0, False, BoundaryDirection.Forward, BoundaryDirection.Forward)

    elif configuration == 2:

        bsplines_patch_util.Reverse(patch1, 2)
        bsplines_patch_util.Reverse(patch1, 1)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.V0, patch2, BoundarySide3D.U0, False, BoundaryDirection.Forward, BoundaryDirection.Reversed)

    elif configuration == 3:

        bsplines_patch_util.Reverse(patch1, 0)
        bsplines_patch_util.Reverse(patch1, 1)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.V0, patch2, BoundarySide3D.U0, False, BoundaryDirection.Reversed, BoundaryDirection.Forward)

    elif configuration == 4:

        bsplines_patch_util.Reverse(patch1, 0)
        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.V1, patch2, BoundarySide3D.U0, False, BoundaryDirection.Reversed, BoundaryDirection.Reversed)

    else:
        raise Exception("Unknown configuration %d" % (configuration))

    return mpatch, [1, 2, 0]

## Create multipatch of two slabs in 3D along the x-direction. The dimension of each cube along x-, y-, and z-direction is L, w, h correspondingly.
## in this setting, the two surface patch matches on u->v and v->w, that means uv_or_vu == True
def CreateTwoSlabs5(configuration = 1, dist = 0.1, L = 1.0, w = 1.0, h = 1.0):

    # create patch 1
    p1 = [0.0, 0.0, 0.0]
    d1 = [0.0, w, 0.0]
    d2 = [0.0, 0.0, h]
    d3 = [L, 0.0, 0.0]
    patch1_ptr = CreateParallelepiped(p1, d1, d2, d3)
    patch1 = patch1_ptr.GetReference()
    patch1.Id = 1
    patch1.LayerIndex = 1

    # create patch 2
    p3 = [L + dist, 0.0, 0.0]
    p4 = [2.0*L + dist, w, h]
    patch2_ptr = CreateSlab(p3, p4)
    patch2 = patch2_ptr.GetReference()
    patch2.Id = 2
    patch2.LayerIndex = 1

    if configuration == 1:

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.W1, patch2, BoundarySide3D.U0, True, BoundaryDirection.Forward, BoundaryDirection.Forward)

    elif configuration == 2:

        bsplines_patch_util.Reverse(patch1, 1)
        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.W0, patch2, BoundarySide3D.U0, True, BoundaryDirection.Forward, BoundaryDirection.Reversed)

    elif configuration == 3:

        bsplines_patch_util.Reverse(patch1, 0)
        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.W0, patch2, BoundarySide3D.U0, True, BoundaryDirection.Reversed, BoundaryDirection.Forward)

    elif configuration == 4:

        bsplines_patch_util.Reverse(patch1, 0)
        bsplines_patch_util.Reverse(patch1, 1)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.W1, patch2, BoundarySide3D.U0, True, BoundaryDirection.Reversed, BoundaryDirection.Reversed)

    else:
        raise Exception("Unknown configuration %d" % (configuration))

    return mpatch, [2, 0, 1]

## Create multipatch of two slabs in 3D along the x-direction. The dimension of each cube along x-, y-, and z-direction is L, w, h correspondingly.
## in this setting, the two surface patch matches on v->v and u->w, that means uv_or_vu == False
def CreateTwoSlabs6(configuration = 1, dist = 0.1, L = 1.0, w = 1.0, h = 1.0):

    # create patch 1
    p1 = [0.0, 0.0, 0.0]
    d1 = [0.0, 0.0, h]
    d2 = [0.0, w, 0.0]
    d3 = [L, 0.0, 0.0]
    patch1_ptr = CreateParallelepiped(p1, d1, d2, d3)
    patch1 = patch1_ptr.GetReference()
    patch1.Id = 1
    patch1.LayerIndex = 1

    # create patch 2
    p3 = [L + dist, 0.0, 0.0]
    p4 = [2.0*L + dist, w, h]
    patch2_ptr = CreateSlab(p3, p4)
    patch2 = patch2_ptr.GetReference()
    patch2.Id = 2
    patch2.LayerIndex = 1

    if configuration == 1:

        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.W0, patch2, BoundarySide3D.U0, False, BoundaryDirection.Forward, BoundaryDirection.Forward)

    elif configuration == 2:

        bsplines_patch_util.Reverse(patch1, 1)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.W1, patch2, BoundarySide3D.U0, False, BoundaryDirection.Forward, BoundaryDirection.Reversed)

    elif configuration == 3:

        bsplines_patch_util.Reverse(patch1, 0)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.W1, patch2, BoundarySide3D.U0, False, BoundaryDirection.Reversed, BoundaryDirection.Forward)

    elif configuration == 4:

        bsplines_patch_util.Reverse(patch1, 0)
        bsplines_patch_util.Reverse(patch1, 1)
        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.W0, patch2, BoundarySide3D.U0, False, BoundaryDirection.Reversed, BoundaryDirection.Reversed)

    else:
        raise Exception("Unknown configuration %d" % (configuration))

    return mpatch, [2, 1, 0]
