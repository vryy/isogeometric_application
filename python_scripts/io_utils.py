""" Create uniform knot vector
"""
def create_uniform_knot_vectors(number, order):

    # create the knot vector
    knot_vec = []

    nu = number
    pu = order

    for i in range(0, pu+1):
        knot_vec.append(0.0)
    for i in range(0, nu-pu-1):
        knot_vec.append((i+1)*1.0 / (nu-pu))
    for i in range(0, pu+1):
        knot_vec.append(1.0)

    return knot_vec

""" Write a single 3D surface patch
    points have to be list in reverse lexicographic order: the points increment first
    in u direction, then v, and then w.
"""
def write_surface_patch_geo_v21(filename, knot_vec_u, knot_vec_v, all_points):
    nu = len(all_points[0])
    nv = len(all_points)

    pu = len(knot_vec_u) - nu - 1
    pv = len(knot_vec_v) - nv - 1

    w = 1.0 # unit weight for all points

    ifile = open(filename, "w")

    ifile.write("# nurbs mesh v.2.1\n")
    ifile.write("2 3 1\n")
    ifile.write("PATCH 1\n")

    ifile.write("%d %d\n" % (pu, pv))
    ifile.write("%d %d\n" % (nu, nv))

    for k in knot_vec_u:
        ifile.write(" %.10f" % (k))
    ifile.write("\n")
    for k in knot_vec_v:
        ifile.write(" %.10f" % (k))
    ifile.write("\n")

    for points in all_points:
        for p in points:
            ifile.write(" %.10f" % (p[0]))
    ifile.write("\n")
    for points in all_points:
        for p in points:
            ifile.write(" %.10f" % (p[1]))
    ifile.write("\n")
    for points in all_points:
        for p in points:
            ifile.write(" %.10f" % (p[2]))
    ifile.write("\n")
    for points in all_points:
        for p in points:
            ifile.write(" %.10f" % (w))
    ifile.write("\n")

    ifile.close()
