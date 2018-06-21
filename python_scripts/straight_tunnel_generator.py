import math
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *

import geometry_factory

### Straight tunnel generator along x-coordinates
### After the model is created, the default order is p_u=2, p_v=2, p_w=1. All the knot vector only contains 0 and 1.
class StraightTunnelGenerator():
    def __init__(self, params):
        self.center = params['center']
        self.round_length = params['round_length']
        self.number_of_slices = params['number_of_slices']
        self.exc_radius = params['exc_radius']
        self.lining_inner_radius = params['lining_inner_radius']
        self.axis = 'x'

        self.generate_excavation = params["generate_excavation"]

        ## dimension to determine the ground
#        self.d1 = params['dimension_d1']
#        self.d2 = params['dimension_d2']
#        self.d3 = params['dimension_d3']
#        self.d4 = params['dimension_d4']
#        self.w1 = params['dimension_w1']
#        self.w2 = params['dimension_w2']
#        self.s = params['exc_inner_cut_ratio']

        self.x_min = -0.5*self.round_length*self.number_of_slices
        self.x_max = -self.x_min

        self.bsplines_patch_util = BSplinesPatchUtility()
        self.patch_id_stepping = 1000
        self.layer_sets = {}

    def CreateModel(self):
        mp = MultiPatch3D()

        number_of_patches_per_slice = 0

        if self.generate_excavation:
            self.layer_sets["excavation"] = []

            for step in range(0, self.number_of_slices):
                x_min = self.x_min + step*self.round_length
                x_max = x_min + self.round_length

                c1 = [x_min, self.center[1], self.center[2]]
                c2 = [x_max, self.center[1], self.center[2]]

                bf1_patches_ptr = self.CreateExcavationBaseSurface(c1, self.exc_radius)
                bf2_patches_ptr = self.CreateExcavationBaseSurface(c2, self.exc_radius)

                if step == 0:
                    number_of_patches_per_slice = number_of_patches_per_slice + len(bf1_patches_ptr)

                for i in range(0, len(bf1_patches_ptr)):
                    bf1_patch_ptr = bf1_patches_ptr[i]
                    bf1_patch = bf1_patch_ptr.GetReference()
                    bf2_patch_ptr = bf2_patches_ptr[i]
                    patch_ptr = self.bsplines_patch_util.CreateConnectedPatch(bf2_patch_ptr, bf1_patch_ptr)
                    patch = patch_ptr.GetReference()
                    patch.Id = step*self.patch_id_stepping + bf1_patch.Id
                    mp.AddPatch(patch_ptr)
                    self.layer_sets["excavation"].append(patch.Id)

                ## make the interface for patches within a slice
                patch1_ptr = mp[step*self.patch_id_stepping + 1]
                patch1 = patch1_ptr.GetReference()
                patch2_ptr = mp[step*self.patch_id_stepping + 2]
                patch2 = patch2_ptr.GetReference()
                patch3_ptr = mp[step*self.patch_id_stepping + 3]
                patch3 = patch3_ptr.GetReference()
                patch4_ptr = mp[step*self.patch_id_stepping + 4]
                patch4 = patch4_ptr.GetReference()
                patch5_ptr = mp[step*self.patch_id_stepping + 5]
                patch5 = patch5_ptr.GetReference()
                patch6_ptr = mp[step*self.patch_id_stepping + 6]
                patch6 = patch6_ptr.GetReference()

                ## TODO check the correctness of the interface

                self.bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.V0, patch3, BoundarySide3D.V1, True, BoundaryDirection.Forward, BoundaryDirection.Forward)
                self.bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.V1, patch2, BoundarySide3D.V0, True, BoundaryDirection.Forward, BoundaryDirection.Forward)

                self.bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.U1, patch4, BoundarySide3D.V1, True, BoundaryDirection.Forward, BoundaryDirection.Forward)
                self.bsplines_patch_util.MakeInterface(patch3, BoundarySide3D.U1, patch4, BoundarySide3D.U0, True, BoundaryDirection.Forward, BoundaryDirection.Forward)

                self.bsplines_patch_util.MakeInterface(patch2, BoundarySide3D.U1, patch5, BoundarySide3D.V1, True, BoundaryDirection.Forward, BoundaryDirection.Forward)
                self.bsplines_patch_util.MakeInterface(patch4, BoundarySide3D.U1, patch5, BoundarySide3D.U0, True, BoundaryDirection.Forward, BoundaryDirection.Forward)

                self.bsplines_patch_util.MakeInterface(patch2, BoundarySide3D.V1, patch6, BoundarySide3D.V1, True, BoundaryDirection.Reversed, BoundaryDirection.Forward)
                self.bsplines_patch_util.MakeInterface(patch5, BoundarySide3D.U1, patch6, BoundarySide3D.U0, True, BoundaryDirection.Forward, BoundaryDirection.Forward)

        ## make the interfaces for patches along x-direction
        for step in range(0, self.number_of_slices-1):
            for i in range(0, number_of_patches_per_slice):
                patch1_ptr = mp[step*self.patch_id_stepping + i+1]
                patch1 = patch1_ptr.GetReference()
                patch2_ptr = mp[(step+1)*self.patch_id_stepping + i+1]
                patch2 = patch2_ptr.GetReference()
                self.bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.W0, patch2, BoundarySide3D.W1, True, BoundaryDirection.Forward, BoundaryDirection.Forward)

        return mp

    ## Create patches in the excavation layer of the left side
    def CreateExcavationBaseSurface(self, center, exc_radius):
        multipatch_refine_util = MultiPatchRefinementUtility()

        ####### create patch 3

        ### create arc 1
        arc1_ptr = geometry_factory.CreateSmallArc(center, 'x', exc_radius, 90.0, 45.0)
        arc1 = arc1_ptr.GetReference()
#        print("arc1:")
#        print(arc1)

        ### create line 1
        p1 = [center[0], center[1], center[2] + 0.5*exc_radius]
        p2 = [center[0], center[1] + 0.5*exc_radius, center[2] + 0.5*exc_radius]
        line1_ptr = geometry_factory.CreateLine(p1, p2, arc1.Order(0))
        line1 = line1_ptr.GetReference()
#        print("line1:")
#        print(line1)

        patch3_ptr = self.bsplines_patch_util.CreateConnectedPatch(arc1, line1)
        multipatch_refine_util.DegreeElevate(patch3_ptr, [0, 1])
        patch3 = patch3_ptr.GetReference()
        patch3.Id = 3

        ######## create patch 1

        ### create line 2
        p1 = [center[0], center[1], center[2]]
        p2 = [center[0], center[1] + 0.5*exc_radius, center[2]]
        line2_ptr = geometry_factory.CreateLine(p1, p2, arc1.Order(0))
        line2 = line2_ptr.GetReference()
#        print("line2")
#        print(line2)

        patch1_ptr = self.bsplines_patch_util.CreateConnectedPatch(line1, line2)
        multipatch_refine_util.DegreeElevate(patch1_ptr, [0, 1])
        patch1 = patch1_ptr.GetReference()
        patch1.Id = 1

        ######## create patch 2

        ### create line 3
        p1 = [center[0], center[1], center[2] - 0.5*exc_radius]
        p2 = [center[0], center[1] + 0.5*exc_radius, center[2] - 0.5*exc_radius]
        line3_ptr = geometry_factory.CreateLine(p1, p2, arc1.Order(0))
        line3 = line3_ptr.GetReference()
#        print("line3")
#        print(line3)

        patch2_ptr = self.bsplines_patch_util.CreateConnectedPatch(line2, line3)
        multipatch_refine_util.DegreeElevate(patch2_ptr, [0, 1])
        patch2 = patch2_ptr.GetReference()
        patch2.Id = 2

        ######## create patch 4

        ### create line 4
        p1 = [center[0], center[1] + 0.5*exc_radius, center[2] + 0.5*exc_radius]
        p2 = [center[0], center[1] + 0.5*exc_radius, center[2]]
        line4_ptr = geometry_factory.CreateLine(p1, p2, arc1.Order(0))
        line4 = line4_ptr.GetReference()
#        print("line4")
#        print(line4)

        ### create arc 2
        arc2_ptr = geometry_factory.CreateSmallArc(center, 'x', exc_radius, 45.0, 0.0)
        arc2 = arc2_ptr.GetReference()
#        print("arc2:")
#        print(arc2)

        patch4_ptr = self.bsplines_patch_util.CreateConnectedPatch(arc2, line4)
        multipatch_refine_util.DegreeElevate(patch4_ptr, [0, 1])
        patch4 = patch4_ptr.GetReference()
        patch4.Id = 4

        ######## create patch 5

        ### create line 5
        p1 = [center[0], center[1] + 0.5*exc_radius, center[2]]
        p2 = [center[0], center[1] + 0.5*exc_radius, center[2] - 0.5*self.exc_radius]
        line5_ptr = geometry_factory.CreateLine(p1, p2, arc1.Order(0))
        line5 = line5_ptr.GetReference()
#        print("line5")
#        print(line5)

        ### create arc 3
        arc3_ptr = geometry_factory.CreateSmallArc(center, 'x', exc_radius, 0.0, -45.0)
        arc3 = arc3_ptr.GetReference()
#        print("arc3:")
#        print(arc3)

        patch5_ptr = self.bsplines_patch_util.CreateConnectedPatch(arc3, line5)
        multipatch_refine_util.DegreeElevate(patch5_ptr, [0, 1])
        patch5 = patch5_ptr.GetReference()
        patch5.Id = 5

        ######## create patch 6

        ### create line 6
        p1 = [center[0], center[1] + 0.5*exc_radius, center[2] - 0.5*self.exc_radius]
        p2 = [center[0], center[1], center[2] - 0.5*self.exc_radius]
        line6_ptr = geometry_factory.CreateLine(p1, p2, arc1.Order(0))
        line6 = line6_ptr.GetReference()
#        print("line6")
#        print(line6)

        ### create arc 4
        arc4_ptr = geometry_factory.CreateSmallArc(center, 'x', exc_radius, -45.0, -90.0)
        arc4 = arc4_ptr.GetReference()
#        print("arc4:")
#        print(arc4)

        patch6_ptr = self.bsplines_patch_util.CreateConnectedPatch(arc4, line6)
        multipatch_refine_util.DegreeElevate(patch6_ptr, [0, 1])
        patch6 = patch6_ptr.GetReference()
        patch6.Id = 6

        return [patch1_ptr, patch2_ptr, patch3_ptr, patch4_ptr, patch5_ptr, patch6_ptr]





