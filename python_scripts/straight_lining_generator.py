import math
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *

import geometry_factory

### Straight lining generator along x-coordinates
class StraightLiningGenerator():
    def __init__(self, params = None):
        self.center = params['center']
        self.x_thick = params['round_length']
        self.rout = params['exc_radius']
        self.rin = params['lining_inner_radius']
        self.axis = 'x'

        self.bsplines_patch_util = BSplinesPatchUtility()

    ### Create one segment at coordinate x and with start angle and end angle
    def CreateSegment(self, x, sangle, eangle):
        cen1 = self.center
        cen1[0] = x
        ring1_ptr = geometry_factory.CreateSmallRing(cen1, self.axis, self.rin, self.rout, sangle, eangle)
        ring1 = ring1_ptr.GetReference()

        cen2 = self.center
        cen2[0] = x + self.x_thick
        ring2_ptr = geometry_factory.CreateSmallRing(cen2, self.axis, self.rin, self.rout, sangle, eangle)
        ring2 = ring2_ptr.GetReference()

        ## create segment patch by connecting the two rings
        segment_patch_ptr = self.bsplines_patch_util.CreateConnectedPatch(ring1, ring2)
        return segment_patch_ptr

    ### Create a complete ring at coordinate x. The number of segment is nsegment + 1 keystone
    ### The angle of the keystone is ksangle and keangle. The angle of the other segment will be determined by ksangle, keangle and nsegment
    def CreateRing(self, x, ksangle, keangle, nsegment):
        ring = []
        keystone_ptr = self.CreateSegment(x, ksangle, keangle)
        ring.append(keystone_ptr)
        sweep = (360.0 - abs(keangle-ksangle)) / nsegment
        if keangle > ksangle:
            dir = 1.0
        else:
            dir = -1.0
        for i in range(0, nsegment):
            sangle = keangle + dir*i*sweep
            eangle = sangle + dir*sweep
            segment_ptr = self.CreateSegment(x, sangle, eangle)
            ring.append(segment_ptr)
        return ring


