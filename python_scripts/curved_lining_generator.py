import math
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *

import geometry_factory

mpatch_export3 = MultiNURBSPatchMatlabExporter()
bsplines_patch_util = BSplinesPatchUtility()

### Straight lining generator along x-coordinates
class CurvedLiningGenerator():
    def __init__(self, params):
        self.bsplines_patch_util = BSplinesPatchUtility()
        self.params = params

    ### Reconstruct the segment based on given data
    ### REF: Bui et al, BIM-based Modelling and Simulation of Segmental Tunnel Lining by means of Isogeometric Analysis
    def CreateSegment(self, segment_data, initial_trans = None):
        ## segment data
        rin = segment_data['inner_radii']
        rout = segment_data['outer_radii']
        start_angle_front = segment_data['start_angle_front']
        angle_front = segment_data['angle_front']
        start_angle_back = segment_data['start_angle_back']
        angle_back = segment_data['angle_back']
        cpoints = segment_data['point']

        ## create the front face
        front_face_ptr = geometry_factory.CreateSmallRing([0.0, 0.0, 0.0], 'z', rin, rout, start_angle_front, start_angle_front+angle_front)
        front_face = front_face_ptr.GetReference()

        ## generate sweep curve
        order = self.params['sweep_curve_order']
        curve_ptr = geometry_factory.CreateCurve(cpoints, order)
        curve = curve_ptr.GetReference()
        curve.Prefix = "curve"
#        mpatch_export3.Export(curve, "curve.m")

        ## generate the local Frenet frame along the trajectory
        npoints = self.params["number_of_sweep_sampling"]
        trans_list = geometry_factory.GenerateLocalFrenetFrame(curve, npoints)

        if initial_trans == None:
            new_trans_list = trans_list
        else:
            # transform the first frame to initial_trans and subsequently other trans to ensure continuity
            new_trans_list = []
            for i in range(0, len(trans_list)):
                trans = Transformation()
                if i == 0:
                    trans.AppendTransformation(initial_trans)
                else:
                    trans.AppendTransformation(trans_list[i])
                    trans.AppendTransformation(trans_list[0].Inverse())
                    trans.AppendTransformation(initial_trans)
                new_trans_list.append(trans)

#        for trans in new_trans_list:
#            print(trans)
#        geometry_factory.ExportLocalFrenetFrameToMatlab(new_trans_list, "frame.m", 2e-1)

        ## generate a list of cut section
        sections_ptr = []
        for i in range(0, npoints):
            start_angle = start_angle_front + i*(start_angle_back-start_angle_front)/(npoints-1)
            end_angle = start_angle_front + angle_front + i*(start_angle_back+angle_back-start_angle_front-angle_front)/(npoints-1)
            section_ptr = geometry_factory.CreateSmallRing([0.0, 0.0, 0.0], 'z', rin, rout, start_angle, end_angle)
            section = section_ptr.GetReference()
            section.Id = i+1
            section.ApplyTransformation(new_trans_list[i])
            # mpatch_export3.Export(section, section.Name() + "_def.m")
            sections_ptr.append(section_ptr)

        ## create segment patch by connecting the two rings
        ## create the sweep volume
        order_w = self.params['sweep_order']
        segment_ptr = bsplines_patch_util.CreateConnectedPatch(sections_ptr, order_w)
        return [segment_ptr, new_trans_list]


