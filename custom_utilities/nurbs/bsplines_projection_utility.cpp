//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 9 Oct 2024 $
//   Revision:            $Revision: 1.0 $
//
//

// System includes

// External includes

// Project includes
#include "bsplines_projection_utility.h"
#include "custom_utilities/nurbs/bsplines_patch_utility.h"

namespace Kratos
{

template<typename TPointType>
int BSplinesProjectionUtility::ComputeNormalProjection(const TPointType& rPoint,
        std::vector<double>& rLocalPoint, TPointType& rGlobalPoint, int& patch_id,
        typename MultiPatch<1>::Pointer pMultiPatch, const int nsampling,
        double TOL, int max_iters,
        int echo_level)
{
    typedef MultiPatch<1> MultiPatchType;
    typedef typename MultiPatchType::patch_ptr_iterator patch_ptr_iterator;

    for (patch_ptr_iterator it = pMultiPatch->Patches().ptr_begin(); it != pMultiPatch->Patches().ptr_end(); ++it)
    {
        // check if the patch contains the inner projection
        const auto control_polygon = BSplinesPatchUtility::ExtractControlPolygon(*it);
        array_1d<double, 3> P1, P2;
        for (auto it = control_polygon.begin(); it != control_polygon.end(); ++it)
        {
            (*it)[0].Copy(P1);
            (*it)[1].Copy(P2);
            // TODO
        }

        // TODO estimate the local coordinates of the projection

        // compute the normal projection
        int error_code = IsogeometricProjectionUtility::ComputeNormalProjection<TPointType, 1>(rPoint, rLocalPoint, rGlobalPoint, *it, TOL, max_iters, echo_level);

        if (error_code == 0)
        {
            patch_id = (*it)->Id();
            return 0;
        }
    }

    KRATOS_ERROR << "TO be implemented"; // TODO to be removed

    patch_id = -1;
    return 1; // can't find the projection point
}

}
