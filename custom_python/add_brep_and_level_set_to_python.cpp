// see brep_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//

// External includes
#include <boost/python.hpp>

// Project includes
#include "custom_python/add_brep_and_level_set_to_python.h"
#include "custom_utilities/curve/patch_curve.h"
#include "custom_utilities/level_set/multipatch_z_level_set.h"
#include "custom_utilities/level_set/multipatch_normal_level_set.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

static PatchCurve<2>::Pointer PatchCurveOnSurface_init(Patch<2>::Pointer pPatch,
        const int idir, const double xi)
{
    auto pCurve = PatchCurve<2>::Pointer(new PatchCurve<2>(pPatch));
    pCurve->SetLocalCoord(0, idir, xi);
    return pCurve;
}

static PatchCurve<3>::Pointer PatchCurveOnVolume_init(Patch<3>::Pointer pPatch,
        const int idir1, const int idir2, const double xi1, const double xi2)
{
    auto pCurve = PatchCurve<3>::Pointer(new PatchCurve<3>(pPatch));
    pCurve->SetLocalCoord(0, idir1, xi1);
    pCurve->SetLocalCoord(1, idir2, xi2);
    return pCurve;
}

void IsogeometricApplication_AddBRepAndLevelSetToPython()
{
    /**************************************************************/
    /*************** EXPORT INTERFACE FOR CURVE *******************/
    /**************************************************************/

    class_<PatchCurve<1>, PatchCurve<1>::Pointer, bases<Curve>, boost::noncopyable>
    ( "PatchCurve", init<Patch<1>::Pointer>() )
    ;

    class_<PatchCurve<2>, PatchCurve<2>::Pointer, bases<Curve>, boost::noncopyable>
    ( "PatchCurveOnSurface", no_init )
    .def("__init__", make_constructor(&PatchCurveOnSurface_init))
    ;

    class_<PatchCurve<3>, PatchCurve<3>::Pointer, bases<Curve>, boost::noncopyable>
    ( "PatchCurveOnVolume", no_init )
    .def("__init__", make_constructor(&PatchCurveOnVolume_init))
    ;

    /**************************************************************/
    /************* EXPORT INTERFACE FOR LEVEL SET *****************/
    /**************************************************************/

    class_<MultiPatchZLevelSet, MultiPatchZLevelSet::Pointer, bases<LevelSet, IsogeometricEcho>, boost::noncopyable>
    ( "MultiPatchZLevelSet", init<typename MultiPatch<2>::Pointer>() )
    .def("SetPredictionSampling", &MultiPatchZLevelSet::SetPredictionSampling)
    .def("SetProjectionTolerance", &MultiPatchZLevelSet::SetProjectionTolerance)
    .def("SetMaxIterations", &MultiPatchZLevelSet::SetMaxIterations)
    ;

    class_<MultiPatchNormalLevelSet<1>, MultiPatchNormalLevelSet<1>::Pointer, bases<LevelSet, IsogeometricEcho>, boost::noncopyable>
    ( "MultiPatchNormalLevelSet1D", init<typename MultiPatch<1>::Pointer>() )
    .def("SetPredictionSampling", &MultiPatchNormalLevelSet<1>::SetPredictionSampling)
    .def("SetProjectionTolerance", &MultiPatchNormalLevelSet<1>::SetProjectionTolerance)
    .def("SetMaxIterations", &MultiPatchNormalLevelSet<1>::SetMaxIterations)
    .def("SetInnerPoint", &MultiPatchNormalLevelSet<1>::SetInnerPoint<array_1d<double, 3> >)
    ;

    class_<MultiPatchNormalLevelSet<2>, MultiPatchNormalLevelSet<2>::Pointer, bases<LevelSet, IsogeometricEcho>, boost::noncopyable>
    ( "MultiPatchNormalLevelSet2D", init<typename MultiPatch<2>::Pointer>() )
    .def("SetPredictionSampling", &MultiPatchNormalLevelSet<2>::SetPredictionSampling)
    .def("SetProjectionTolerance", &MultiPatchNormalLevelSet<2>::SetProjectionTolerance)
    .def("SetMaxIterations", &MultiPatchNormalLevelSet<2>::SetMaxIterations)
    .def("SetCurve", &MultiPatchNormalLevelSet<2>::SetCurve)
    ;
}

}  // namespace Python.

}  // namespace Kratos.
