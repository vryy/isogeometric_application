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
#include "custom_utilities/level_set/multipatch_z_level_set.h"
#include "custom_utilities/multipatch.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

void IsogeometricApplication_AddBRepAndLevelSetToPython()
{
    /**************************************************************/
    /************* EXPORT INTERFACE FOR LEVEL SET *****************/
    /**************************************************************/

    class_<MultiPatchZLevelSet, MultiPatchZLevelSet::Pointer, bases<LevelSet, IsogeometricEcho>, boost::noncopyable>
    ( "MultiPatchZLevelSet", init<typename MultiPatch<2>::Pointer>() )
    .def("SetPredictionSampling", &MultiPatchZLevelSet::SetPredictionSampling)
    .def("SetProjectionTolerance", &MultiPatchZLevelSet::SetProjectionTolerance)
    .def("SetMaxIterations", &MultiPatchZLevelSet::SetMaxIterations)
    .def(self_ns::str(self))
    ;
}

}  // namespace Python.

}  // namespace Kratos.
