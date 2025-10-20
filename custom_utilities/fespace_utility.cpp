//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 19 Oct 2025 $
//
//

#include "custom_utilities/fespace_utility.h"
#include "custom_utilities/weighted_fespace.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"

namespace Kratos
{

template<int TDim, typename TLocalCoordinateType>
typename FESpace<TDim, TLocalCoordinateType>::Pointer FESpaceUtility<TDim, TLocalCoordinateType>::CreateEmptyFESpace(const std::string& type)
{
    const std::string wof_name = ("WeightedFESpace" + std::to_string(TDim) + "D_over_BSplinesFESpace" + std::to_string(TDim) + "D");

    if (type == BSplinesFESpace<TDim, TLocalCoordinateType>::StaticType())
    {
        return BSplinesFESpace<TDim, TLocalCoordinateType>::Create();
    }
    else if (type == WeightedFESpace<TDim, TLocalCoordinateType, double>::StaticType())
    {
        return WeightedFESpace<TDim, TLocalCoordinateType, double>::Create(FESpace<TDim, TLocalCoordinateType>::Create(), {});
    }
    else if (type == wof_name) // hack
    {
        return WeightedFESpace<TDim, TLocalCoordinateType, double>::Create(BSplinesFESpace<TDim, TLocalCoordinateType>::Create(), {});
    }
    else
        return FESpace<TDim, TLocalCoordinateType>::Create();
}

/* template class instantiation */
template class FESpaceUtility<0, double>;
template class FESpaceUtility<1, double>;
template class FESpaceUtility<2, double>;
template class FESpaceUtility<3, double>;

} // end namespace Kratos
