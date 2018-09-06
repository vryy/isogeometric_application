//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 19 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes
#include "custom_utilities/nurbs/bsplines_indexing_utility.h"

namespace Kratos
{

template<>
std::size_t BSplinesIndexingUtility::Index<1>(const std::vector<std::size_t>& I, const std::vector<std::size_t>& N)
{
    return BSplinesIndexingUtility_Helper::Index1D(I[0], N[0]);
}

template<>
std::size_t BSplinesIndexingUtility::Index<2>(const std::vector<std::size_t>& I, const std::vector<std::size_t>& N)
{
    return BSplinesIndexingUtility_Helper::Index2D(I[0], I[1], N[0], N[1]);
}

template<>
std::size_t BSplinesIndexingUtility::Index<3>(const std::vector<std::size_t>& I, const std::vector<std::size_t>& N)
{
    return BSplinesIndexingUtility_Helper::Index3D(I[0], I[1], I[2], N[0], N[1], N[2]);
}

template<>
std::vector<std::size_t> BSplinesIndexingUtility::IndexArray<1>(const std::size_t& I, const std::vector<std::size_t>& N)
{
    return BSplinesIndexingUtility_Helper::IndexArray1D(I, N[0]);
}

template<>
std::vector<std::size_t> BSplinesIndexingUtility::IndexArray<2>(const std::size_t& I, const std::vector<std::size_t>& N)
{
    return BSplinesIndexingUtility_Helper::IndexArray2D(I, N[0], N[1]);
}

template<>
std::vector<std::size_t> BSplinesIndexingUtility::IndexArray<3>(const std::size_t& I, const std::vector<std::size_t>& N)
{
    return BSplinesIndexingUtility_Helper::IndexArray3D(I, N[0], N[1], N[2]);
}

} // namespace Kratos.

