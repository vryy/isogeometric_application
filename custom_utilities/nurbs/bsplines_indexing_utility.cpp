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
std::vector<std::size_t> BSplinesIndexingUtility::IndexArray<1>(std::size_t I, const std::vector<std::size_t>& N)
{
    return BSplinesIndexingUtility_Helper::IndexArray1D(I, N[0]);
}

template<>
std::vector<std::size_t> BSplinesIndexingUtility::IndexArray<2>(std::size_t I, const std::vector<std::size_t>& N)
{
    return BSplinesIndexingUtility_Helper::IndexArray2D(I, N[0], N[1]);
}

template<>
std::vector<std::size_t> BSplinesIndexingUtility::IndexArray<3>(std::size_t I, const std::vector<std::size_t>& N)
{
    return BSplinesIndexingUtility_Helper::IndexArray3D(I, N[0], N[1], N[2]);
}

void BSplinesIndexingUtility::Transform(std::vector<std::size_t>& func_indices, const std::vector<std::size_t>& size_info,
        const bool uv_or_vu, const BoundaryDirection dir1, const BoundaryDirection dir2)
{
    if (!uv_or_vu)
    {
        // transpose the function indices
        std::vector<std::size_t> new_func_indices(func_indices.size());

        for (std::size_t i = 0; i < size_info[0]; ++i)
        {
            for (std::size_t j = 0; j < size_info[1]; ++j)
            {
                const auto k = BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, size_info[0], size_info[1]);
                const auto kt = BSplinesIndexingUtility_Helper::Index2D(j + 1, i + 1, size_info[1], size_info[0]);
                new_func_indices[kt] = func_indices[k];
            }
        }

        func_indices = new_func_indices;
    }

    // reverse in each direction if needed
    if (dir2 == BoundaryDirection::_REVERSED_)
    {
        for (std::size_t i = 0; i < size_info[0]; ++i)
        {
            std::vector<std::size_t> tmp(size_info[1]);
            for (std::size_t j = 0; j < size_info[1]; ++j)
            {
                tmp[j] = func_indices[BSplinesIndexingUtility_Helper::Index2D(i + 1, size_info[1] - j, size_info[0], size_info[1])];
            }
            for (std::size_t j = 0; j < size_info[1]; ++j)
            {
                func_indices[BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, size_info[0], size_info[1])] = tmp[j];
            }
        }
    }

    if (dir1 == BoundaryDirection::_REVERSED_)
    {
        for (std::size_t j = 0; j < size_info[1]; ++j)
        {
            std::vector<std::size_t> tmp(size_info[0]);
            for (std::size_t i = 0; i < size_info[0]; ++i)
            {
                tmp[i] = func_indices[BSplinesIndexingUtility_Helper::Index2D(size_info[0] - i, j + 1, size_info[0], size_info[1])];
            }
            for (std::size_t i = 0; i < size_info[0]; ++i)
            {
                func_indices[BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, size_info[0], size_info[1])] = tmp[i];
            }
        }
    }
}

} // namespace Kratos.
