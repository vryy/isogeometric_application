//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 19 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_INDEXING_UTILITY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_INDEXING_UTILITY_H_INCLUDED

// System includes
#include <vector>
#include <algorithm>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/iga_define.h"

namespace Kratos
{

struct BSplinesIndexingUtility_Helper
{
    /// Compute the BSplines index in 1D.
    static inline std::size_t Index1D(std::size_t i, std::size_t n1)
    {
        return i - 1;
    }

    /// Compute the BSplines index in 2D
    static inline std::size_t Index2D(std::size_t i, std::size_t j,
                                      std::size_t n1, std::size_t n2)
    {
        return (j - 1) * n1 + (i - 1);
    }

    /// Compute the BSplines index in 3D
    static inline std::size_t Index3D(std::size_t i, std::size_t j, std::size_t k,
                                      std::size_t n1, std::size_t n2, std::size_t n3)
    {
        return ((k - 1) * n2 + (j - 1)) * n1 + (i - 1);
    }

    /// Compute the BSplines index array in 1D
    static inline std::vector<std::size_t> IndexArray1D(std::size_t i,
            std::size_t n1)
    {
        return std::vector<std::size_t> {(i - 1) % n1 + 1};
    }

    /// Compute the BSplines index array in 2D
    static inline std::vector<std::size_t> IndexArray2D(std::size_t i,
            std::size_t n1, std::size_t n2)
    {
        return std::vector<std::size_t> {(i - 1) % n1 + 1, (i - 1) / n1 + 1};
    }

    /// Compute the BSplines index array in 3D
    static inline std::vector<std::size_t> IndexArray3D(std::size_t i,
            std::size_t n1, std::size_t n2, std::size_t n3)
    {
        return std::vector<std::size_t> {(i - 1) % n1 + 1, ((i - 1) / n1) % n2 + 1, ((std::size_t)((i - 1) / n1)) / n2 + 1};
    }
};

template<int TDim, class TContainerType, class TIndexContainerType>
struct BSplinesIndexingUtility_Reverse_Helper
{
    static void Reverse(TContainerType& values, const TIndexContainerType& sizes, std::size_t idir)
    {
        return;
    }
};

template<int TDim, class TContainerType, class TIndexContainerType>
struct BSplinesIndexingUtility_Transpose_Helper
{
    static void Transpose(TContainerType& values, const TIndexContainerType& sizes, std::size_t idir, std::size_t jdir)
    {
        return;
    }
};

/**
This class provides sub-routines to index the BSplines basis function in 1D, 2D, 3D.
The base index is assumed to be 1.
 */
class BSplinesIndexingUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BSplinesIndexingUtility);

    /// Type definition

    /// Default constructor
    BSplinesIndexingUtility() {}

    /// Destructor
    virtual ~BSplinesIndexingUtility() {}

    /// Generic function to compute the BSplines index
    template<int TDim>
    static std::size_t Index(const std::vector<std::size_t>& I, const std::vector<std::size_t>& N)
    {
        return -1;
    }

    /// Generic function to compute the BSplines index array
    template<int TDim>
    static std::vector<std::size_t> IndexArray(std::size_t I, const std::vector<std::size_t>& N)
    {
        return std::vector<std::size_t> {};
    }

    /// Reverse an array in specific direction
    /// The array is a {dim}D-tensor which is addressed by using IndexArray
    template<int TDim, class TContainerType, class TIndexContainerType>
    static void Reverse(TContainerType& values, const TIndexContainerType& sizes, std::size_t idir)
    {
        BSplinesIndexingUtility_Reverse_Helper<TDim, TContainerType, TIndexContainerType>::Reverse(values, sizes, idir);
    }

    /// Transform the indices with parameter mapping for 2D surface patch
    static void Transform(std::vector<std::size_t>& func_indices, const std::vector<std::size_t>& size_info,
            const bool uv_or_vu, const BoundaryDirection dir1, const BoundaryDirection dir2);

    /// Transpose an array in specific direction
    /// The array is a {dim}D-tensor which is addressed by using IndexArray
    template<int TDim, class TContainerType, class TIndexContainerType>
    static void Transpose(TContainerType& values, const TIndexContainerType& sizes, std::size_t idir, std::size_t jdir)
    {
        BSplinesIndexingUtility_Transpose_Helper<TDim, TContainerType, TIndexContainerType>::Transpose(values, sizes, idir, jdir);
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "BSplinesIndexingUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const BSplinesIndexingUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

template<class TContainerType, class TIndexContainerType>
struct BSplinesIndexingUtility_Reverse_Helper<1, TContainerType, TIndexContainerType>
{
    static void Reverse(TContainerType& values, const TIndexContainerType& sizes, std::size_t idir)
    {
        if (idir == 0)
        {
            std::reverse(values.begin(), values.end());
        }
        else
            KRATOS_ERROR << "Invalid direction " << idir;
    }
};

template<class TContainerType, class TIndexContainerType>
struct BSplinesIndexingUtility_Reverse_Helper<2, TContainerType, TIndexContainerType>
{
    static void Reverse(TContainerType& values, const TIndexContainerType& sizes, std::size_t idir)
    {
        if (idir == 0)
        {
            for (std::size_t j = 0; j < sizes[1]; ++j)
            {
                std::reverse(values.begin() + j * sizes[0], values.begin() + (j + 1)*sizes[0]);
            }
        }
        else if (idir == 1)
        {
            std::size_t loc;

            for (std::size_t i = 0; i < sizes[0]; ++i)
            {
                // Note: here we touch and reverse each row

                // extract the value
                TContainerType Temp(sizes[1]);
                for (std::size_t j = 0; j < sizes[1]; ++j)
                {
                    loc = BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, sizes[0], sizes[1]);
                    Temp[j] = values[loc];
                }

                // assign the reverse value
                for (std::size_t j = 0; j < sizes[1]; ++j)
                {
                    loc = BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, sizes[0], sizes[1]);
                    values[loc] = Temp[sizes[1] - 1 - j];
                }
            }
        }
        else
            KRATOS_ERROR << "Invalid direction " << idir;
    }
};

template<class TContainerType, class TIndexContainerType>
struct BSplinesIndexingUtility_Reverse_Helper<3, TContainerType, TIndexContainerType>
{
    static void Reverse(TContainerType& values, const TIndexContainerType& sizes, std::size_t idir)
    {
        if (idir == 0)
        {
            for (std::size_t j = 0; j < sizes[1]; ++j)
            {
                for (std::size_t k = 0; k < sizes[2]; ++k)
                {
                    std::reverse(values.begin() + (k * sizes[1] + j)*sizes[0], values.begin() + (k * sizes[1] + j + 1)*sizes[0]);
                }
            }
        }
        else if (idir == 1)
        {
            std::size_t loc;

            for (std::size_t i = 0; i < sizes[0]; ++i)
            {
                for (std::size_t k = 0; k < sizes[2]; ++k)
                {
                    // Note: here we touch and reverse each row

                    // extract the value
                    TContainerType Temp(sizes[1]);
                    for (std::size_t j = 0; j < sizes[1]; ++j)
                    {
                        loc = BSplinesIndexingUtility_Helper::Index3D(i + 1, j + 1, k + 1, sizes[0], sizes[1], sizes[2]);
                        Temp[j] = values[loc];
                    }

                    // assign the reverse value
                    for (std::size_t j = 0; j < sizes[1]; ++j)
                    {
                        loc = BSplinesIndexingUtility_Helper::Index3D(i + 1, j + 1, k + 1, sizes[0], sizes[1], sizes[2]);
                        values[loc] = Temp[sizes[1] - 1 - j];
                    }
                }
            }
        }
        else if (idir == 2)
        {
            std::size_t loc;

            for (std::size_t i = 0; i < sizes[0]; ++i)
            {
                for (std::size_t j = 0; j < sizes[1]; ++j)
                {
                    // Note: here we touch and reverse each row

                    // extract the value
                    TContainerType Temp(sizes[2]);
                    for (std::size_t k = 0; k < sizes[2]; ++k)
                    {
                        loc = BSplinesIndexingUtility_Helper::Index3D(i + 1, j + 1, k + 1, sizes[0], sizes[1], sizes[2]);
                        Temp[k] = values[loc];
                    }

                    // assign the reverse value
                    for (std::size_t k = 0; k < sizes[2]; ++k)
                    {
                        loc = BSplinesIndexingUtility_Helper::Index3D(i + 1, j + 1, k + 1, sizes[0], sizes[1], sizes[2]);
                        values[loc] = Temp[sizes[2] - 1 - k];
                    }
                }
            }
        }
        else
            KRATOS_ERROR << "Invalid direction " << idir;
    }
};

template<class TContainerType, class TIndexContainerType>
struct BSplinesIndexingUtility_Transpose_Helper<1, TContainerType, TIndexContainerType>
{
    static void Transpose(TContainerType& values, const TIndexContainerType& sizes, std::size_t idir, std::size_t jdir)
    {
        KRATOS_ERROR << "Transpose is not relevant for 1D";
    }
};

template<class TContainerType, class TIndexContainerType>
struct BSplinesIndexingUtility_Transpose_Helper<2, TContainerType, TIndexContainerType>
{
    static void Transpose(TContainerType& values, const TIndexContainerType& sizes, std::size_t idir, std::size_t jdir)
    {
        if (idir == jdir)
            return; // DO NOTHING

        if ((idir != 0) && (idir != 1))
            KRATOS_ERROR << "Invalid direction " << idir;

        if ((jdir != 0) && (jdir != 1))
            KRATOS_ERROR << "Invalid direction " << jdir;

        // copy the transpose alues to a new array
        TContainerType new_values(sizes[0]*sizes[1]);
        std::size_t loc, new_loc;

        for (std::size_t i = 0; i < sizes[0]; ++i)
        {
            for (std::size_t j = 0; j < sizes[1]; ++j)
            {
                loc = BSplinesIndexingUtility_Helper::Index2D(i + 1, j + 1, sizes[0], sizes[1]);
                new_loc = BSplinesIndexingUtility_Helper::Index2D(j + 1, i + 1, sizes[1], sizes[0]);

                new_values[new_loc] = values[loc];
            }
        }

        // overwrite the old array with new array
        for (std::size_t k = 0; k < new_values.size(); ++k)
        {
            values[k] = new_values[k];
        }
    }
};

template<class TContainerType, class TIndexContainerType>
struct BSplinesIndexingUtility_Transpose_Helper<3, TContainerType, TIndexContainerType>
{
    static void Transpose(TContainerType& values, const TIndexContainerType& sizes, std::size_t idir, std::size_t jdir)
    {
        if (idir == jdir)
            return; // DO NOTHING

        // TODO
        KRATOS_ERROR << "Transpose is not implemented for 3D";
    }
};

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_INDEXING_UTILITY_H_INCLUDED defined
