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

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{

struct BSplinesIndexingUtility_Helper
{
    /// Compute the BSplines index in 1D.
    static inline std::size_t Index1D(const std::size_t& i, const std::size_t& n1)
    {
        return i-1;
    }

    /// Compute the BSplines index in 2D
    static inline std::size_t Index2D(const std::size_t& i, const std::size_t& j,
        const std::size_t& n1, const std::size_t& n2)
    {
        return (j-1)*n1 + (i-1);
    }

    /// Compute the BSplines index in 3D
    static inline std::size_t Index3D(const std::size_t& i, const std::size_t& j, const std::size_t& k,
        const std::size_t& n1, const std::size_t& n2, const std::size_t& n3)
    {
        return ((k-1)*n2 + (j-1))*n1 + (i-1);
    }

    /// Compute the BSplines index array in 1D
    static inline std::vector<std::size_t> IndexArray1D(const std::size_t& i,
        const std::size_t& n1)
    {
        return std::vector<std::size_t>{(i-1) % n1 + 1};
    }

    /// Compute the BSplines index array in 2D
    static inline std::vector<std::size_t> IndexArray2D(const std::size_t& i,
        const std::size_t& n1, const std::size_t& n2)
    {
        return std::vector<std::size_t>{(i-1) % n1 + 1, (i-1) / n1 + 1};
    }

    /// Compute the BSplines index array in 3D
    static inline std::vector<std::size_t> IndexArray3D(const std::size_t& i,
        const std::size_t& n1, const std::size_t& n2, const std::size_t& n3)
    {
        return std::vector<std::size_t>{(i-1) % n1 + 1, ((i-1) / n1) % n2 + 1, ((std::size_t)((i-1) / n1)) / n2 + 1};
    }
};

template<int TDim, class TContainerType, class TIndexContainerType>
struct BSplinesIndexingUtility_Reverse_Helper
{
    static void Reverse(TContainerType& values, const TIndexContainerType& sizes, const std::size_t& idir)
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
    static std::vector<std::size_t> IndexArray(const std::size_t& I, const std::vector<std::size_t>& N)
    {
        return std::vector<std::size_t>{};
    }

    /// Reverse an array in specific direction
    /// The array is a {dim}D-tensor which is addressed by using IndexArray
    template<int TDim, class TContainerType, class TIndexContainerType>
    static void Reverse(TContainerType& values, const TIndexContainerType& sizes, const std::size_t& idir)
    {
        BSplinesIndexingUtility_Reverse_Helper<TDim, TContainerType, TIndexContainerType>::Reverse(values, sizes, idir);
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
    static void Reverse(TContainerType& values, const TIndexContainerType& sizes, const std::size_t& idir)
    {
        if (idir == 0)
            std::reverse(values.begin(), values.end());
    }
};

template<class TContainerType, class TIndexContainerType>
struct BSplinesIndexingUtility_Reverse_Helper<2, TContainerType, TIndexContainerType>
{
    static void Reverse(TContainerType& values, const TIndexContainerType& sizes, const std::size_t& idir)
    {
        if (idir == 0)
        {
            for (std::size_t j = 0; j < sizes[1]; ++j)
            {
                std::reverse(values.begin() + j*sizes[0], values.begin() + (j+1)*sizes[0]);
            }
        }
        else if (idir == 1)
        {
            std::size_t loc;

            for (std::size_t i = 0; i < sizes[0]; ++i)
            {
                // extract the value
                TContainerType Temp(sizes[1]);
                for (std::size_t j = 0; j < sizes[1]; ++j)
                {
                    loc = BSplinesIndexingUtility_Helper::Index2D(i+1, j+1, sizes[0], sizes[1]);
                    Temp[j] = values[loc];
                }

                // assign the reverse value
                for (std::size_t j = 0; j < sizes[1]; ++j)
                {
                    loc = BSplinesIndexingUtility_Helper::Index2D(i+1, j+1, sizes[0], sizes[1]);
                    values[loc] = Temp[sizes[1]-1-j];
                }
            }
        }
    }
};

template<class TContainerType, class TIndexContainerType>
struct BSplinesIndexingUtility_Reverse_Helper<3, TContainerType, TIndexContainerType>
{
    static void Reverse(TContainerType& values, const TIndexContainerType& sizes, const std::size_t& idir)
    {
        if (idir == 0)
        {
            for (std::size_t j = 0; j < sizes[1]; ++j)
                for (std::size_t k = 0; k < sizes[2]; ++k)
                    std::reverse(values.begin() + (k*sizes[1] + j)*sizes[0], values.begin() + (k*sizes[1] + j + 1)*sizes[0]);
        }
        else if (idir == 1)
        {
            std::size_t loc;

            for (std::size_t i = 0; i < sizes[0]; ++i)
            {
                for (std::size_t k = 0; k < sizes[2]; ++k)
                {
                    // extract the value
                    TContainerType Temp(sizes[1]);
                    for (std::size_t j = 0; j < sizes[1]; ++j)
                    {
                        loc = BSplinesIndexingUtility_Helper::Index3D(i+1, j+1, k+1, sizes[0], sizes[1], sizes[2]);
                        Temp[j] = values[loc];
                    }

                    // assign the reverse value
                    for (std::size_t j = 0; j < sizes[1]; ++j)
                    {
                        loc = BSplinesIndexingUtility_Helper::Index3D(i+1, j+1, k+1, sizes[0], sizes[1], sizes[2]);
                        values[loc] = Temp[sizes[1]-1-j];
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
                    // extract the value
                    TContainerType Temp(sizes[2]);
                    for (std::size_t k = 0; k < sizes[2]; ++k)
                    {
                        loc = BSplinesIndexingUtility_Helper::Index3D(i+1, j+1, k+1, sizes[0], sizes[1], sizes[2]);
                        Temp[k] = values[loc];
                    }

                    // assign the reverse value
                    for (std::size_t k = 0; k < sizes[2]; ++k)
                    {
                        loc = BSplinesIndexingUtility_Helper::Index3D(i+1, j+1, k+1, sizes[0], sizes[1], sizes[2]);
                        values[loc] = Temp[sizes[2]-1-k];
                    }
                }
            }
        }
    }
};

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_INDEXING_UTILITY_H_INCLUDED defined
