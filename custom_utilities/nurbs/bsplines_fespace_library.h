//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 8 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_FESPACE_LIBRARY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_FESPACE_LIBRARY_H_INCLUDED

// System includes
#include <vector>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"

namespace Kratos
{

/**
This class is a library to generate typical BSplines patch for computational mechanics benchmarks.
 */
class BSplinesFESpaceLibrary
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BSplinesFESpaceLibrary);

    /// Type definition
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /// Default constructor
    BSplinesFESpaceLibrary() {}

    /// Destructor
    virtual ~BSplinesFESpaceLibrary() {}

    /// Create the primitive open knot vector with order p
    /// The primitive knot vector is the knot vector of only 0 and 1.
    static knot_container_t CreatePrimitiveOpenKnotVector(const std::size_t& order)
    {
        knot_container_t knot_vector;
        for (std::size_t i = 0; i < order+1; ++i)
            knot_vector.pCreateKnot(0.0);
        for (std::size_t i = 0; i < order+1; ++i)
            knot_vector.pCreateKnot(1.0);
        return knot_vector;
    }

    /// Create the uniform open knot vector with order p and n basis functions
    static knot_container_t CreateUniformOpenKnotVector(const std::size_t& number, const std::size_t& order)
    {
        knot_container_t knot_vector;
        for (std::size_t i = 0; i < order+1; ++i)
            knot_vector.pCreateKnot(0.0);
        std::size_t m = number - order;
        for (std::size_t i = 0; i < m-1; ++i)
            knot_vector.pCreateKnot(((double)(i+1))/m);
        for (std::size_t i = 0; i < order+1; ++i)
            knot_vector.pCreateKnot(1.0);
        return knot_vector;
    }

    /// Generate regular BSplines patch. For 2D, it's rectangle and for 3D it's a cube.
    /// The knot vector only contains 0 and 1, i.e [0 0 0 ... 1 1 1].
    /// All the weights are 1.
    template<int TDim>
    static typename BSplinesFESpace<TDim>::Pointer CreatePrimitiveFESpace(const std::vector<std::size_t>& Orders)
    {
        typename BSplinesFESpace<TDim>::Pointer pFESpace = typename BSplinesFESpace<TDim>::Pointer(new BSplinesFESpace<TDim>());

        for (std::size_t dim = 0; dim < TDim; ++dim)
        {
            knot_container_t knot_vector = CreatePrimitiveOpenKnotVector(Orders[dim]);
            pFESpace->SetKnotVector(dim, knot_vector);
            pFESpace->SetInfo(dim, Orders[dim]+1, Orders[dim]);
        }

        pFESpace->ResetFunctionIndices();

        return pFESpace;
    }

    /// Generate regular BSplines patch. For 2D, it's rectangle and for 3D it's a cube.
    /// The knot vector is assumed open and uniform.
    /// All the weights are 1.
    template<int TDim>
    static typename BSplinesFESpace<TDim>::Pointer CreateUniformFESpace(const std::vector<std::size_t>& Numbers,
        const std::vector<std::size_t>& Orders)
    {
        typename BSplinesFESpace<TDim>::Pointer pFESpace = typename BSplinesFESpace<TDim>::Pointer(new BSplinesFESpace<TDim>());

        for (std::size_t dim = 0; dim < TDim; ++dim)
        {
            knot_container_t knot_vector = CreateUniformOpenKnotVector(Numbers[dim], Orders[dim]);
            pFESpace->SetKnotVector(dim, knot_vector);
            pFESpace->SetInfo(dim, Numbers[dim], Orders[dim]);
        }

        pFESpace->ResetFunctionIndices();

        return pFESpace;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "BSplinesFESpaceLibrary";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const BSplinesFESpaceLibrary& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_FESPACE_LIBRARY_H_INCLUDED defined

