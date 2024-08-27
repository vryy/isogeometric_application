//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_PATCH_UTILITY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_PATCH_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/multipatch.h"
#include "custom_utilities/patch_interface.h"

namespace Kratos
{

/**
THis class supports some operations on B-Splines patch
 */
class BSplinesPatchUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BSplinesPatchUtility);

    /// Default constructor
    BSplinesPatchUtility() {}

    /// Destructor
    virtual ~BSplinesPatchUtility() {}

    /// Construct a higher dimension patch by connecting two patches with the straight B-Splines curve. The order of the connection curve is 1.
    /// To have higher order of the connection one needs to elevate the degree.
    /// Right now, the two sub-patches must have same parameters (knot vectors) and are B-Splines.
    template<int TDim>
    static typename Patch<TDim>::Pointer CreateLoftPatch(typename Patch < TDim - 1 >::Pointer pPatch1, typename Patch < TDim - 1 >::Pointer pPatch2);

    /// Construct a higher dimension patch by connecting multiple patches with the B-Splines curve.
    /// The knot vector of the curve is adjusted accordingly with the given order and the number of patches.
    /// The knot vector will be by default uniform.
    /// To have higher order of the connection one needs to elevate the degree.
    /// Right now, the two sub-patches must have same parameters (knot vectors) and are B-Splines.
    template<int TDim>
    static typename Patch<TDim>::Pointer CreateLoftPatch(std::vector < typename Patch < TDim - 1 >::Pointer > pPatches, int order);

    /// Reverse the B-Splines patch in specific direction
    template<int TDim>
    static void Reverse(typename Patch<TDim>::Pointer pPatch, std::size_t idir);

    /// Reverse the B-Splines patch in specific direction
    template<int TDim>
    static void ReverseImpl(typename Patch<TDim>::Pointer pPatch, std::size_t idir, std::set<std::size_t>& reversed_patches);

    /// Get the dimension of underlying NURBS in geo file
    static int GetDimensionOfGeo(const std::string& fn);

    /// Create the B-Splines patch from geo file
    /// This function is kept for backward compatibility. New user should use MultiNURBSPatchGeoImporter instead.
    template<int TDim>
    static typename Patch<TDim>::Pointer CreatePatchFromGeo(const std::string& fn);

    /// Create the B-Splines multipatch from geo file
    /// This function is kept for backward compatibility. New user should use MultiNURBSPatchGeoImporter instead.
    template<int TDim>
    static typename MultiPatch<TDim>::Pointer CreateMultiPatchFromGeo(const std::string& fn);

    /// Dummy function to silence the compiler
    static void MakeInterface2D(typename Patch<1>::Pointer pPatch1, const BoundarySide side1,
                                typename Patch<1>::Pointer pPatch2, const BoundarySide side2, const BoundaryDirection direction);

    /// Make the interface between two patches in 2D
    static void MakeInterface2D(typename Patch<2>::Pointer pPatch1, const BoundarySide side1,
                                typename Patch<2>::Pointer pPatch2, const BoundarySide side2, const BoundaryDirection direction);

    /// Dummy function to silence the compiler
    static void MakeInterface2D(typename Patch<3>::Pointer pPatch1, const BoundarySide side1,
                                typename Patch<3>::Pointer pPatch2, const BoundarySide side2, const BoundaryDirection direction);

    /// Dummy function to silence the compiler
    static void MakeInterface3D(typename Patch<1>::Pointer pPatch1, const BoundarySide side1,
                                typename Patch<1>::Pointer pPatch2, const BoundarySide side2, const bool uv_or_vu,
                                const BoundaryDirection direction1, const BoundaryDirection direction2);

    /// Dummy function to silence the compiler
    static void MakeInterface3D(typename Patch<2>::Pointer pPatch1, const BoundarySide side1,
                                typename Patch<2>::Pointer pPatch2, const BoundarySide side2, const bool uv_or_vu,
                                const BoundaryDirection direction1, const BoundaryDirection direction2);

    /// Make the interface between two patches in 3D
    static void MakeInterface3D(typename Patch<3>::Pointer pPatch1, const BoundarySide side1,
                                typename Patch<3>::Pointer pPatch2, const BoundarySide side2, const bool uv_or_vu,
                                const BoundaryDirection direction1, const BoundaryDirection direction2);

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "BSplinesPatchUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const BSplinesPatchUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_BSPLINES_INDEXING_UTILITY_H_INCLUDED defined
