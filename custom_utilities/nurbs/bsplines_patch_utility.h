//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Nov 2017 $
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
class KRATOS_API(ISOGEOMETRIC_APPLICATION) BSplinesPatchUtility
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
    template<class TPatchType>
    static typename TPatchType::Pointer CreateLoftPatch(typename TPatchType::BoundaryPatchType::Pointer pPatch1, typename TPatchType::BoundaryPatchType::Pointer pPatch2);

    /// Construct a higher dimension patch by connecting multiple patches with the B-Splines curve.
    /// The knot vector of the curve is adjusted accordingly with the given order and the number of patches.
    /// The knot vector will be by default uniform.
    /// To have higher order of the connection one needs to elevate the degree.
    /// Right now, the sub-patches must have same parameters (knot vectors) and are B-Splines.
    /// Since the connecting curve is B-Splines, it is generally assumed that the sub-patches are not the cut section of the lofting volume,
    /// except the case that the connecting curve is a line
    template<class TPatchType>
    static typename TPatchType::Pointer CreateLoftPatch(std::vector<typename TPatchType::BoundaryPatchType::Pointer> pPatches, int order);

    /// Reverse the B-Splines patch in specific direction
    template<typename TPatchPointerType>
    static void Reverse(TPatchPointerType pPatch, std::size_t idir);

    /// Reverse the B-Splines patch in specific direction
    template<typename TPatchPointerType>
    static void ReverseImpl(TPatchPointerType pPatch, std::size_t idir, std::set<std::size_t>& reversed_patches);

    /// Transpose the 2D B-Splines patch
    template<class TPatchPointerType>
    static void Transpose2D(TPatchPointerType pPatch);

    /// Transpose the 3D B-Splines patch
    template<class TPatchPointerType>
    static void Transpose3D(TPatchPointerType pPatch, std::size_t idir, std::size_t jdir);

    /// Transpose the B-Splines patch
    template<class TPatchType>
    static void Transpose(TPatchType& rPatch, std::size_t idir, std::size_t jdir);

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

    /// Make the interface between two patches in 1D
    template<typename TPatchPointerType>
    static void MakeInterface1D(TPatchPointerType pPatch1, const BoundarySide side1,
                                TPatchPointerType pPatch2, const BoundarySide side2);

    /// Make the interface between two patches in 2D
    template<typename TPatchPointerType>
    static void MakeInterface2D(TPatchPointerType pPatch1, const BoundarySide side1,
                                TPatchPointerType pPatch2, const BoundarySide side2, const BoundaryDirection direction);

    /// Make the interface between two patches in 3D
    template<typename TPatchPointerType>
    static void MakeInterface3D(TPatchPointerType pPatch1, const BoundarySide side1,
                                TPatchPointerType pPatch2, const BoundarySide side2, const bool uv_or_vu,
                                const BoundaryDirection direction1, const BoundaryDirection direction2);

    /// Search and create interfaces automatically for multipatch
    template<class TMultiPatchType>
    static void CreateInterfaces(typename TMultiPatchType::Pointer pMultiPatch);

    /// Check for knot duplication for all the patches. This shall be used when we want no reduction in continuity.
    template<class TMultiPatchType>
    static void CheckRepeatedKnot(typename TMultiPatchType::Pointer pMultiPatch);

    /// Extract the control polygon of a 1D patch
    static std::vector<std::array<typename Patch<1>::ControlPointType, 2> > ExtractControlPolygon(typename Patch<1>::ConstPointer pPatch);

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
