//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 13 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_UTILITY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_UTILITY_H_INCLUDED

// System includes
#include <vector>
#include <fstream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_utilities/isogeometric_utility.h"
#include "custom_utilities/grid_function.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/patch_interface.h"
#include "custom_utilities/multipatch.h"

namespace Kratos
{

/**
This class is a library to generate typical NURBS patch for computational mechanics benchmarks.
 */
class MultiPatchUtility : public IsogeometricUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiPatchUtility);

    /// Type definition

    /// Default constructor
    MultiPatchUtility() {}

    /// Destructor
    virtual ~MultiPatchUtility() {}

    /// Create new patch from a FESpace and wrap it with pointer
    #ifdef SD_APP_FORWARD_COMPATIBILITY
    template<int TDim>
    static iga::Wrapper<Patch<TDim>, typename Patch<TDim>::Pointer> CreatePatchPointer(const std::size_t& Id, typename FESpace<TDim>::Pointer pFESpace)
    {
        return iga::Wrapper<Patch<TDim>, typename Patch<TDim>::Pointer>(typename Patch<TDim>::Pointer(new Patch<TDim>(Id, pFESpace)));
    }
    #else
    template<int TDim>
    static typename Patch<TDim>::Pointer CreatePatchPointer(const std::size_t& Id, typename FESpace<TDim>::Pointer pFESpace)
    {
        return typename Patch<TDim>::Pointer(new Patch<TDim>(Id, pFESpace));
    }
    #endif

    /// Make a simple interface between two patches
    /// For BSplines patch, one shall use BSplinesPatchUtility::MakeInterfacexD instead
    template<int TDim>
    static void MakeInterface(typename Patch<TDim>::Pointer pPatch1, const BoundarySide& side1,
            typename Patch<TDim>::Pointer pPatch2, const BoundarySide& side2)
    {
        typename PatchInterface<TDim>::Pointer pInterface12;
        typename PatchInterface<TDim>::Pointer pInterface21;

        pInterface12 = iga::make_shared<PatchInterface<TDim> >(pPatch1, side1, pPatch2, side2);
        pInterface21 = iga::make_shared<PatchInterface<TDim> >(pPatch2, side2, pPatch1, side1);

        pInterface12->SetOtherInterface(pInterface21);
        pInterface21->SetOtherInterface(pInterface12);

        pPatch1->AddInterface(pInterface12);
        pPatch2->AddInterface(pInterface21);
    }

    /// Construct the 12 edge patches of a 3D patch
    static std::vector<Patch<1>::Pointer> ConstructEdgePatches(Patch<3>::ConstPointer pPatch)
    {
        Patch<2>::Pointer pTopPatch = pPatch->ConstructBoundaryPatch(_BTOP_);
        Patch<2>::Pointer pBottomPatch = pPatch->ConstructBoundaryPatch(_BBOTTOM_);
        Patch<2>::Pointer pLeftPatch = pPatch->ConstructBoundaryPatch(_BLEFT_);
        Patch<2>::Pointer pRightPatch = pPatch->ConstructBoundaryPatch(_BRIGHT_);

        std::vector<Patch<1>::Pointer> pEdgePatches;

        pEdgePatches.push_back(pBottomPatch->ConstructBoundaryPatch(_BLEFT_));
        pEdgePatches.push_back(pBottomPatch->ConstructBoundaryPatch(_BRIGHT_));
        pEdgePatches.push_back(pBottomPatch->ConstructBoundaryPatch(_BBOTTOM_));
        pEdgePatches.push_back(pBottomPatch->ConstructBoundaryPatch(_BTOP_));

        pEdgePatches.push_back(pTopPatch->ConstructBoundaryPatch(_BLEFT_));
        pEdgePatches.push_back(pTopPatch->ConstructBoundaryPatch(_BRIGHT_));
        pEdgePatches.push_back(pTopPatch->ConstructBoundaryPatch(_BBOTTOM_));
        pEdgePatches.push_back(pTopPatch->ConstructBoundaryPatch(_BTOP_));

        pEdgePatches.push_back(pLeftPatch->ConstructBoundaryPatch(_BLEFT_));
        pEdgePatches.push_back(pLeftPatch->ConstructBoundaryPatch(_BRIGHT_));
        pEdgePatches.push_back(pRightPatch->ConstructBoundaryPatch(_BLEFT_));
        pEdgePatches.push_back(pRightPatch->ConstructBoundaryPatch(_BRIGHT_));

        return pEdgePatches;
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MultiPatchUtility";
    }

    void PrintData(std::ostream& rOStream) const override
    {
    }

}; // end class MultiPatchUtility


/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const MultiPatchUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_UTILITY_H_INCLUDED defined

