//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 26 Apr 2018 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_BENDING_STRIP_UTILITY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_BENDING_STRIP_UTILITY_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/patch_interface.h"
#include "custom_utilities/nurbs/bending_strip_nurbs_patch.h"

namespace Kratos
{

/**
Utility class to for bending strip patch.
 */
class BendingStripUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BendingStripUtility);

    /// Type definition
    typedef ControlPoint<double> ControlPointType;

    /// Default constructor
    BendingStripUtility() {}

    /// Destructor
    virtual ~BendingStripUtility() {}

    template<int TDim>
    static typename Patch<TDim>::Pointer CreateBendingStripNURBSPatch(
        const std::size_t& Id,
        typename Patch<TDim>::Pointer pPatch1, const BoundarySide& side1,
        typename Patch<TDim>::Pointer pPatch2, const BoundarySide& side2,
        const int& Order)
    {
        return typename Patch<TDim>::Pointer( new BendingStripNURBSPatch<TDim>(Id, pPatch1, side1, pPatch2, side2, Order) );
    }

    template<int TDim>
    static typename Patch<TDim>::Pointer CreateBendingStripNURBSPatch(
        const std::size_t& Id,
        typename Patch<TDim>::Pointer pPatch1, const BoundarySide& side1,
        typename Patch<TDim>::Pointer pPatch2, const BoundarySide& side2,
        const std::vector<int>& Orders)
    {
        return typename Patch<TDim>::Pointer( new BendingStripNURBSPatch<TDim>(Id, pPatch1, side1, pPatch2, side2, Orders) );
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "BendingStripUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const BendingStripUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_BENDING_STRIP_UTILITY_H_INCLUDED defined

