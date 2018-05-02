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
#include <vector>

// External includes
#include "boost/algorithm/string.hpp"

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"
#include "custom_utilities/nurbs/bsplines_fespace_library.h"
#include "custom_utilities/import_export/multi_nurbs_patch_geo_importer.h"
#include "custom_utilities/patch.h"

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
    static typename Patch<TDim>::Pointer CreateConnectedPatch(typename Patch<TDim-1>::Pointer pPatch1, typename Patch<TDim-1>::Pointer pPatch2)
    {
        // check prerequisites
        if (pPatch1->pFESpace()->Type() != BSplinesFESpace<TDim-1>::StaticType())
            KRATOS_THROW_ERROR(std::logic_error, "Patch 1 is not B-Splines patch", "")

        if (pPatch2->pFESpace()->Type() != BSplinesFESpace<TDim-1>::StaticType())
            KRATOS_THROW_ERROR(std::logic_error, "Patch 2 is not B-Splines patch", "")

        // check the FESpace
        if ( !pPatch1->pFESpace()->IsCompatible(*(pPatch2->pFESpace())) )
        {
            KRATOS_THROW_ERROR(std::logic_error, "The two patches are not compatible", "")
        }

        // create the new FESpace
        typename BSplinesFESpace<TDim-1>::Pointer pFESpace1 = boost::dynamic_pointer_cast<BSplinesFESpace<TDim-1> >(pPatch1->pFESpace());
        typename BSplinesFESpace<TDim>::Pointer pNewFESpace = BSplinesFESpace<TDim>::Create();
        for (std::size_t dim = 0; dim < TDim-1; ++dim)
        {
            pNewFESpace->SetKnotVector(dim, pFESpace1->KnotVector(dim));
            pNewFESpace->SetInfo(dim, pFESpace1->Number(dim), pFESpace1->Order(dim));
        }

        std::size_t connect_order = 1;
        typename BSplinesFESpace<TDim>::knot_container_t new_knot_vector = BSplinesFESpaceLibrary::CreatePrimitiveOpenKnotVector(connect_order);
        pNewFESpace->SetKnotVector(TDim-1, new_knot_vector);
        pNewFESpace->SetInfo(TDim-1, connect_order+1, connect_order);

        // create the new patch
        typename Patch<TDim>::Pointer pNewPatch = typename Patch<TDim>::Pointer(new Patch<TDim>(-1, pNewFESpace));

        // create the new control point grid
        typedef typename Patch<TDim>::ControlPointType ControlPointType;
        typename StructuredControlGrid<TDim-1, ControlPointType>::Pointer pControlPointGrid1
            = boost::dynamic_pointer_cast<StructuredControlGrid<TDim-1, ControlPointType> >(pPatch1->pControlPointGridFunction()->pControlGrid());
        typename ControlGrid<ControlPointType>::Pointer pControlPointGrid2 = pPatch2->pControlPointGridFunction()->pControlGrid();

        //// make a size check, it is not necessary anyway
        if (pControlPointGrid1->size() != pControlPointGrid2->size())
            KRATOS_THROW_ERROR(std::logic_error, "The size of two control point grid are not the same", "")

        // assign data to the new control point grid
        std::vector<std::size_t> new_sizes(TDim);
        for (std::size_t dim = 0; dim < TDim-1; ++dim)
            new_sizes[dim] = pControlPointGrid1->Size(dim);
        new_sizes[TDim-1] = 2;
        typename StructuredControlGrid<TDim, ControlPointType>::Pointer pNewControlPointGrid = StructuredControlGrid<TDim, ControlPointType>::Create(new_sizes);
        for (std::size_t i = 0; i < pControlPointGrid1->size(); ++i)
        {
            pNewControlPointGrid->SetData(i, pControlPointGrid1->GetData(i));
            pNewControlPointGrid->SetData(i + pControlPointGrid1->size(), pControlPointGrid2->GetData(i));
        }
        pNewControlPointGrid->SetName(pControlPointGrid1->Name());

        // assign new control point grid to new patch
        pNewPatch->CreateControlPointGridFunction(pNewControlPointGrid);

        // TODO create other grid function data

        // enumerate the first time
        std::size_t start = 0;
        pNewPatch->pFESpace()->ResetFunctionIndices();
        start = pNewPatch->pFESpace()->Enumerate(start);

        return pNewPatch;
    }

    /// Get the dimension of underlying NURBS in geo file
    static int GetDimensionOfGeo(const std::string& fn)
    {
        return GetDimensionOfGeoHelper(fn);
    }

    /// Create the B-Splines patch from geo file
    /// This function is kept for backward compatibility. New user should use MultiNURBSPatchGeoImporter instead.
    template<int TDim>
    static typename Patch<TDim>::Pointer CreatePatchFromGeo(const std::string& fn)
    {
        MultiNURBSPatchGeoImporter<TDim> dummy;
        return dummy.ImportSingle(fn);
    }

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

