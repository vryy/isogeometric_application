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
#include "custom_utilities/multipatch.h"
#include "custom_utilities/patch_interface.h"
#include "custom_utilities/nurbs/bsplines_patch_interface.h"

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

    /// Create the B-Splines multipatch from geo file
    /// This function is kept for backward compatibility. New user should use MultiNURBSPatchGeoImporter instead.
    template<int TDim>
    static typename MultiPatch<TDim>::Pointer CreateMultiPatchFromGeo(const std::string& fn)
    {
        MultiNURBSPatchGeoImporter<TDim> dummy;
        return dummy.Import(fn);
    }

    /// Dummy function to silence the compiler
    static void MakeInterface2D(typename Patch<1>::Pointer pPatch1, const BoundarySide& side1,
            typename Patch<1>::Pointer pPatch2, const BoundarySide& side2, const BoundaryDirection& direction)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not reallistic in 1D")
    }

    /// Make the interface between two patches in 2D
    static void MakeInterface2D(typename Patch<2>::Pointer pPatch1, const BoundarySide& side1,
            typename Patch<2>::Pointer pPatch2, const BoundarySide& side2, const BoundaryDirection& direction)
    {
        typename FESpace<1>::Pointer pBFESpace1 = pPatch1->pFESpace()->ConstructBoundaryFESpace(side1);

        typename FESpace<1>::Pointer pBFESpace2;

        std::map<std::size_t, std::size_t> local_parameter_map;
        std::vector<BoundaryDirection> directions = {direction};
        pBFESpace2 = pPatch2->pFESpace()->ConstructBoundaryFESpace(side2, local_parameter_map, directions);

        if( (*pBFESpace1) == (*pBFESpace2) )
        {
            typename PatchInterface<2>::Pointer pInterface12;
            typename PatchInterface<2>::Pointer pInterface21;

            pInterface12 = boost::make_shared<BSplinesPatchInterface<2> >(pPatch1, side1, pPatch2, side2, direction);
            pInterface21 = boost::make_shared<BSplinesPatchInterface<2> >(pPatch2, side2, pPatch1, side1, direction);

            pInterface12->SetOtherInterface(pInterface21);
            pInterface21->SetOtherInterface(pInterface12);

            pPatch1->AddInterface(pInterface12);
            pPatch2->AddInterface(pInterface21);
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "The interface is not created because the two patch's boundaries are not conformed.", "")
    }

    /// Dummy function to silence the compiler
    static void MakeInterface2D(typename Patch<3>::Pointer pPatch1, const BoundarySide& side1,
            typename Patch<3>::Pointer pPatch2, const BoundarySide& side2, const BoundaryDirection& direction)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not reallistic in 3D")
    }

    /// Dummy function to silence the compiler
    static void MakeInterface3D(typename Patch<1>::Pointer pPatch1, const BoundarySide& side1,
            typename Patch<1>::Pointer pPatch2, const BoundarySide& side2, const bool& uv_or_vu,
            const BoundaryDirection& direction1, const BoundaryDirection& direction2)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not reallistic in 1D")
    }

    /// Dummy function to silence the compiler
    static void MakeInterface3D(typename Patch<2>::Pointer pPatch1, const BoundarySide& side1,
            typename Patch<2>::Pointer pPatch2, const BoundarySide& side2, const bool& uv_or_vu,
            const BoundaryDirection& direction1, const BoundaryDirection& direction2)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not reallistic in 2D")
    }

    /// Make the interface between two patches in 3D
    static void MakeInterface3D(typename Patch<3>::Pointer pPatch1, const BoundarySide& side1,
            typename Patch<3>::Pointer pPatch2, const BoundarySide& side2, const bool& uv_or_vu,
            const BoundaryDirection& direction1, const BoundaryDirection& direction2)
    {
        typename FESpace<2>::Pointer pBFESpace1 = pPatch1->pFESpace()->ConstructBoundaryFESpace(side1);

        typename FESpace<2>::Pointer pBFESpace2;

        std::map<std::size_t, std::size_t> local_parameter_map;
        if (uv_or_vu)
        {
            local_parameter_map[0] = 0;
            local_parameter_map[1] = 1;
        }
        else
        {
            local_parameter_map[0] = 1;
            local_parameter_map[1] = 0;
        }
        std::vector<BoundaryDirection> directions = {direction1, direction2};
        pBFESpace2 = pPatch2->pFESpace()->ConstructBoundaryFESpace(side2, local_parameter_map, directions);

        if( (*pBFESpace1) == (*pBFESpace2) )
        {
            typename PatchInterface<3>::Pointer pInterface12;
            typename PatchInterface<3>::Pointer pInterface21;

            pInterface12 = boost::make_shared<BSplinesPatchInterface<3> >(pPatch1, side1, pPatch2, side2, uv_or_vu, direction1, direction2);
            pInterface21 = boost::make_shared<BSplinesPatchInterface<3> >(pPatch2, side2, pPatch1, side1, uv_or_vu, direction1, direction2);

            pInterface12->SetOtherInterface(pInterface21);
            pInterface21->SetOtherInterface(pInterface12);

            pPatch1->AddInterface(pInterface12);
            pPatch2->AddInterface(pInterface21);
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "The interface is not created because the two patch's boundaries are not conformed.", "")
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

