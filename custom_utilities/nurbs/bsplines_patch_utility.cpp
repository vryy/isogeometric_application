//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_utilities/import_export/multi_nurbs_patch_geo_importer.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"
#include "custom_utilities/nurbs/bsplines_fespace_library.h"
#include "custom_utilities/nurbs/bsplines_patch_interface.h"
#include "custom_utilities/nurbs/bsplines_patch_utility.h"

namespace Kratos
{

template<int TDim>
typename Patch<TDim>::Pointer BSplinesPatchUtility::CreateLoftPatch(typename Patch < TDim - 1 >::Pointer pPatch1, typename Patch < TDim - 1 >::Pointer pPatch2)
{
    std::vector < typename Patch < TDim - 1 >::Pointer > pPatches = {pPatch1, pPatch2};
    return CreateLoftPatch<TDim>(pPatches, 1);
}

template<int TDim>
typename Patch<TDim>::Pointer BSplinesPatchUtility::CreateLoftPatch(std::vector < typename Patch < TDim - 1 >::Pointer > pPatches, int order)
{
    if (pPatches.size() == 0)
    {
        return typename Patch<TDim>::Pointer(new Patch<TDim>(-1));
    }

    // check prerequisites
    for (std::size_t i = 0; i < pPatches.size(); ++i)
    {
        if (pPatches[i]->pFESpace()->Type() != BSplinesFESpace < TDim - 1 >::StaticType())
        {
            KRATOS_ERROR << "Patch " << pPatches[i]->Name() << " is not B-Splines patch";
        }

        // check the FESpace
        if (i != 0)
        {
            if ( !pPatches[0]->pFESpace()->IsCompatible(*(pPatches[i]->pFESpace())) )
            {
                KRATOS_ERROR << "The patch " << pPatches[i]->Name() << " is not compatible with patch 0";
            }
        }
    }

    // create the new FESpace
    typename BSplinesFESpace < TDim - 1 >::Pointer pFESpace0 = iga::dynamic_pointer_cast < BSplinesFESpace < TDim - 1 > > (pPatches[0]->pFESpace());
    if (pFESpace0 == NULL)
        KRATOS_THROW_ERROR(std::runtime_error, "The cast to BSplinesFESpace is failed.", "")
        typename BSplinesFESpace<TDim>::Pointer pNewFESpace = BSplinesFESpace<TDim>::Create();
    for (std::size_t dim = 0; dim < TDim - 1; ++dim)
    {
        pNewFESpace->SetKnotVector(dim, pFESpace0->KnotVector(dim));
        pNewFESpace->SetInfo(dim, pFESpace0->Number(dim), pFESpace0->Order(dim));
    }

    typename BSplinesFESpace<TDim>::knot_container_t new_knot_vector = BSplinesFESpaceLibrary::CreateUniformOpenKnotVector(pPatches.size(), order);
    pNewFESpace->SetKnotVector(TDim - 1, new_knot_vector);
    pNewFESpace->SetInfo(TDim - 1, pPatches.size(), order);

    // create the new patch
    typename Patch<TDim>::Pointer pNewPatch = typename Patch<TDim>::Pointer(new Patch<TDim>(-1, pNewFESpace));

    // create the new control point grid
    typedef typename Patch<TDim>::ControlPointType ControlPointType;
    typename StructuredControlGrid < TDim - 1, ControlPointType >::Pointer pControlPointGrid0
        = iga::dynamic_pointer_cast < StructuredControlGrid < TDim - 1, ControlPointType > > (pPatches[0]->pControlPointGridFunction()->pControlGrid());
    if (pControlPointGrid0 == NULL)
    {
        KRATOS_ERROR << "The cast to StructuredControlGrid is failed.";
    }

    //// make a size check, it is not necessary anyway
    for (std::size_t j = 1; j < pPatches.size(); ++j)
    {
        typename ControlGrid<ControlPointType>::Pointer pControlPointGrid = pPatches[j]->pControlPointGridFunction()->pControlGrid();
        if (pControlPointGrid0->size() != pControlPointGrid->size())
        {
            KRATOS_ERROR << "The size of the control point grid is not the same as the first control grid";
        }
    }

    // assign data to the new control point grid
    std::vector<std::size_t> new_sizes(TDim);
    for (std::size_t dim = 0; dim < TDim - 1; ++dim)
    {
        new_sizes[dim] = pControlPointGrid0->Size(dim);
    }
    new_sizes[TDim - 1] = pPatches.size();
    typename StructuredControlGrid<TDim, ControlPointType>::Pointer pNewControlPointGrid = StructuredControlGrid<TDim, ControlPointType>::Create(new_sizes);
    for (std::size_t j = 0; j < pPatches.size(); ++j)
    {
        typename ControlGrid<ControlPointType>::Pointer pControlPointGrid = pPatches[j]->pControlPointGridFunction()->pControlGrid();
        for (std::size_t i = 0; i < pControlPointGrid0->size(); ++i)
        {
            pNewControlPointGrid->SetData(i + j * pControlPointGrid0->size(), pControlPointGrid->GetData(i));
        }
    }
    pNewControlPointGrid->SetName(pControlPointGrid0->Name());

    // assign new control point grid to new patch
    pNewPatch->CreateControlPointGridFunction(pNewControlPointGrid);

    // TODO create other grid function data

    // reset the function indices
    pNewPatch->pFESpace()->ResetFunctionIndices();

    // // enumerate the first time
    //std::size_t start = 0;
    //start = pNewPatch->pFESpace()->Enumerate(start);

    return pNewPatch;
}

template<int TDim>
void BSplinesPatchUtility::Reverse(typename Patch<TDim>::Pointer pPatch, std::size_t idir)
{
    std::set<std::size_t> reversed_patches;
    ReverseImpl<TDim>(pPatch, idir, reversed_patches);
}

template<int TDim>
void BSplinesPatchUtility::ReverseImpl(typename Patch<TDim>::Pointer pPatch, std::size_t idir, std::set<std::size_t>& reversed_patches)
{
    if (reversed_patches.find(pPatch->Id()) != reversed_patches.end())
    {
        return;
    }

    std::cout << "Patch " << pPatch->Name() << " will be reversed in direction " << idir << std::endl;

    if (pPatch->pFESpace()->Type() != BSplinesFESpace<TDim>::StaticType())
    {
        KRATOS_ERROR << "Patch " << pPatch->Name() << " is not B-Splines patch. Reverse can't be done.";
    }

    // reverse the FESPace
    typename BSplinesFESpace<TDim>::Pointer pFESpace = iga::dynamic_pointer_cast<BSplinesFESpace<TDim> >(pPatch->pFESpace());
    pFESpace->Reverse(idir);

    // reverse the structured control grid
    typedef typename Patch<TDim>::ControlPointType ControlPointType;
    typename StructuredControlGrid<TDim, ControlPointType>::Pointer pControlPointGrid =
        iga::dynamic_pointer_cast<StructuredControlGrid<TDim, ControlPointType> >(pPatch->pControlPointGridFunction()->pControlGrid());
    if (pControlPointGrid != NULL)
    {
        pControlPointGrid->Reverse(idir);
    }
    else
    {
        KRATOS_ERROR << "The control point grid is not structured";
    }

    typedef typename Patch<TDim>::DoubleGridFunctionContainerType DoubleGridFunctionContainerType;
    DoubleGridFunctionContainerType DoubleGridFunctions_ = pPatch->DoubleGridFunctions();
    for (typename DoubleGridFunctionContainerType::const_iterator it = DoubleGridFunctions_.begin();
            it != DoubleGridFunctions_.end(); ++it)
    {
        typename StructuredControlGrid<TDim, double>::Pointer pControlValueGrid =
            iga::dynamic_pointer_cast<StructuredControlGrid<TDim, double> >((*it)->pControlGrid());
        if (pControlValueGrid != NULL)
        {
            pControlValueGrid->Reverse(idir);
        }
        else
        {
            KRATOS_ERROR << "The control value grid " << (*it)->pControlGrid()->Name() << " is not structured";
        }
    }

    typedef typename Patch<TDim>::Array1DGridFunctionContainerType Array1DGridFunctionContainerType;
    Array1DGridFunctionContainerType Array1DGridFunctions_ = pPatch->Array1DGridFunctions();
    for (typename Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions_.begin();
            it != Array1DGridFunctions_.end(); ++it)
    {
        typename StructuredControlGrid<TDim, array_1d<double, 3> >::Pointer pControlValueGrid =
            iga::dynamic_pointer_cast<StructuredControlGrid<TDim, array_1d<double, 3> > >((*it)->pControlGrid());
        if ((*it)->pControlGrid()->Name() == "CONTROL_POINT_COORDINATES") { continue; }
        if (pControlValueGrid != NULL)
        {
            pControlValueGrid->Reverse(idir);
        }
        else
        {
            KRATOS_ERROR << "The control value grid " << (*it)->pControlGrid()->Name() << " is not structured";
        }
    }

    typedef typename Patch<TDim>::VectorGridFunctionContainerType VectorGridFunctionContainerType;
    VectorGridFunctionContainerType VectorGridFunctions_ = pPatch->VectorGridFunctions();
    for (typename VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions_.begin();
            it != VectorGridFunctions_.end(); ++it)
    {
        typename StructuredControlGrid<TDim, Vector>::Pointer pControlValueGrid =
            iga::dynamic_pointer_cast<StructuredControlGrid<TDim, Vector> >((*it)->pControlGrid());
        if (pControlValueGrid != NULL)
        {
            pControlValueGrid->Reverse(idir);
        }
        else
        {
            KRATOS_ERROR << "The control value grid " << (*it)->pControlGrid()->Name() << " is not structured";
        }
    }

    // add to the list of reversed patches
    reversed_patches.insert(pPatch->Id());

    // look for neighbours and reverse
    if (TDim > 1)
    {
        for (std::size_t i = 0; i < pPatch->NumberOfInterfaces(); ++i)
        {
            typename BSplinesPatchInterface<TDim>::Pointer pInterface = iga::dynamic_pointer_cast<BSplinesPatchInterface<TDim> >(pPatch->pInterface(i));

            if (pInterface != NULL)
            {
                if (TDim == 2)
                {
                    std::size_t idir1 = static_cast<std::size_t>(ParameterDirection<2>::Get_(pInterface->Side1()));
                    std::size_t idir2 = static_cast<std::size_t>(ParameterDirection<2>::Get_(pInterface->Side2()));

                    if (idir1 == idir)
                    {
                        // reverse the neighbour patch
                        ReverseImpl<TDim>(pInterface->pPatch2(), idir2, reversed_patches);
                    }
                    else
                    {
                        pInterface->FlipSide1();
                        pInterface->pOtherInterface()->FlipSide2();
                    }
                }
                else if (TDim == 3)
                {
                    std::vector<int> dirs1 = ParameterDirection<3>::Get(pInterface->Side1());
                    std::vector<int> dirs2 = ParameterDirection<3>::Get(pInterface->Side2());

                    if (static_cast<std::size_t>(dirs1[0]) == idir)
                    {
                        std::size_t idir2 = dirs2[pInterface->LocalParameterMapping(0)];

                        // reverse the neighbour patch
                        ReverseImpl<TDim>(pInterface->pPatch2(), idir2, reversed_patches);
                    }
                    else if (static_cast<std::size_t>(dirs1[1]) == idir)
                    {
                        std::size_t idir2 = dirs2[pInterface->LocalParameterMapping(1)];

                        // reverse the neighbour patch
                        ReverseImpl<TDim>(pInterface->pPatch2(), idir2, reversed_patches);
                    }
                    else
                    {
                        pInterface->FlipSide1();
                        pInterface->pOtherInterface()->FlipSide2();
                    }
                }
            }
            else
            {
                KRATOS_ERROR << "The interface is not B-Splines patch interface";
            }
        }
    }
}

int BSplinesPatchUtility::GetDimensionOfGeo(const std::string& fn)
{
    return GetDimensionOfGeoHelper(fn);
}

template<int TDim>
typename Patch<TDim>::Pointer BSplinesPatchUtility::CreatePatchFromGeo(const std::string& fn)
{
    MultiNURBSPatchGeoImporter<TDim> dummy;
    return dummy.ImportSingle(fn);
}

template<int TDim>
typename MultiPatch<TDim>::Pointer BSplinesPatchUtility::CreateMultiPatchFromGeo(const std::string& fn)
{
    MultiNURBSPatchGeoImporter<TDim> dummy;
    return dummy.Import(fn);
}

void BSplinesPatchUtility::MakeInterface1D(typename Patch<1>::Pointer pPatch1, const BoundarySide side1,
        typename Patch<1>::Pointer pPatch2, const BoundarySide side2)
{
    typename FESpace<0>::Pointer pBFESpace1 = pPatch1->pFESpace()->ConstructBoundaryFESpace(side1);

    typename FESpace<0>::Pointer pBFESpace2 = pPatch2->pFESpace()->ConstructBoundaryFESpace(side2);

    if ( (*pBFESpace1) == (*pBFESpace2) )
    {
        typename PatchInterface<1>::Pointer pInterface12;
        typename PatchInterface<1>::Pointer pInterface21;

        pInterface12 = iga::make_shared<BSplinesPatchInterface<1> >(pPatch1, side1, pPatch2, side2);
        pInterface21 = iga::make_shared<BSplinesPatchInterface<1> >(pPatch2, side2, pPatch1, side1);

        pInterface12->SetOtherInterface(pInterface21);
        pInterface21->SetOtherInterface(pInterface12);

        pPatch1->AddInterface(pInterface12);
        pPatch2->AddInterface(pInterface21);
    }
    else
    {
        KRATOS_ERROR << "The interface is not created because the two patch's boundaries are not conformed.";
    }
}

void BSplinesPatchUtility::MakeInterface1D(typename Patch<2>::Pointer pPatch1, const BoundarySide side1,
        typename Patch<2>::Pointer pPatch2, const BoundarySide side2)
{
    KRATOS_ERROR << "Irrelevant for 2D";
}

void BSplinesPatchUtility::MakeInterface1D(typename Patch<3>::Pointer pPatch1, const BoundarySide side1,
        typename Patch<3>::Pointer pPatch2, const BoundarySide side2)
{
    KRATOS_ERROR << "Irrelevant for 3D";
}

void BSplinesPatchUtility::MakeInterface2D(typename Patch<1>::Pointer pPatch1, const BoundarySide side1,
        typename Patch<1>::Pointer pPatch2, const BoundarySide side2, const BoundaryDirection direction)
{
    KRATOS_ERROR << "Irrelevant for 1D";
}

void BSplinesPatchUtility::MakeInterface2D(typename Patch<2>::Pointer pPatch1, const BoundarySide side1,
        typename Patch<2>::Pointer pPatch2, const BoundarySide side2, const BoundaryDirection direction)
{
    typename FESpace<1>::Pointer pBFESpace1 = pPatch1->pFESpace()->ConstructBoundaryFESpace(side1);

    typename FESpace<1>::Pointer pBFESpace2;

    std::map<std::size_t, std::size_t> local_parameter_map;
    std::vector<BoundaryDirection> directions = {direction};
    pBFESpace2 = pPatch2->pFESpace()->ConstructBoundaryFESpace(side2, local_parameter_map, directions);

    if ( (*pBFESpace1) == (*pBFESpace2) )
    {
        typename PatchInterface<2>::Pointer pInterface12;
        typename PatchInterface<2>::Pointer pInterface21;

        pInterface12 = iga::make_shared<BSplinesPatchInterface<2> >(pPatch1, side1, pPatch2, side2, direction);
        pInterface21 = iga::make_shared<BSplinesPatchInterface<2> >(pPatch2, side2, pPatch1, side1, direction);

        pInterface12->SetOtherInterface(pInterface21);
        pInterface21->SetOtherInterface(pInterface12);

        pPatch1->AddInterface(pInterface12);
        pPatch2->AddInterface(pInterface21);
    }
    else
    {
        KRATOS_ERROR << "The interface is not created because the two patch's boundaries are not conformed.";
    }
}

void BSplinesPatchUtility::MakeInterface2D(typename Patch<3>::Pointer pPatch1, const BoundarySide side1,
        typename Patch<3>::Pointer pPatch2, const BoundarySide side2, const BoundaryDirection direction)
{
    KRATOS_ERROR << "Irrelevant for 3D";
}

void BSplinesPatchUtility::MakeInterface3D(typename Patch<1>::Pointer pPatch1, const BoundarySide side1,
        typename Patch<1>::Pointer pPatch2, const BoundarySide side2, const bool uv_or_vu,
        const BoundaryDirection direction1, const BoundaryDirection direction2)
{
    KRATOS_ERROR << "Irrelevant for 1D";
}

void BSplinesPatchUtility::MakeInterface3D(typename Patch<2>::Pointer pPatch1, const BoundarySide side1,
        typename Patch<2>::Pointer pPatch2, const BoundarySide side2, const bool uv_or_vu,
        const BoundaryDirection direction1, const BoundaryDirection direction2)
{
    KRATOS_ERROR << "Irrelevant for 2D";
}

void BSplinesPatchUtility::MakeInterface3D(typename Patch<3>::Pointer pPatch1, const BoundarySide side1,
        typename Patch<3>::Pointer pPatch2, const BoundarySide side2, const bool uv_or_vu,
        const BoundaryDirection direction1, const BoundaryDirection direction2)
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

    if ( (*pBFESpace1) == (*pBFESpace2) )
    {
        typename PatchInterface<3>::Pointer pInterface12;
        typename PatchInterface<3>::Pointer pInterface21;

        pInterface12 = iga::make_shared<BSplinesPatchInterface<3> >(pPatch1, side1, pPatch2, side2, uv_or_vu, direction1, direction2);
        pInterface21 = iga::make_shared<BSplinesPatchInterface<3> >(pPatch2, side2, pPatch1, side1, uv_or_vu, direction1, direction2);

        pInterface12->SetOtherInterface(pInterface21);
        pInterface21->SetOtherInterface(pInterface12);

        pPatch1->AddInterface(pInterface12);
        pPatch2->AddInterface(pInterface21);
    }
    else
    {
        KRATOS_ERROR << "The interface is not created because the two patch's boundaries are not conformed.";
    }
}

/// template instantiation
template typename Patch<2>::Pointer BSplinesPatchUtility::CreateLoftPatch<2>(typename Patch<1>::Pointer pPatch1, typename Patch<1>::Pointer pPatch2);
template typename Patch<3>::Pointer BSplinesPatchUtility::CreateLoftPatch<3>(typename Patch<2>::Pointer pPatch1, typename Patch<2>::Pointer pPatch2);
template typename Patch<2>::Pointer BSplinesPatchUtility::CreateLoftPatch<2>(std::vector<typename Patch<1>::Pointer> pPatches, int order);
template typename Patch<3>::Pointer BSplinesPatchUtility::CreateLoftPatch<3>(std::vector<typename Patch<2>::Pointer> pPatches, int order);
template typename Patch<2>::Pointer BSplinesPatchUtility::CreatePatchFromGeo<2>(const std::string& fn);
template typename Patch<3>::Pointer BSplinesPatchUtility::CreatePatchFromGeo<3>(const std::string& fn);
template void BSplinesPatchUtility::Reverse<1>(typename Patch<1>::Pointer pPatch, std::size_t idir);
template void BSplinesPatchUtility::Reverse<2>(typename Patch<2>::Pointer pPatch, std::size_t idir);
template void BSplinesPatchUtility::Reverse<3>(typename Patch<3>::Pointer pPatch, std::size_t idir);

} // namespace Kratos
