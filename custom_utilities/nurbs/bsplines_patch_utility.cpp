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

template<class TPatchType>
typename TPatchType::Pointer BSplinesPatchUtility::CreateLoftPatch(typename TPatchType::BoundaryPatchType::Pointer pPatch1, typename TPatchType::BoundaryPatchType::Pointer pPatch2)
{
    std::vector<typename TPatchType::BoundaryPatchType::Pointer> pPatches = {pPatch1, pPatch2};
    return CreateLoftPatch<TPatchType>(pPatches, 1);
}

template<class TPatchType>
typename TPatchType::Pointer BSplinesPatchUtility::CreateLoftPatch(std::vector<typename TPatchType::BoundaryPatchType::Pointer> pPatches, int order)
{
    constexpr int Dim = TPatchType::Dim;

    if (pPatches.size() == 0)
    {
        return typename TPatchType::Pointer(new TPatchType(-1));
    }

    // check prerequisites
    for (std::size_t i = 0; i < pPatches.size(); ++i)
    {
        if (pPatches[i]->pFESpace()->Type() != BSplinesFESpace<Dim-1>::StaticType())
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
    typename BSplinesFESpace<Dim-1>::Pointer pFESpace0 = iga::dynamic_pointer_cast<BSplinesFESpace<Dim-1> > (pPatches[0]->pFESpace());
    if (pFESpace0 == nullptr)
        KRATOS_ERROR << "The cast to BSplinesFESpace is failed.";
    typename BSplinesFESpace<Dim>::Pointer pNewFESpace = BSplinesFESpace<Dim>::Create();
    for (std::size_t dim = 0; dim < Dim - 1; ++dim)
    {
        pNewFESpace->SetKnotVector(dim, pFESpace0->KnotVector(dim));
        pNewFESpace->SetInfo(dim, pFESpace0->Number(dim), pFESpace0->Order(dim));
    }

    typename BSplinesFESpace<Dim>::knot_container_t new_knot_vector = BSplinesFESpaceLibrary::CreateUniformOpenKnotVector(pPatches.size(), order);
    pNewFESpace->SetKnotVector(Dim - 1, new_knot_vector);
    pNewFESpace->SetInfo(Dim - 1, pPatches.size(), order);

    // create the new patch
    typename TPatchType::Pointer pNewPatch = typename TPatchType::Pointer(new TPatchType(-1, pNewFESpace));

    // create the new control point grid
    typedef typename TPatchType::ControlPointType ControlPointType;
    typename StructuredControlGrid < Dim - 1, ControlPointType >::Pointer pControlPointGrid0
        = iga::dynamic_pointer_cast < StructuredControlGrid < Dim - 1, ControlPointType > > (pPatches[0]->pControlPointGridFunction()->pControlGrid());
    if (pControlPointGrid0 == nullptr)
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
    std::vector<std::size_t> new_sizes(Dim);
    for (std::size_t dim = 0; dim < Dim - 1; ++dim)
    {
        new_sizes[dim] = pControlPointGrid0->Size(dim);
    }
    new_sizes[Dim - 1] = pPatches.size();
    typename StructuredControlGrid<Dim, ControlPointType>::Pointer pNewControlPointGrid = StructuredControlGrid<Dim, ControlPointType>::Create(new_sizes);
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

template<typename TPatchPointerType>
void BSplinesPatchUtility::Reverse(TPatchPointerType pPatch, std::size_t idir)
{
    std::set<std::size_t> reversed_patches;
    ReverseImpl(pPatch, idir, reversed_patches);
}

template<typename TPatchPointerType>
void BSplinesPatchUtility::ReverseImpl(TPatchPointerType pPatch, std::size_t idir, std::set<std::size_t>& reversed_patches)
{
    typedef typename TPatchPointerType::element_type PatchType;
    typedef typename PatchType::DataType DataType;
    typedef typename PatchType::VectorType VectorType;

    constexpr int Dim = PatchType::Dim;

    typedef BSplinesPatchInterface<Dim, typename PatchType::LocalCoordinateType,
            typename PatchType::CoordinateType, typename PatchType::DataType> BSplinesPatchInterfaceType;

    if (reversed_patches.find(pPatch->Id()) != reversed_patches.end())
    {
        return;
    }

    // std::cout << "Patch " << pPatch->Name() << " will be reversed in direction " << idir << std::endl;

    if (pPatch->pFESpace()->Type() != BSplinesFESpace<Dim>::StaticType())
    {
        KRATOS_ERROR << "Patch " << pPatch->Name() << " is not B-Splines patch. Reverse can't be done.";
    }

    // reverse the FESPace
    typename BSplinesFESpace<Dim>::Pointer pFESpace = iga::dynamic_pointer_cast<BSplinesFESpace<Dim> >(pPatch->pFESpace());
    pFESpace->Reverse(idir);

    // reverse the structured control grid
    typedef typename PatchType::ControlPointType ControlPointType;
    typename StructuredControlGrid<Dim, ControlPointType>::Pointer pControlPointGrid =
        iga::dynamic_pointer_cast<StructuredControlGrid<Dim, ControlPointType> >(pPatch->pControlPointGridFunction()->pControlGrid());
    if (pControlPointGrid != nullptr)
    {
        pControlPointGrid->Reverse(idir);
    }
    else
    {
        KRATOS_ERROR << "The control point grid is not structured";
    }

    typedef typename PatchType::DoubleGridFunctionContainerType DoubleGridFunctionContainerType;
    DoubleGridFunctionContainerType DoubleGridFunctions_ = pPatch->DoubleGridFunctions();
    for (typename DoubleGridFunctionContainerType::const_iterator it = DoubleGridFunctions_.begin();
            it != DoubleGridFunctions_.end(); ++it)
    {
        typename StructuredControlGrid<Dim, DataType>::Pointer pControlValueGrid =
            iga::dynamic_pointer_cast<StructuredControlGrid<Dim, DataType> >((*it)->pControlGrid());
        if (pControlValueGrid != nullptr)
        {
            pControlValueGrid->Reverse(idir);
        }
        else
        {
            KRATOS_ERROR << "The control value grid " << (*it)->pControlGrid()->Name() << " is not structured";
        }
    }

    typedef typename PatchType::Array1DGridFunctionContainerType Array1DGridFunctionContainerType;
    Array1DGridFunctionContainerType Array1DGridFunctions_ = pPatch->Array1DGridFunctions();
    for (typename Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions_.begin();
            it != Array1DGridFunctions_.end(); ++it)
    {
        typename StructuredControlGrid<Dim, array_1d<DataType, 3> >::Pointer pControlValueGrid =
            iga::dynamic_pointer_cast<StructuredControlGrid<Dim, array_1d<DataType, 3> > >((*it)->pControlGrid());
        if ((*it)->pControlGrid()->Name() == "CONTROL_POINT_COORDINATES") { continue; }
        if (pControlValueGrid != nullptr)
        {
            pControlValueGrid->Reverse(idir);
        }
        else
        {
            KRATOS_ERROR << "The control value grid " << (*it)->pControlGrid()->Name() << " is not structured";
        }
    }

    typedef typename PatchType::VectorGridFunctionContainerType VectorGridFunctionContainerType;
    VectorGridFunctionContainerType VectorGridFunctions_ = pPatch->VectorGridFunctions();
    for (typename VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions_.begin();
            it != VectorGridFunctions_.end(); ++it)
    {
        typename StructuredControlGrid<Dim, VectorType>::Pointer pControlValueGrid =
            iga::dynamic_pointer_cast<StructuredControlGrid<Dim, VectorType> >((*it)->pControlGrid());
        if (pControlValueGrid != nullptr)
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
    if constexpr (Dim > 1)
    {
        for (std::size_t i = 0; i < pPatch->NumberOfInterfaces(); ++i)
        {
            typename BSplinesPatchInterfaceType::Pointer pInterface = iga::dynamic_pointer_cast<BSplinesPatchInterfaceType>(pPatch->pInterface(i));

            if (pInterface != nullptr)
            {
                if constexpr (Dim == 2)
                {
                    std::size_t idir1 = static_cast<std::size_t>(ParameterDirection<2>::Get_(pInterface->Side1()));
                    std::size_t idir2 = static_cast<std::size_t>(ParameterDirection<2>::Get_(pInterface->Side2()));

                    if (idir1 == idir)
                    {
                        // the reversed direction aligns with the direction on the interface, therefore
                        // we need to reverse the neighbour patch
                        ReverseImpl(pInterface->pPatch2(), idir2, reversed_patches);
                    }
                    else
                    {
                        // the reversed direction does not align with the direction on the interface.
                        // For 2D, it shall be perpendicular to the direction of the interface.
                        // Therefore, we need to flip the side of the interface..
                        pInterface->FlipSide1();
                        pInterface->pOtherInterface()->FlipSide2();
                    }
                }
                else if constexpr (Dim == 3)
                {
                    std::vector<int> dirs1 = ParameterDirection<3>::Get(pInterface->Side1());
                    std::vector<int> dirs2 = ParameterDirection<3>::Get(pInterface->Side2());

                    if (static_cast<std::size_t>(dirs1[0]) == idir)
                    {
                        std::size_t idir2 = dirs2[pInterface->LocalParameterMapping(0)];

                        // reverse the neighbour patch
                        ReverseImpl(pInterface->pPatch2(), idir2, reversed_patches);
                    }
                    else if (static_cast<std::size_t>(dirs1[1]) == idir)
                    {
                        std::size_t idir2 = dirs2[pInterface->LocalParameterMapping(1)];

                        // reverse the neighbour patch
                        ReverseImpl(pInterface->pPatch2(), idir2, reversed_patches);
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

template<class TPatchPointerType>
void BSplinesPatchUtility::Transpose2D(TPatchPointerType pPatch)
{
    typedef typename TPatchPointerType::element_type PatchType;

    if constexpr (PatchType::Dim == 2)
        Transpose<PatchType>(*pPatch, 0, 1);
    else
        KRATOS_ERROR << "Not a 2D patch";
}

template<class TPatchPointerType>
void BSplinesPatchUtility::Transpose3D(TPatchPointerType pPatch, std::size_t idir, std::size_t jdir)
{
    typedef typename TPatchPointerType::element_type PatchType;

    if constexpr (PatchType::Dim == 3)
        Transpose<PatchType>(*pPatch, idir, jdir);
    else
        KRATOS_ERROR << "Not a 3D patch";
}

template<class TPatchType>
void BSplinesPatchUtility::Transpose(TPatchType& rPatch, std::size_t idir, std::size_t jdir)
{
    typedef TPatchType PatchType;
    typedef typename PatchType::LocalCoordinateType LocalCoordinateType;
    typedef typename PatchType::CoordinateType CoordinateType;
    typedef typename PatchType::DataType DataType;
    typedef typename PatchType::VectorType VectorType;

    constexpr int Dim = PatchType::Dim;

    typedef BSplinesFESpace<Dim, LocalCoordinateType> BSplinesFESpaceType;
    typedef BSplinesPatchInterface<Dim, LocalCoordinateType, CoordinateType, DataType> BSplinesPatchInterfaceType;

    if (rPatch.pFESpace()->Type() != BSplinesFESpaceType::StaticType())
    {
        KRATOS_ERROR << "Patch " << rPatch.Name() << " is not B-Splines patch. Transpose can't be done.";
    }

    // transpose the FESPace
    typename BSplinesFESpaceType::Pointer pFESpace = iga::dynamic_pointer_cast<BSplinesFESpaceType>(rPatch.pFESpace());
    pFESpace->Transpose(idir, jdir);

    // transpose the structured control grid
    typedef typename PatchType::ControlPointType ControlPointType;
    typename StructuredControlGrid<Dim, ControlPointType>::Pointer pControlPointGrid =
        iga::dynamic_pointer_cast<StructuredControlGrid<Dim, ControlPointType> >(rPatch.pControlPointGridFunction()->pControlGrid());
    if (pControlPointGrid != nullptr)
    {
        pControlPointGrid->Transpose(idir, jdir);
    }
    else
    {
        KRATOS_ERROR << "The control point grid is not structured";
    }

    typedef typename PatchType::DoubleGridFunctionContainerType DoubleGridFunctionContainerType;
    DoubleGridFunctionContainerType DoubleGridFunctions_ = rPatch.DoubleGridFunctions();
    for (typename DoubleGridFunctionContainerType::const_iterator it = DoubleGridFunctions_.begin();
            it != DoubleGridFunctions_.end(); ++it)
    {
        typename StructuredControlGrid<Dim, DataType>::Pointer pControlValueGrid =
            iga::dynamic_pointer_cast<StructuredControlGrid<Dim, DataType> >((*it)->pControlGrid());
        if (pControlValueGrid != nullptr)
        {
            pControlValueGrid->Transpose(idir, jdir);
        }
        else
        {
            KRATOS_ERROR << "The control value grid " << (*it)->pControlGrid()->Name() << " is not structured";
        }
    }

    typedef typename PatchType::Array1DGridFunctionContainerType Array1DGridFunctionContainerType;
    Array1DGridFunctionContainerType Array1DGridFunctions_ = rPatch.Array1DGridFunctions();
    for (typename Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions_.begin();
            it != Array1DGridFunctions_.end(); ++it)
    {
        typename StructuredControlGrid<Dim, array_1d<DataType, 3> >::Pointer pControlValueGrid =
            iga::dynamic_pointer_cast<StructuredControlGrid<Dim, array_1d<DataType, 3> > >((*it)->pControlGrid());
        if ((*it)->pControlGrid()->Name() == "CONTROL_POINT_COORDINATES") { continue; }
        if (pControlValueGrid != nullptr)
        {
            pControlValueGrid->Transpose(idir, jdir);
        }
        else
        {
            KRATOS_ERROR << "The control value grid " << (*it)->pControlGrid()->Name() << " is not structured";
        }
    }

    typedef typename PatchType::VectorGridFunctionContainerType VectorGridFunctionContainerType;
    VectorGridFunctionContainerType VectorGridFunctions_ = rPatch.VectorGridFunctions();
    for (typename VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions_.begin();
            it != VectorGridFunctions_.end(); ++it)
    {
        typename StructuredControlGrid<Dim, VectorType>::Pointer pControlValueGrid =
            iga::dynamic_pointer_cast<StructuredControlGrid<Dim, VectorType> >((*it)->pControlGrid());
        if (pControlValueGrid != nullptr)
        {
            pControlValueGrid->Transpose(idir, jdir);
        }
        else
        {
            KRATOS_ERROR << "The control value grid " << (*it)->pControlGrid()->Name() << " is not structured";
        }
    }

    // check for any neighbours and adapt accordingly
    if constexpr (Dim > 1)
    {
        for (std::size_t i = 0; i < rPatch.NumberOfInterfaces(); ++i)
        {
            typename BSplinesPatchInterfaceType::Pointer pInterface = iga::dynamic_pointer_cast<BSplinesPatchInterfaceType>(rPatch.pInterface(i));

            if (pInterface != nullptr)
            {
                if constexpr (Dim == 2)
                {
                    std::size_t idir1 = static_cast<std::size_t>(ParameterDirection<2>::Get_(pInterface->Side1()));

                    BoundarySide new_side;
                    if (idir1 == idir)
                    {
                        new_side = ParameterDirection<2>::GetSide(jdir);
                    }
                    else if (idir1 == jdir)
                    {
                        new_side = ParameterDirection<2>::GetSide(idir);
                    }

                    pInterface->SetSide1(new_side);
                    pInterface->pOtherInterface()->SetSide2(new_side);
                }
                else if constexpr (Dim == 3)
                {
                    // TODO
                    KRATOS_ERROR << "Transpose on the interface is not yet implemented for 3D";
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

template<typename TPatchPointerType>
void BSplinesPatchUtility::MakeInterface1D(TPatchPointerType pPatch1, const BoundarySide side1,
        TPatchPointerType pPatch2, const BoundarySide side2)
{
    typedef typename TPatchPointerType::element_type PatchType;
    constexpr int Dim = PatchType::Dim;

    typedef BSplinesPatchInterface<1, typename PatchType::LocalCoordinateType,
            typename PatchType::CoordinateType, typename PatchType::DataType> BSplinesPatchInterfaceType;

    if constexpr (Dim == 1)
    {
        auto pBFESpace1 = pPatch1->pFESpace()->ConstructBoundaryFESpace(side1);

        auto pBFESpace2 = pPatch2->pFESpace()->ConstructBoundaryFESpace(side2);

        if ( (*pBFESpace1) == (*pBFESpace2) )
        {
            typename PatchType::PatchInterfaceType::Pointer pInterface12;
            typename PatchType::PatchInterfaceType::Pointer pInterface21;

            pInterface12 = iga::make_shared<BSplinesPatchInterfaceType>(pPatch1, side1, pPatch2, side2);
            pInterface21 = iga::make_shared<BSplinesPatchInterfaceType>(pPatch2, side2, pPatch1, side1);

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
    else
        KRATOS_ERROR << "Irrelevant for " << Dim << "D";
}

template<typename TPatchPointerType>
void BSplinesPatchUtility::MakeInterface2D(TPatchPointerType pPatch1, const BoundarySide side1,
        TPatchPointerType pPatch2, const BoundarySide side2, const BoundaryDirection direction)
{
    typedef typename TPatchPointerType::element_type PatchType;

    constexpr int Dim = PatchType::Dim;

    if constexpr (Dim == 2)
    {
        typedef BSplinesPatchInterface<2, typename PatchType::LocalCoordinateType,
                typename PatchType::CoordinateType, typename PatchType::DataType> BSplinesPatchInterfaceType;

        auto pBFESpace1 = pPatch1->pFESpace()->ConstructBoundaryFESpace(side1);

        typename PatchType::BoundaryPatchType::FESpaceType::Pointer pBFESpace2;

        std::map<std::size_t, std::size_t> local_parameter_map;
        std::vector<BoundaryDirection> directions = {direction};
        pBFESpace2 = pPatch2->pFESpace()->ConstructBoundaryFESpace(side2, local_parameter_map, directions);

        if ( (*pBFESpace1) == (*pBFESpace2) )
        {
            typename PatchType::PatchInterfaceType::Pointer pInterface12;
            typename PatchType::PatchInterfaceType::Pointer pInterface21;

            /*
             * The logic of making interface here is simple. It's all about relative position to each other.
             * If patch 1 boundary see patch 2 boundary in the reversed direction, so is the patch 2 boundary.
             * Therefore, the direction information of both two interfaces is the same.
             */

            pInterface12 = iga::make_shared<BSplinesPatchInterfaceType>(pPatch1, side1, pPatch2, side2, direction);
            pInterface21 = iga::make_shared<BSplinesPatchInterfaceType>(pPatch2, side2, pPatch1, side1, direction);

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
    else
        KRATOS_ERROR << "Irrelevant for " << Dim << "D";
}

template<typename TPatchPointerType>
void BSplinesPatchUtility::MakeInterface3D(TPatchPointerType pPatch1, const BoundarySide side1,
        TPatchPointerType pPatch2, const BoundarySide side2, const bool uv_or_vu,
        const BoundaryDirection direction1, const BoundaryDirection direction2)
{
    typedef typename TPatchPointerType::element_type PatchType;

    constexpr int Dim = PatchType::Dim;

    if constexpr (Dim == 3)
    {
        typedef BSplinesPatchInterface<3, typename PatchType::LocalCoordinateType,
                typename PatchType::CoordinateType, typename PatchType::DataType> BSplinesPatchInterfaceType;

        auto pBFESpace1 = pPatch1->pFESpace()->ConstructBoundaryFESpace(side1);

        typename PatchType::BoundaryPatchType::FESpaceType::Pointer pBFESpace2;

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
            typename PatchType::PatchInterfaceType::Pointer pInterface12;
            typename PatchType::PatchInterfaceType::Pointer pInterface21;

            pInterface12 = iga::make_shared<BSplinesPatchInterfaceType>(pPatch1, side1, pPatch2, side2, uv_or_vu, direction1, direction2);

            if (uv_or_vu)
            {
                pInterface21 = iga::make_shared<BSplinesPatchInterfaceType>(pPatch2, side2, pPatch1, side1, uv_or_vu, direction1, direction2);
            }
            else
            {
                // if the local parameter space is swapped, then the direction seeing from the other interface must be reversed
                pInterface21 = iga::make_shared<BSplinesPatchInterfaceType>(pPatch2, side2, pPatch1, side1, uv_or_vu, direction2, direction1);
            }

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
    else
        KRATOS_ERROR << "Irrelevant for " << Dim << "D";
}

template<class TMultiPatchType>
void BSplinesPatchUtility::CreateInterfaces(typename TMultiPatchType::Pointer pMultiPatch)
{
    constexpr int Dim = TMultiPatchType::PatchType::Dim;

    if constexpr (Dim == 1)
    {
        // TODO
        KRATOS_ERROR << "To be implemented";
    }
    else if constexpr (Dim == 2)
    {
        // TODO

        KRATOS_ERROR << "To be implemented";

        typedef typename TMultiPatchType::patch_const_iterator patch_const_iterator;

        // const std::vector<BoundarySide> sides = {}

        // Calculation sequence:
        //  +   compute the average of control points on each side of each patch
        //  +   check interface matching by matching the control point average
        //  +   for each match, check further criterion:
        //      -   all control points on the interface patch must match
        //      -   the knot vector is the same -> the matching direction is forward
        //      -   the knot vector is symmetric around the pivot-> the matching direction is reversed

        // std::map<std::size_t, std::map<std::size_t, PointType> > PatchBoundaryCenters;

        for (patch_const_iterator it = pMultiPatch->begin(); it != pMultiPatch->end(); ++it)
        {

        }
    }
    else if constexpr (Dim == 3)
    {
        // TODO

        KRATOS_ERROR << "To be implemented";
    }
}

template<class TMultiPatchType>
void BSplinesPatchUtility::CheckRepeatedKnot(typename TMultiPatchType::Pointer pMultiPatch)
{
    constexpr int Dim = TMultiPatchType::PatchType::Dim;

    typedef typename TMultiPatchType::PatchType::LocalCoordinateType LocalCoordinateType;
    typedef typename TMultiPatchType::PatchType::CoordinateType CoordinateType;

    typedef BSplinesFESpace<Dim, LocalCoordinateType> BSplinesFESpaceType;

    for (auto it = pMultiPatch->begin(); it != pMultiPatch->end(); ++it)
    {
        typename BSplinesFESpaceType::Pointer pFESpace = iga::dynamic_pointer_cast<BSplinesFESpaceType>(it->pFESpace());

        if (pFESpace == nullptr)
            KRATOS_ERROR << "The underlying FESpace is not BSplinesFESpace";

        for (int i = 0; i < Dim; ++i)
        {
            const auto& knot_vector = pFESpace->KnotVector(i);
            const auto numbers = knot_vector.GetKnotRepeatedNumber();
            for (int i = 1; i < numbers.size() - 1; ++i)
            {
                if (numbers[i] > 1)
                {
                    knot_vector.PrintInfo(std::cout); std::cout << std::endl;
                    KRATOS_WATCH_STD_CON(numbers)
                    KRATOS_ERROR << "There are duplicated knots";
                }
            }
        }
    }
}

std::vector<std::array<typename Patch<1>::ControlPointType, 2> > BSplinesPatchUtility::ExtractControlPolygon(typename Patch<1>::ConstPointer pPatch)
{
    if (pPatch->pFESpace()->Type() != BSplinesFESpace<1>::StaticType())
    {
        KRATOS_ERROR << "Patch " << pPatch->Name() << " is not B-Splines patch.";
    }

    // extract the structured control grid
    typedef typename Patch<1>::ControlPointType ControlPointType;
    typename StructuredControlGrid<1, ControlPointType>::ConstPointer pControlPointGrid =
        iga::dynamic_pointer_cast<const StructuredControlGrid<1, ControlPointType> >(pPatch->pControlPointGridFunction()->pControlGrid());

    std::vector<std::array<ControlPointType, 2> > control_polygon;
    if (pControlPointGrid->size() > 0)
    {
        for (std::size_t i = 0; i < pControlPointGrid->size() - 1; ++i)
        {
            std::array<ControlPointType, 2> line;
            line[0] = pControlPointGrid->GetData(i);
            line[1] = pControlPointGrid->GetData(i+1);
            control_polygon.push_back(line);
        }
    }

    return std::move(control_polygon);
}

/// template instantiation
template typename PatchSelector<2>::RealPatch::Pointer BSplinesPatchUtility::CreatePatchFromGeo<2>(const std::string& fn);
template typename PatchSelector<3>::RealPatch::Pointer BSplinesPatchUtility::CreatePatchFromGeo<3>(const std::string& fn);

//

#define BSplinesPatchUtility_Template_Instantiate(PatchName, MultiPatchName) \
template typename PatchSelector<2>::PatchName::Pointer BSplinesPatchUtility::CreateLoftPatch<PatchSelector<2>::PatchName>(typename PatchSelector<1>::PatchName::Pointer pPatch1, typename PatchSelector<1>::PatchName::Pointer pPatch2);    \
template typename PatchSelector<3>::PatchName::Pointer BSplinesPatchUtility::CreateLoftPatch<PatchSelector<3>::PatchName>(typename PatchSelector<2>::PatchName::Pointer pPatch1, typename PatchSelector<2>::PatchName::Pointer pPatch2);    \
template typename PatchSelector<2>::PatchName::Pointer BSplinesPatchUtility::CreateLoftPatch<PatchSelector<2>::PatchName>(std::vector<typename PatchSelector<1>::PatchName::Pointer> pPatches, int order);  \
template typename PatchSelector<3>::PatchName::Pointer BSplinesPatchUtility::CreateLoftPatch<PatchSelector<3>::PatchName>(std::vector<typename PatchSelector<2>::PatchName::Pointer> pPatches, int order);  \
\
template void BSplinesPatchUtility::MakeInterface1D<PatchSelector<1>::PatchName::Pointer>(typename PatchSelector<1>::PatchName::Pointer pPatch1, const BoundarySide side1, typename PatchSelector<1>::PatchName::Pointer pPatch2, const BoundarySide side2);    \
template void BSplinesPatchUtility::MakeInterface2D<PatchSelector<2>::PatchName::Pointer>(typename PatchSelector<2>::PatchName::Pointer pPatch1, const BoundarySide side1, typename PatchSelector<2>::PatchName::Pointer pPatch2, const BoundarySide side2, const BoundaryDirection direction); \
template void BSplinesPatchUtility::MakeInterface3D<PatchSelector<3>::PatchName::Pointer>(typename PatchSelector<3>::PatchName::Pointer pPatch1, const BoundarySide side1, typename PatchSelector<3>::PatchName::Pointer pPatch2, const BoundarySide side2, const bool uv_or_vu, const BoundaryDirection direction1, const BoundaryDirection direction2);   \
\
template void BSplinesPatchUtility::Reverse<PatchSelector<1>::PatchName::Pointer>(typename PatchSelector<1>::PatchName::Pointer pPatch, std::size_t idir);   \
template void BSplinesPatchUtility::Reverse<PatchSelector<2>::PatchName::Pointer>(typename PatchSelector<2>::PatchName::Pointer pPatch, std::size_t idir);   \
template void BSplinesPatchUtility::Reverse<PatchSelector<3>::PatchName::Pointer>(typename PatchSelector<3>::PatchName::Pointer pPatch, std::size_t idir);   \
template void BSplinesPatchUtility::Transpose2D<PatchSelector<2>::PatchName::Pointer>(PatchSelector<2>::PatchName::Pointer pPatch); \
template void BSplinesPatchUtility::Transpose3D<PatchSelector<3>::PatchName::Pointer>(PatchSelector<3>::PatchName::Pointer pPatch, std::size_t idir, std::size_t jdir); \
template void BSplinesPatchUtility::Transpose<PatchSelector<2>::PatchName>(PatchSelector<2>::PatchName& rPatch, std::size_t idir, std::size_t jdir);    \
template void BSplinesPatchUtility::Transpose<PatchSelector<3>::PatchName>(PatchSelector<3>::PatchName& rPatch, std::size_t idir, std::size_t jdir);    \
template void BSplinesPatchUtility::CreateInterfaces<PatchSelector<1>::MultiPatchName>(typename PatchSelector<1>::MultiPatchName::Pointer pMultiPatch); \
template void BSplinesPatchUtility::CreateInterfaces<PatchSelector<2>::MultiPatchName>(typename PatchSelector<2>::MultiPatchName::Pointer pMultiPatch); \
template void BSplinesPatchUtility::CreateInterfaces<PatchSelector<3>::MultiPatchName>(typename PatchSelector<3>::MultiPatchName::Pointer pMultiPatch); \
template void BSplinesPatchUtility::CheckRepeatedKnot<PatchSelector<1>::MultiPatchName>(typename PatchSelector<1>::MultiPatchName::Pointer pMultiPatch);    \
template void BSplinesPatchUtility::CheckRepeatedKnot<PatchSelector<2>::MultiPatchName>(typename PatchSelector<2>::MultiPatchName::Pointer pMultiPatch);    \
template void BSplinesPatchUtility::CheckRepeatedKnot<PatchSelector<3>::MultiPatchName>(typename PatchSelector<3>::MultiPatchName::Pointer pMultiPatch);    \

BSplinesPatchUtility_Template_Instantiate(RealPatch, RealMultiPatch);
BSplinesPatchUtility_Template_Instantiate(ComplexPatch, ComplexMultiPatch);

} // namespace Kratos
