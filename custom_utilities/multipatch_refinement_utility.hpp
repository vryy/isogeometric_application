//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 8 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_REFINEMENT_UTILITY_HPP_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_REFINEMENT_UTILITY_HPP_INCLUDED

// Project includes
#include "custom_utilities/multipatch_refinement_utility.h"

namespace Kratos
{

/// Insert the knots to the NURBS patch and make it compatible across neighbors
template<int TDim>
void MultiPatchRefinementUtility::InsertKnots(typename Patch<TDim>::Pointer& pPatch, std::set<std::size_t>& refined_patches, const std::vector<std::vector<double> >& ins_knots)
{
    if (pPatch->pFESpace()->Type() != BSplinesFESpace<TDim>::StaticType())
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "only support the NURBS patch")

    if (std::find(refined_patches.begin(), refined_patches.end(), pPatch->Id()) == refined_patches.end())
    {
        // create new patch with same Id
        typename Patch<TDim>::Pointer pNewPatch = typename Patch<TDim>::Pointer(new Patch<TDim>(pPatch->Id()));

        // compute the transformation matrix
        Matrix T;
        std::vector<std::vector<double> > new_knots(TDim);

        typename BSplinesFESpace<TDim>::Pointer pFESpace = boost::dynamic_pointer_cast<BSplinesFESpace<TDim> >(pPatch->pFESpace());
        typename BSplinesFESpace<TDim>::Pointer pNewFESpace = typename BSplinesFESpace<TDim>::Pointer(new BSplinesFESpace<TDim>());

        std::vector<std::size_t> new_size(TDim);

        this->ComputeBsplinesKnotInsertionCoefficients<TDim>(T, new_knots, pFESpace, ins_knots);

        for (std::size_t dim = 0; dim < TDim; ++dim)
        {
            new_size[dim] = new_knots[dim].size() - pPatch->Order(dim) - 1;
            pNewFESpace->SetKnotVector(dim, new_knots[dim]);
            pNewFESpace->SetInfo(dim, new_size[dim], pPatch->Order(dim));
        }

        pNewFESpace->ResetFunctionIndices();

        // KRATOS_WATCH(T)

        // set the new FESpace
        pNewPatch->SetFESpace(pNewFESpace);

        // transform and transfer the control points
        typename ControlGrid<ControlPoint<double> >::Pointer pNewControlPoints = typename ControlGrid<ControlPoint<double> >::Pointer (new StructuredControlGrid<TDim, ControlPoint<double> >(new_size));
        ControlGridUtility::Transform<ControlPoint<double>, Matrix>(T, *(pPatch->pControlPointGridFunction()->pControlGrid()), *pNewControlPoints);
        pNewControlPoints->SetName(pPatch->pControlPointGridFunction()->pControlGrid()->Name());
        pNewPatch->CreateControlPointGridFunction(pNewControlPoints);
//        KRATOS_WATCH(*pNewControlPoints)

        std::vector<double> old_weights = pPatch->GetControlWeights();
        std::vector<double> new_weights = pNewPatch->GetControlWeights();

        typename Patch<TDim>::DoubleGridFunctionContainerType DoubleGridFunctions_ = pPatch->DoubleGridFunctions();

        typename Patch<TDim>::Array1DGridFunctionContainerType Array1DGridFunctions_ = pPatch->Array1DGridFunctions();

        typename Patch<TDim>::VectorGridFunctionContainerType VectorGridFunctions_ = pPatch->VectorGridFunctions();

        // transfer the grid function
        // here to transfer correctly we apply a two-step process:
        // + firstly the old control values is multiplied with weight to make it weighted control values
        // + secondly the control values will be transferred
        // + the new control values will be divided by the new weight to make it unweighted
        for (typename Patch<TDim>::DoubleGridFunctionContainerType::const_iterator it = DoubleGridFunctions_.begin();
                it != DoubleGridFunctions_.end(); ++it)
        {
            typename ControlGrid<double>::Pointer pNewDoubleControlGrid = typename ControlGrid<double>::Pointer (new StructuredControlGrid<TDim, double>(new_size));
            ControlGridUtility::Transform<double, Matrix>(T, old_weights, *((*it)->pControlGrid()), new_weights, *pNewDoubleControlGrid);
            pNewDoubleControlGrid->SetName((*it)->pControlGrid()->Name());
            pNewPatch->template CreateGridFunction<double>(pNewDoubleControlGrid);
        }

        for (typename Patch<TDim>::Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions_.begin();
                it != Array1DGridFunctions_.end(); ++it)
        {
            typename ControlGrid<array_1d<double, 3> >::Pointer pNewArray1DControlGrid = typename ControlGrid<array_1d<double, 3> >::Pointer (new StructuredControlGrid<TDim, array_1d<double, 3> >(new_size));
            ControlGridUtility::Transform<array_1d<double, 3>, Matrix>(T, old_weights, *((*it)->pControlGrid()), new_weights, *pNewArray1DControlGrid);
            pNewArray1DControlGrid->SetName((*it)->pControlGrid()->Name());
            pNewPatch->template CreateGridFunction<array_1d<double, 3> >(pNewArray1DControlGrid);
        }

        for (typename Patch<TDim>::VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions_.begin();
                it != VectorGridFunctions_.end(); ++it)
        {
            typename ControlGrid<Vector>::Pointer pNewVectorControlGrid = typename ControlGrid<Vector>::Pointer (new StructuredControlGrid<TDim, Vector>(new_size));
            ControlGridUtility::Transform<Vector, Matrix>(T, old_weights, *((*it)->pControlGrid()), new_weights, *pNewVectorControlGrid);
            pNewVectorControlGrid->SetName((*it)->pControlGrid()->Name());
            pNewPatch->template CreateGridFunction<Vector>(pNewVectorControlGrid);
        }

        // mark refined patch
        refined_patches.insert(pPatch->Id());

        // transfer the inserted knots to neighbors
        std::vector<std::vector<double> > neib_ins_knots(TDim);

        for (std::size_t i = 0; i < pPatch->NumberOfInterfaces(); ++i)
        {
            typename PatchInterface<TDim>::Pointer pInterface = pPatch->pInterface(i);
            typename Patch<TDim>::Pointer pNeighbor = pInterface->pPatch2();

            if (pNeighbor->pFESpace()->Type() != BSplinesFESpace<TDim>::StaticType())
                KRATOS_THROW_ERROR(std::logic_error, "The FESpace of the neighbor is not BSplinesFESpace", "")

            if (TDim == 2)
            {
                std::size_t dir1 = ParameterDirection<2>::Get(pInterface->Side1());
                std::size_t dir2 = ParameterDirection<2>::Get(pInterface->Side2());

                if (pInterface->Rotation() == _ROTATE_0_)
                    neib_ins_knots[dir2] = ins_knots[dir1];
                else if (pInterface->Rotation() == _ROTATE_180_)
                    neib_ins_knots[dir2] = KnotArray1D<double>::ReverseKnots(ins_knots[dir1]);
            }
            else if (TDim == 3)
            {
                // TODO
                KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "in 3D is not yet implemented")
            }

            InsertKnots<TDim>(pNeighbor, refined_patches, neib_ins_knots);

            pInterface->SetPatch1(pNewPatch);
            pInterface->SetPatch2(pNeighbor);
            pInterface->pOtherInterface()->SetPatch2(pNewPatch);

            pNewPatch->AddInterface(pInterface);
        }

        // get the parent multipatch
        typename MultiPatch<TDim>::Pointer pMultiPatch = pPatch->pParentMultiPatch();

        if (pMultiPatch != NULL)
        {
            // set the parent multipatch
            pNewPatch->pSetParentMultiPatch(pMultiPatch);

            // remove this patch from multipatch
            pMultiPatch->Patches().erase(pPatch->Id());
        }

        // swap
        std::cout << __FUNCTION__ << ": " << pPatch << " is swapped with ";
        pPatch.swap(pNewPatch);
        std::cout << pPatch << std::endl;

        if (pMultiPatch != NULL)
        {
            // replace the corresponding patch in multipatch
            pMultiPatch->Patches().push_back(pPatch);
            pMultiPatch->Patches().Unique();
        }

        std::cout << __FUNCTION__ << " completed for patch " << pPatch->Id() << std::endl;
    }
}


/// Degree elevation for the NURBS patch and make it compatible across neighbors
template<int TDim>
void MultiPatchRefinementUtility::DegreeElevate(typename Patch<TDim>::Pointer& pPatch, std::set<std::size_t>& refined_patches, const std::vector<std::size_t>& order_increment)
{
    if (pPatch->pFESpace()->Type() != BSplinesFESpace<TDim>::StaticType())
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "only support the NURBS patch")

    if (std::find(refined_patches.begin(), refined_patches.end(), pPatch->Id()) == refined_patches.end())
    {
        // create new patch with same Id
        typename Patch<TDim>::Pointer pNewPatch = typename Patch<TDim>::Pointer(new Patch<TDim>(pPatch->Id()));

        // elevate the degree and initialize new patch
        typename BSplinesFESpace<TDim>::Pointer pFESpace = boost::dynamic_pointer_cast<BSplinesFESpace<TDim> >(pPatch->pFESpace());
        typename BSplinesFESpace<TDim>::Pointer pNewFESpace = typename BSplinesFESpace<TDim>::Pointer(new BSplinesFESpace<TDim>());

        std::vector<std::vector<double> > new_knots(TDim);

        std::vector<std::size_t> new_size(TDim);
        for (std::size_t i = 0; i < TDim; ++i)
            new_size[i] = pFESpace->Number(i);

        typename StructuredControlGrid<TDim, ControlPoint<double> >::Pointer pControlPoints
            = boost::dynamic_pointer_cast<StructuredControlGrid<TDim, ControlPoint<double> > >(pPatch->pControlPointGridFunction()->pControlGrid());

        typename StructuredControlGrid<TDim, ControlPoint<double> >::Pointer pNewControlPoints
            = typename StructuredControlGrid<TDim, ControlPoint<double> >::Pointer(new StructuredControlGrid<TDim, ControlPoint<double> >(new_size)); // note here that the size is just temporary, it will be raised later on.

        this->ComputeBsplinesDegreeElevation<TDim, ControlPoint<double> >(*pControlPoints, *pFESpace, order_increment, *pNewControlPoints, new_knots);

        for (std::size_t dim = 0; dim < TDim; ++dim)
        {
            new_size[dim] = new_knots[dim].size() - pFESpace->Order(dim) - order_increment[dim] - 1;
            pNewFESpace->SetKnotVector(dim, new_knots[dim]);
            pNewFESpace->SetInfo(dim, new_size[dim], pFESpace->Order(dim) + order_increment[dim]);
        }

        pNewFESpace->ResetFunctionIndices();

        // set the new FESpace
        pNewPatch->SetFESpace(pNewFESpace);

        pNewControlPoints->SetName(pPatch->pControlPointGridFunction()->pControlGrid()->Name());
        pNewPatch->CreateControlPointGridFunction(pNewControlPoints);

        // raise the order for other control grids
        typename Patch<TDim>::DoubleGridFunctionContainerType DoubleGridFunctions_ = pPatch->DoubleGridFunctions();

        typename Patch<TDim>::Array1DGridFunctionContainerType Array1DGridFunctions_ = pPatch->Array1DGridFunctions();

        typename Patch<TDim>::VectorGridFunctionContainerType VectorGridFunctions_ = pPatch->VectorGridFunctions();

        // transfer the grid function
        // here to transfer correctly we apply a two-step process:
        // + firstly the old control values is multiplied with weight to make it weighted control values
        // + secondly the control values will be transferred
        // + the new control values will be divided by the new weight to make it unweighted
        for (typename Patch<TDim>::DoubleGridFunctionContainerType::const_iterator it = DoubleGridFunctions_.begin();
                it != DoubleGridFunctions_.end(); ++it)
        {
            typename StructuredControlGrid<TDim, double>::Pointer pNewDoubleControlGrid = typename StructuredControlGrid<TDim, double>::Pointer(new StructuredControlGrid<TDim, double>(new_size));
            typename StructuredControlGrid<TDim, double>::Pointer pDoubleControlGrid = boost::dynamic_pointer_cast<StructuredControlGrid<TDim, double> >((*it)->pControlGrid());
            this->ComputeBsplinesDegreeElevation<TDim, double>(*pDoubleControlGrid, *pFESpace, order_increment, *pNewDoubleControlGrid, new_knots);
            pNewDoubleControlGrid->SetName((*it)->pControlGrid()->Name());
            pNewPatch->template CreateGridFunction<double>(pNewDoubleControlGrid);
        }

        for (typename Patch<TDim>::Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions_.begin();
                it != Array1DGridFunctions_.end(); ++it)
        {
            typename StructuredControlGrid<TDim, array_1d<double, 3> >::Pointer pNewArray1DControlGrid = typename StructuredControlGrid<TDim, array_1d<double, 3> >::Pointer(new StructuredControlGrid<TDim, array_1d<double, 3> >(new_size));
            typename StructuredControlGrid<TDim, array_1d<double, 3> >::Pointer pArray1DControlGrid = boost::dynamic_pointer_cast<StructuredControlGrid<TDim, array_1d<double, 3> > >((*it)->pControlGrid());
            this->ComputeBsplinesDegreeElevation<TDim, array_1d<double, 3> >(*pArray1DControlGrid, *pFESpace, order_increment, *pNewArray1DControlGrid, new_knots);
            pNewArray1DControlGrid->SetName((*it)->pControlGrid()->Name());
            pNewPatch->template CreateGridFunction<array_1d<double, 3> >(pNewArray1DControlGrid);
        }

        for (typename Patch<TDim>::VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions_.begin();
                it != VectorGridFunctions_.end(); ++it)
        {
            typename StructuredControlGrid<TDim, Vector>::Pointer pNewVectorControlGrid = typename StructuredControlGrid<TDim, Vector>::Pointer(new StructuredControlGrid<TDim, Vector>(new_size));
            typename StructuredControlGrid<TDim, Vector>::Pointer pVectorControlGrid = boost::dynamic_pointer_cast<StructuredControlGrid<TDim, Vector> >((*it)->pControlGrid());
            this->ComputeBsplinesDegreeElevation<TDim, Vector>(*pVectorControlGrid, *pFESpace, order_increment, *pNewVectorControlGrid, new_knots);
            pNewVectorControlGrid->SetName((*it)->pControlGrid()->Name());
            pNewPatch->template CreateGridFunction<Vector>(pNewVectorControlGrid);
        }

        // mark refined patch
        refined_patches.insert(pPatch->Id());

        // transfer the order increment to neighbors
        std::vector<std::size_t> neib_order_increment(TDim);

        for (std::size_t i = 0; i < pPatch->NumberOfInterfaces(); ++i)
        {
            typename PatchInterface<TDim>::Pointer pInterface = pPatch->pInterface(i);
            typename Patch<TDim>::Pointer pNeighbor = pInterface->pPatch2();

            if (pNeighbor->pFESpace()->Type() != BSplinesFESpace<TDim>::StaticType())
                KRATOS_THROW_ERROR(std::logic_error, "The FESpace of the neighbor is not BSplinesFESpace", "")

            if (TDim == 2)
            {
                std::size_t dir1 = ParameterDirection<2>::Get(pInterface->Side1());
                std::size_t dir2 = ParameterDirection<2>::Get(pInterface->Side2());

                neib_order_increment[dir2] = order_increment[dir1];
            }
            else if (TDim == 3)
            {
                // TODO
                KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "in 3D is not yet implemented")
            }

            DegreeElevate<TDim>(pNeighbor, refined_patches, neib_order_increment);

            pInterface->SetPatch1(pNewPatch);
            pInterface->SetPatch2(pNeighbor);
            pInterface->pOtherInterface()->SetPatch2(pNewPatch);

            pNewPatch->AddInterface(pInterface);
        }

        // get the parent multipatch
        typename MultiPatch<TDim>::Pointer pMultiPatch = pPatch->pParentMultiPatch();

        if (pMultiPatch != NULL)
        {
            // set the parent multipatch
            pNewPatch->pSetParentMultiPatch(pMultiPatch);

            // remove this patch from multipatch
            pMultiPatch->Patches().erase(pPatch->Id());
        }

        // swap
        std::cout << __FUNCTION__ << ": " << pPatch << " is swapped with ";
        pPatch.swap(pNewPatch);
        std::cout << pPatch << std::endl;

        if (pMultiPatch != NULL)
        {
            // replace the corresponding patch in multipatch
            pMultiPatch->Patches().push_back(pPatch);
            pMultiPatch->Patches().Unique();
        }

        std::cout << __FUNCTION__ << " completed for patch " << pPatch->Id() << std::endl;
    }
}


template<>
struct ComputeBsplinesKnotInsertionCoefficients_Helper<1>
{
    static void Compute(Matrix& T,
        std::vector<std::vector<double> >& new_knots,
        typename BSplinesFESpace<1>::Pointer& pFESpace,
        const std::vector<std::vector<double> >& ins_knots)
    {
        BSplineUtils::ComputeBsplinesKnotInsertionCoefficients1D(T,
                new_knots[0],
                pFESpace->Order(0),
                pFESpace->KnotVector(0),
                ins_knots[0]);
    }
};

template<>
struct ComputeBsplinesKnotInsertionCoefficients_Helper<2>
{
    static void Compute(Matrix& T,
        std::vector<std::vector<double> >& new_knots,
        typename BSplinesFESpace<2>::Pointer& pFESpace,
        const std::vector<std::vector<double> >& ins_knots)
    {
        BSplineUtils::ComputeBsplinesKnotInsertionCoefficients2D(T,
                new_knots[0], new_knots[1],
                pFESpace->Order(0), pFESpace->Order(1),
                pFESpace->KnotVector(0), pFESpace->KnotVector(1),
                ins_knots[0], ins_knots[1]);
    }
};

template<>
struct ComputeBsplinesKnotInsertionCoefficients_Helper<3>
{
    static void Compute(Matrix& T,
        std::vector<std::vector<double> >& new_knots,
        typename BSplinesFESpace<3>::Pointer& pFESpace,
        const std::vector<std::vector<double> >& ins_knots)
    {
        BSplineUtils::ComputeBsplinesKnotInsertionCoefficients3D(T,
                new_knots[0], new_knots[1], new_knots[2],
                pFESpace->Order(0), pFESpace->Order(1), pFESpace->Order(2),
                pFESpace->KnotVector(0), pFESpace->KnotVector(1), pFESpace->KnotVector(2),
                ins_knots[0], ins_knots[1], ins_knots[2]);
    }
};

template<typename TDataType>
struct ComputeBsplinesDegreeElevation_Helper<1, TDataType>
{
    static void Compute(const StructuredControlGrid<1, TDataType>& ControlValues,
        const BSplinesFESpace<1>& rFESpace,
        const std::vector<std::size_t>& order_increment,
        StructuredControlGrid<1, TDataType>& NewControlValues,
        std::vector<std::vector<double> >& new_knots)
    {
        TDataType null_control_value(0.0);

        BSplineUtils::ComputeBsplinesDegreeElevation1D(rFESpace.Order(0),
                ControlValues,
                rFESpace.KnotVector(0),
                order_increment[0],
                NewControlValues,
                new_knots[0],
                null_control_value);
    }
};

template<typename TDataType>
struct ComputeBsplinesDegreeElevation_Helper<2, TDataType>
{
    static void Compute(const StructuredControlGrid<2, TDataType>& ControlValues,
        const BSplinesFESpace<2>& rFESpace,
        const std::vector<std::size_t>& order_increment,
        StructuredControlGrid<2, TDataType>& NewControlValues,
        std::vector<std::vector<double> >& new_knots)
    {
        TDataType null_control_value(0.0);

        BSplineUtils::ComputeBsplinesDegreeElevation2D(rFESpace.Order(0), rFESpace.Order(1),
                ControlValues,
                rFESpace.KnotVector(0), rFESpace.KnotVector(1),
                order_increment[0], order_increment[1],
                NewControlValues,
                new_knots[0], new_knots[1],
                null_control_value);
    }
};

template<typename TDataType>
struct ComputeBsplinesDegreeElevation_Helper<3, TDataType>
{
    static void Compute(const StructuredControlGrid<3, TDataType>& ControlValues,
        const BSplinesFESpace<3>& rFESpace,
        const std::vector<std::size_t>& order_increment,
        StructuredControlGrid<3, TDataType>& NewControlValues,
        std::vector<std::vector<double> >& new_knots)
    {
        TDataType null_control_value(0.0);

        BSplineUtils::ComputeBsplinesDegreeElevation3D(rFESpace.Order(0), rFESpace.Order(1), rFESpace.Order(2),
                ControlValues,
                rFESpace.KnotVector(0), rFESpace.KnotVector(1), rFESpace.KnotVector(2),
                order_increment[0], order_increment[1], order_increment[2],
                NewControlValues,
                new_knots[0], new_knots[1], new_knots[2],
                null_control_value);
    }
};

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_REFINEMENT_UTILITY_HPP_INCLUDED defined

