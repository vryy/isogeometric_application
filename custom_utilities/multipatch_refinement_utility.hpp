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

//#ifdef DEBUG_INS_KNOTS

namespace Kratos
{

/// Insert the knots to the NURBS patch and make it compatible across neighbors
template<int TDim>
void MultiPatchRefinementUtility::InsertKnots(typename Patch<TDim>::Pointer& pPatch,
    std::map<std::size_t, std::vector<int> >& refined_patches,
    const std::vector<std::vector<double> >& ins_knots,
    std::map<std::size_t, Matrix>& trans_mats,
    bool record_trans_mat)
{
    if (pPatch->pFESpace()->Type() != BSplinesFESpace<TDim>::StaticType())
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "only support the NURBS patch")

    #ifdef DEBUG_INS_KNOTS
    std::cout << "ins_knots:";
    for (unsigned int i = 0; i < ins_knots.size(); ++i)
    {
        for (unsigned int j = 0; j < ins_knots[i].size(); ++j)
        {
            std::cout << " " << ins_knots[i][j];
        }
        std::cout << ";";
    }
    std::cout << std::endl;
    #endif

    bool to_refine = false;
    if (refined_patches.find( pPatch->Id() ) == refined_patches.end())
    // the patch is not yet refined
    {
        to_refine = true;
    }
    else
    // the patch is refined but it has some unrefined directions
    {
        for (unsigned int i = 0; i < TDim; ++i)
        {
            if ( (ins_knots[i].size() != 0) && (refined_patches[pPatch->Id()][i] == 0) )
                to_refine = true;
        }
    }

    if (to_refine)
    {
        // create new patch with same Id
        typename Patch<TDim>::Pointer pNewPatch = typename Patch<TDim>::Pointer(new Patch<TDim>(pPatch->Id()));

        // compute the transformation matrix
        std::vector<std::vector<double> > new_knots(TDim);

        typename BSplinesFESpace<TDim>::Pointer pFESpace = boost::dynamic_pointer_cast<BSplinesFESpace<TDim> >(pPatch->pFESpace());
        if (pFESpace == NULL)
            KRATOS_THROW_ERROR(std::runtime_error, "The cast to BSplinesFESpace is failed.", "")
        typename BSplinesFESpace<TDim>::Pointer pNewFESpace = typename BSplinesFESpace<TDim>::Pointer(new BSplinesFESpace<TDim>());

        Matrix T;
        this->ComputeBsplinesKnotInsertionCoefficients<TDim>(T, new_knots, pFESpace, ins_knots);

        std::vector<std::size_t> new_size(TDim);
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

        // transfer the grid function
        // here to transfer correctly we apply a two-step process:
        // + firstly the old control values is multiplied with weight to make it weighted control values
        // + secondly the control values will be transferred
        // + the new control values will be divided by the new weight to make it unweighted

        std::vector<double> old_weights = pPatch->GetControlWeights();
        std::vector<double> new_weights = pNewPatch->GetControlWeights();

        if (record_trans_mat)
        {
            Matrix M = trans(T);

            for (std::size_t i = 0; i < new_weights.size(); ++i)
                row(M, i) /= new_weights[i];

            for (std::size_t j = 0; j < old_weights.size(); ++j)
                column(M, j) *= old_weights[j];

            trans_mats[pPatch->Id()] = M;
        }

        typename Patch<TDim>::DoubleGridFunctionContainerType DoubleGridFunctions_ = pPatch->DoubleGridFunctions();

        typename Patch<TDim>::Array1DGridFunctionContainerType Array1DGridFunctions_ = pPatch->Array1DGridFunctions();

        typename Patch<TDim>::VectorGridFunctionContainerType VectorGridFunctions_ = pPatch->VectorGridFunctions();

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
            if ((*it)->pControlGrid()->Name() == "CONTROL_POINT_COORDINATES") continue;
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
        if (refined_patches.find(pPatch->Id()) == refined_patches.end())
        {
            refined_patches[pPatch->Id()].resize(TDim);
            for (unsigned int i = 0; i < TDim; ++i)
                refined_patches[pPatch->Id()][i] = 0;
        }
        for (unsigned int i = 0; i < TDim; ++i)
        {
            if (ins_knots[i].size() != 0)
            {
                refined_patches[pPatch->Id()][i] = 1;
            }
        }

        for (std::size_t ii = 0; ii < pPatch->NumberOfInterfaces(); ++ii)
        {
            typename BSplinesPatchInterface<TDim>::Pointer pInterface = boost::dynamic_pointer_cast<BSplinesPatchInterface<TDim> >(pPatch->pInterface(ii));
            if (pInterface == NULL)
                KRATOS_THROW_ERROR(std::runtime_error, "The cast to BSplinesPatchInterface is failed", "")
            typename Patch<TDim>::Pointer pNeighbor = pInterface->pPatch2();

            if (pNeighbor->pFESpace()->Type() != BSplinesFESpace<TDim>::StaticType())
                KRATOS_THROW_ERROR(std::logic_error, "The FESpace of the neighbor is not BSplinesFESpace", "")

            // transfer the inserted knots to neighbors
            std::vector<std::vector<double> > neib_ins_knots(TDim);
            for (unsigned int i = 0; i < TDim; ++i)
                neib_ins_knots[i].resize(0);

            if (TDim == 2)
            {
                int dir1 = ParameterDirection<2>::Get_(pInterface->Side1());
                int dir2 = ParameterDirection<2>::Get_(pInterface->Side2());

                neib_ins_knots[dir2] = KnotArray1D<double>::CloneKnots(ins_knots[dir1], pInterface->Direction(0));
                #ifdef DEBUG_INS_KNOTS
                std::cout << "Propagate [";
                for (unsigned int i = 0; i < neib_ins_knots[dir2].size(); ++i)
                    std::cout << ", " << neib_ins_knots[dir2][i];
                std::cout << "] from Patch " << pPatch->Id() << " dir " << dir1;
                std::cout << " to " << pNeighbor->Id() << " dir " << dir2 << std::endl;
                #endif
            }
            else if (TDim == 3)
            {
                std::vector<int> param_dirs_1 = ParameterDirection<3>::Get(pInterface->Side1());
                std::vector<int> param_dirs_2 = ParameterDirection<3>::Get(pInterface->Side2());

                neib_ins_knots[param_dirs_2[ pInterface->LocalParameterMapping(0) ] ] = KnotArray1D<double>::CloneKnots(ins_knots[param_dirs_1[0]], pInterface->Direction(0));
                neib_ins_knots[param_dirs_2[ pInterface->LocalParameterMapping(1) ] ] = KnotArray1D<double>::CloneKnots(ins_knots[param_dirs_1[1]], pInterface->Direction(1));
            }

            InsertKnots<TDim>(pNeighbor, refined_patches, neib_ins_knots, trans_mats, record_trans_mat);

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
void MultiPatchRefinementUtility::DegreeElevate(typename Patch<TDim>::Pointer& pPatch, std::map<std::size_t, std::vector<int> >& refined_patches, const std::vector<std::size_t>& order_increment)
{
    if (pPatch->pFESpace()->Type() != BSplinesFESpace<TDim>::StaticType())
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "only support the NURBS patch")

    bool to_refine = false;
    if (refined_patches.find( pPatch->Id() ) == refined_patches.end())
    // the patch is not yet refined
    {
        to_refine = true;
    }
    else
    // the patch is refined but it has some unrefined directions
    {
        for (unsigned int i = 0; i < TDim; ++i)
        {
            if ( (order_increment[i] != 0) && (refined_patches[pPatch->Id()][i] == 0) )
                to_refine = true;
        }
    }

    if (to_refine)
    {
        // create new patch with same Id
        typename Patch<TDim>::Pointer pNewPatch = typename Patch<TDim>::Pointer(new Patch<TDim>(pPatch->Id()));

        // elevate the degree and initialize new patch
        typename BSplinesFESpace<TDim>::Pointer pFESpace = boost::dynamic_pointer_cast<BSplinesFESpace<TDim> >(pPatch->pFESpace());
        if (pFESpace == NULL)
            KRATOS_THROW_ERROR(std::runtime_error, "The cast to BSplinesFESpace is failed.", "")
        typename BSplinesFESpace<TDim>::Pointer pNewFESpace = typename BSplinesFESpace<TDim>::Pointer(new BSplinesFESpace<TDim>());

        std::vector<std::vector<double> > new_knots(TDim);

        std::vector<std::size_t> new_size(TDim);
        for (std::size_t i = 0; i < TDim; ++i)
            new_size[i] = pFESpace->Number(i);

        typename StructuredControlGrid<TDim, ControlPoint<double> >::Pointer pControlPoints
            = boost::dynamic_pointer_cast<StructuredControlGrid<TDim, ControlPoint<double> > >(pPatch->pControlPointGridFunction()->pControlGrid());
        if (pControlPoints == NULL)
            KRATOS_THROW_ERROR(std::runtime_error, "The cast to StructuredControlGrid is failed.", "")

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
            if (pDoubleControlGrid == NULL)
                KRATOS_THROW_ERROR(std::runtime_error, "The cast to StructuredControlGrid is failed.", "")
            this->ComputeBsplinesDegreeElevation<TDim, double>(*pDoubleControlGrid, *pFESpace, order_increment, *pNewDoubleControlGrid, new_knots);
            pNewDoubleControlGrid->SetName((*it)->pControlGrid()->Name());
            pNewPatch->template CreateGridFunction<double>(pNewDoubleControlGrid);
        }

        for (typename Patch<TDim>::Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions_.begin();
                it != Array1DGridFunctions_.end(); ++it)
        {
            if ((*it)->pControlGrid()->Name() == "CONTROL_POINT_COORDINATES") continue;
            typename StructuredControlGrid<TDim, array_1d<double, 3> >::Pointer pNewArray1DControlGrid = typename StructuredControlGrid<TDim, array_1d<double, 3> >::Pointer(new StructuredControlGrid<TDim, array_1d<double, 3> >(new_size));
            typename StructuredControlGrid<TDim, array_1d<double, 3> >::Pointer pArray1DControlGrid = boost::dynamic_pointer_cast<StructuredControlGrid<TDim, array_1d<double, 3> > >((*it)->pControlGrid());
            if (pArray1DControlGrid == NULL)
                KRATOS_THROW_ERROR(std::runtime_error, "The cast to StructuredControlGrid is failed.", "")
            this->ComputeBsplinesDegreeElevation<TDim, array_1d<double, 3> >(*pArray1DControlGrid, *pFESpace, order_increment, *pNewArray1DControlGrid, new_knots);
            pNewArray1DControlGrid->SetName((*it)->pControlGrid()->Name());
            pNewPatch->template CreateGridFunction<array_1d<double, 3> >(pNewArray1DControlGrid);
        }

        for (typename Patch<TDim>::VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions_.begin();
                it != VectorGridFunctions_.end(); ++it)
        {
            typename StructuredControlGrid<TDim, Vector>::Pointer pNewVectorControlGrid = typename StructuredControlGrid<TDim, Vector>::Pointer(new StructuredControlGrid<TDim, Vector>(new_size));
            typename StructuredControlGrid<TDim, Vector>::Pointer pVectorControlGrid = boost::dynamic_pointer_cast<StructuredControlGrid<TDim, Vector> >((*it)->pControlGrid());
            if (pVectorControlGrid == NULL)
                KRATOS_THROW_ERROR(std::runtime_error, "The cast to StructuredControlGrid is failed.", "")
            this->ComputeBsplinesDegreeElevation<TDim, Vector>(*pVectorControlGrid, *pFESpace, order_increment, *pNewVectorControlGrid, new_knots);
            pNewVectorControlGrid->SetName((*it)->pControlGrid()->Name());
            pNewPatch->template CreateGridFunction<Vector>(pNewVectorControlGrid);
        }

        // mark refined patch
        if (refined_patches.find(pPatch->Id()) == refined_patches.end())
        {
            refined_patches[pPatch->Id()].resize(TDim);
            for (unsigned int i = 0; i < TDim; ++i)
                refined_patches[pPatch->Id()][i] = 0;
        }
        for (unsigned int i = 0; i < TDim; ++i)
        {
            if (order_increment[i] != 0)
            {
                refined_patches[pPatch->Id()][i] = 1;
            }
        }

        for (std::size_t ii = 0; ii < pPatch->NumberOfInterfaces(); ++ii)
        {
            typename BSplinesPatchInterface<TDim>::Pointer pInterface = boost::dynamic_pointer_cast<BSplinesPatchInterface<TDim> >(pPatch->pInterface(ii));
            if (pInterface == NULL)
                KRATOS_THROW_ERROR(std::runtime_error, "The cast to BSplinesPatchInterface is failed", "")
            typename Patch<TDim>::Pointer pNeighbor = pInterface->pPatch2();

            if (pNeighbor->pFESpace()->Type() != BSplinesFESpace<TDim>::StaticType())
                KRATOS_THROW_ERROR(std::logic_error, "The FESpace of the neighbor is not BSplinesFESpace", "")

            // transfer the order increment to neighbors
            std::vector<std::size_t> neib_order_increment(TDim);
            for (unsigned int i = 0; i < TDim; ++i)
                neib_order_increment[i] = 0;

            if (TDim == 2)
            {
                int dir1 = ParameterDirection<2>::Get_(pInterface->Side1());
                int dir2 = ParameterDirection<2>::Get_(pInterface->Side2());

                neib_order_increment[dir2] = order_increment[dir1];
            }
            else if (TDim == 3)
            {
                std::vector<int> param_dirs_1 = ParameterDirection<3>::Get(pInterface->Side1());
                std::vector<int> param_dirs_2 = ParameterDirection<3>::Get(pInterface->Side2());

                neib_order_increment[param_dirs_2[ pInterface->LocalParameterMapping(0) ] ] = order_increment[param_dirs_1[0]];
                neib_order_increment[param_dirs_2[ pInterface->LocalParameterMapping(1) ] ] = order_increment[param_dirs_1[1]];
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

#undef DEBUG_INS_KNOTS

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_REFINEMENT_UTILITY_HPP_INCLUDED defined

