//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 8 Nov 2017 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_REFINEMENT_UTILITY_HPP_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_REFINEMENT_UTILITY_HPP_INCLUDED

// #define DEBUG_INS_KNOTS
// #define ENABLE_VERBOSE_REFINEMENT

namespace Kratos
{

/// Insert the knots to the NURBS patch and make it compatible across neighbors
template<class TPatchType>
void MultiPatchRefinementUtility::InsertKnots(typename TPatchType::Pointer& pPatch,
    std::map<std::size_t, std::vector<int> >& refined_patches,
    const std::vector<std::vector<typename TPatchType::LocalCoordinateType> >& ins_knots,
    std::map<std::size_t, Matrix>& trans_mats,
    bool record_trans_mat)
{
    constexpr int Dim = TPatchType::Dim;

    typedef typename TPatchType::DataType DataType;
    typedef typename MatrixVectorTypeSelector<DataType>::VectorType VectorType;
    typedef typename TPatchType::LocalCoordinateType LocalCoordinateType;
    typedef typename TPatchType::CoordinateType CoordinateType;
    typedef BSplinesFESpace<Dim, LocalCoordinateType> BSplinesFESpaceType;
    typedef BSplinesPatchInterface<Dim, LocalCoordinateType, CoordinateType, DataType> BSplinesPatchInterfaceType;

    if (pPatch->pFESpace()->Type() != BSplinesFESpaceType::StaticType())
    {
        KRATOS_ERROR << "Only support the NURBS patch";
    }

    bool to_refine = false;
    if (refined_patches.find( pPatch->Id() ) == refined_patches.end())
    // the patch is not yet refined
    {
        for (unsigned int i = 0; i < Dim; ++i)
        {
            if ( ins_knots[i].size() != 0 )
                to_refine = true;
        }
    }
    else
    // the patch is refined but it has some unrefined directions
    {
        for (unsigned int i = 0; i < Dim; ++i)
        {
            if ( (ins_knots[i].size() != 0) && (refined_patches[pPatch->Id()][i] == 0) )
                to_refine = true;
        }
    }

    if (to_refine)
    {
        #ifdef ENABLE_VERBOSE_REFINEMENT
        std::cout << "Patch " << pPatch->Id() << " will be refined" << std::endl;
        #endif

        #ifdef DEBUG_INS_KNOTS
        std::cout << " ins_knots:";
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

        // create new patch with same Id
        typename TPatchType::Pointer pNewPatch = typename TPatchType::Pointer(new TPatchType(pPatch->Id()));
        pNewPatch->SetPrefix(pPatch->Prefix());
        pNewPatch->SetLayerIndex(pPatch->LayerIndex());

        // compute the transformation matrix
        std::vector<std::vector<double> > new_knots(Dim);

        typename BSplinesFESpaceType::Pointer pFESpace = iga::dynamic_pointer_cast<BSplinesFESpaceType>(pPatch->pFESpace());
        if (pFESpace == nullptr)
            KRATOS_ERROR << "The cast to BSplinesFESpace is failed.";
        typename BSplinesFESpaceType::Pointer pNewFESpace = typename BSplinesFESpaceType::Pointer(new BSplinesFESpaceType());

        Matrix T;
        this->ComputeBsplinesKnotInsertionCoefficients<Dim>(T, new_knots, pFESpace, ins_knots);

        std::vector<std::size_t> new_size(Dim);
        for (std::size_t dim = 0; dim < Dim; ++dim)
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
        typename ControlGrid<ControlPoint<double> >::Pointer pNewControlPoints = typename ControlGrid<ControlPoint<double> >::Pointer (new StructuredControlGrid<Dim, ControlPoint<double> >(new_size));
        ControlGridUtility::Transform<ControlPoint<double>, Matrix>(T, *(pPatch->pControlPointGridFunction()->pControlGrid()), *pNewControlPoints);
        pNewControlPoints->SetName(pPatch->pControlPointGridFunction()->pControlGrid()->Name());
        pNewPatch->CreateControlPointGridFunction(pNewControlPoints);

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

        typename TPatchType::DoubleGridFunctionContainerType DoubleGridFunctions_ = pPatch->DoubleGridFunctions();

        typename TPatchType::Array1DGridFunctionContainerType Array1DGridFunctions_ = pPatch->Array1DGridFunctions();

        typename TPatchType::VectorGridFunctionContainerType VectorGridFunctions_ = pPatch->VectorGridFunctions();

        for (typename TPatchType::DoubleGridFunctionContainerType::const_iterator it = DoubleGridFunctions_.begin();
                it != DoubleGridFunctions_.end(); ++it)
        {
            typename ControlGrid<DataType>::Pointer pNewDoubleControlGrid = typename ControlGrid<DataType>::Pointer (new StructuredControlGrid<Dim, DataType>(new_size));
            ControlGridUtility::Transform<DataType, Matrix>(T, old_weights, *((*it)->pControlGrid()), new_weights, *pNewDoubleControlGrid);
            pNewDoubleControlGrid->SetName((*it)->pControlGrid()->Name());
            pNewPatch->template CreateGridFunction<DataType>(pNewDoubleControlGrid);
        }

        for (typename TPatchType::Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions_.begin();
                it != Array1DGridFunctions_.end(); ++it)
        {
            if ((*it)->pControlGrid()->Name() == "CONTROL_POINT_COORDINATES") continue;
            typename ControlGrid<array_1d<DataType, 3> >::Pointer pNewArray1DControlGrid = typename ControlGrid<array_1d<DataType, 3> >::Pointer (new StructuredControlGrid<Dim, array_1d<DataType, 3> >(new_size));
            ControlGridUtility::Transform<array_1d<DataType, 3>, Matrix>(T, old_weights, *((*it)->pControlGrid()), new_weights, *pNewArray1DControlGrid);
            pNewArray1DControlGrid->SetName((*it)->pControlGrid()->Name());
            pNewPatch->template CreateGridFunction<array_1d<DataType, 3> >(pNewArray1DControlGrid);
        }

        for (typename TPatchType::VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions_.begin();
                it != VectorGridFunctions_.end(); ++it)
        {
            typename ControlGrid<VectorType>::Pointer pNewVectorControlGrid = typename ControlGrid<VectorType>::Pointer (new StructuredControlGrid<Dim, VectorType>(new_size));
            ControlGridUtility::Transform<VectorType, Matrix>(T, old_weights, *((*it)->pControlGrid()), new_weights, *pNewVectorControlGrid);
            pNewVectorControlGrid->SetName((*it)->pControlGrid()->Name());
            pNewPatch->template CreateGridFunction<VectorType>(pNewVectorControlGrid);
        }

        // mark refined patch
        if (refined_patches.find(pPatch->Id()) == refined_patches.end())
        {
            refined_patches[pPatch->Id()].resize(Dim);
            for (unsigned int i = 0; i < Dim; ++i)
                refined_patches[pPatch->Id()][i] = 0;
        }
        for (unsigned int i = 0; i < Dim; ++i)
        {
            if (ins_knots[i].size() != 0)
            {
                refined_patches[pPatch->Id()][i] = 1;
            }
        }

        for (typename TPatchType::interface_iterator it = pPatch->InterfaceBegin(); it != pPatch->InterfaceEnd(); ++it)
        {
            typename BSplinesPatchInterfaceType::Pointer pInterface = iga::dynamic_pointer_cast<BSplinesPatchInterfaceType>(*it);
            if (pInterface == nullptr)
                KRATOS_ERROR << "The cast to BSplinesPatchInterface is failed";
            typename TPatchType::Pointer pNeighbor = pInterface->pPatch2();

            if (pNeighbor->pFESpace()->Type() != BSplinesFESpaceType::StaticType())
                KRATOS_ERROR << "The FESpace of the neighbor is not BSplinesFESpace";

            // transfer the inserted knots to neighbors
            std::vector<std::vector<double> > neib_ins_knots(Dim);
            for (unsigned int i = 0; i < Dim; ++i)
                neib_ins_knots[i].resize(0);

            if constexpr (Dim == 2)
            {
                int dir1 = ParameterDirection<2>::Get_(pInterface->Side1());
                int dir2 = ParameterDirection<2>::Get_(pInterface->Side2());

                neib_ins_knots[dir2] = KnotArray1D<double>::CloneKnotsWithPivot(new_knots[dir1].back(), ins_knots[dir1], pInterface->Direction(0));
                #ifdef DEBUG_INS_KNOTS
                std::cout << "Propagate [";
                for (unsigned int i = 0; i < neib_ins_knots[dir2].size(); ++i)
                    std::cout << ", " << neib_ins_knots[dir2][i];
                std::cout << "] from Patch " << pPatch->Id() << " dir " << dir1;
                std::cout << " to " << pNeighbor->Id() << " dir " << dir2 << std::endl;
                #endif
            }
            else if constexpr (Dim == 3)
            {
                std::vector<int> param_dirs_1 = ParameterDirection<3>::Get(pInterface->Side1());
                std::vector<int> param_dirs_2 = ParameterDirection<3>::Get(pInterface->Side2());

                neib_ins_knots[param_dirs_2[ pInterface->LocalParameterMapping(0) ] ] = KnotArray1D<double>::CloneKnotsWithPivot(new_knots[param_dirs_1[0]].back(), ins_knots[param_dirs_1[0]], pInterface->Direction(0));
                neib_ins_knots[param_dirs_2[ pInterface->LocalParameterMapping(1) ] ] = KnotArray1D<double>::CloneKnotsWithPivot(new_knots[param_dirs_1[1]].back(), ins_knots[param_dirs_1[1]], pInterface->Direction(1));
            }

            #ifdef ENABLE_VERBOSE_REFINEMENT
            std::cout << "Neighbor patch " << pNeighbor->Id() << " of patch " << pPatch->Id() << " is accounted" << std::endl;
            #endif

            InsertKnots<TPatchType>(pNeighbor, refined_patches, neib_ins_knots, trans_mats, record_trans_mat);

            pInterface->SetPatch1(pNewPatch);
            pInterface->SetPatch2(pNeighbor);
            pInterface->pOtherInterface()->SetPatch2(pNewPatch);

            pNewPatch->AddInterface(pInterface);
        }

        // get the parent multipatch
        auto pMultiPatch = pPatch->pParentMultiPatch();

        if (pMultiPatch != nullptr)
        {
            // set the parent multipatch
            pNewPatch->pSetParentMultiPatch(pMultiPatch);

            // remove this patch from multipatch
            pMultiPatch->Patches().erase(pPatch->Id());
        }

        // swap
        #ifdef ENABLE_VERBOSE_REFINEMENT
        std::cout << __FUNCTION__ << ": " << pPatch << " is swapped with ";
        #endif
        pPatch.swap(pNewPatch);
        #ifdef ENABLE_VERBOSE_REFINEMENT
        std::cout << pPatch << std::endl;
        #endif

        if (pMultiPatch != nullptr)
        {
            // replace the corresponding patch in multipatch
            pMultiPatch->Patches().push_back(pPatch);
            pMultiPatch->Patches().Unique();
        }

        #ifdef ENABLE_VERBOSE_REFINEMENT
        std::cout << __FUNCTION__ << " completed for patch " << pPatch->Id() << std::endl;
        #endif
    }
}

/// Degree elevation for the NURBS patch and make it compatible across neighbors
template<class TPatchType>
void MultiPatchRefinementUtility::DegreeElevate(typename TPatchType::Pointer& pPatch, std::map<std::size_t, std::vector<int> >& refined_patches, const std::vector<std::size_t>& order_increment)
{
    constexpr int Dim = TPatchType::Dim;

    typedef typename TPatchType::DataType DataType;
    typedef typename TPatchType::LocalCoordinateType LocalCoordinateType;
    typedef typename TPatchType::CoordinateType CoordinateType;
    typedef BSplinesFESpace<Dim, LocalCoordinateType> BSplinesFESpaceType;
    typedef BSplinesPatchInterface<Dim, LocalCoordinateType, CoordinateType, DataType> BSplinesPatchInterfaceType;

    if (pPatch->pFESpace()->Type() != BSplinesFESpaceType::StaticType())
        KRATOS_ERROR << "Only support the NURBS patch";

    bool to_refine = false;
    if (refined_patches.find( pPatch->Id() ) == refined_patches.end())
    // the patch is not yet refined
    {
        to_refine = true;
    }
    else
    // the patch is refined but it has some unrefined directions
    {
        for (unsigned int i = 0; i < Dim; ++i)
        {
            if ( (order_increment[i] != 0) && (refined_patches[pPatch->Id()][i] == 0) )
                to_refine = true;
        }
    }

    if (to_refine)
    {
        // create new patch with same Id
        typename TPatchType::Pointer pNewPatch = typename TPatchType::Pointer(new TPatchType(pPatch->Id()));
        pNewPatch->SetPrefix(pPatch->Prefix());
        pNewPatch->SetLayerIndex(pPatch->LayerIndex());

        // elevate the degree and initialize new patch
        typename BSplinesFESpaceType::Pointer pFESpace = iga::dynamic_pointer_cast<BSplinesFESpaceType>(pPatch->pFESpace());
        if (pFESpace == nullptr)
            KRATOS_ERROR << "The cast to BSplinesFESpace is failed.";
        typename BSplinesFESpaceType::Pointer pNewFESpace = typename BSplinesFESpaceType::Pointer(new BSplinesFESpaceType());

        std::vector<std::vector<double> > new_knots(Dim);

        std::vector<std::size_t> new_size(Dim);
        for (std::size_t i = 0; i < Dim; ++i)
            new_size[i] = pFESpace->Number(i);

        typename StructuredControlGrid<Dim, ControlPoint<double> >::Pointer pControlPoints
            = iga::dynamic_pointer_cast<StructuredControlGrid<Dim, ControlPoint<double> > >(pPatch->pControlPointGridFunction()->pControlGrid());
        if (pControlPoints == nullptr)
            KRATOS_ERROR << "The cast to StructuredControlGrid is failed.";

        typename StructuredControlGrid<Dim, ControlPoint<double> >::Pointer pNewControlPoints
            = typename StructuredControlGrid<Dim, ControlPoint<double> >::Pointer(new StructuredControlGrid<Dim, ControlPoint<double> >(new_size)); // note here that the size is just temporary, it will be raised later on.

        this->ComputeBsplinesDegreeElevation<Dim, ControlPoint<double> >(*pControlPoints, *pFESpace, order_increment, *pNewControlPoints, new_knots);

        for (std::size_t dim = 0; dim < Dim; ++dim)
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
        typename TPatchType::DoubleGridFunctionContainerType DoubleGridFunctions_ = pPatch->DoubleGridFunctions();

        typename TPatchType::Array1DGridFunctionContainerType Array1DGridFunctions_ = pPatch->Array1DGridFunctions();

        typename TPatchType::VectorGridFunctionContainerType VectorGridFunctions_ = pPatch->VectorGridFunctions();

        // transfer the grid function
        // here to transfer correctly we apply a two-step process:
        // + firstly the old control values is multiplied with weight to make it weighted control values
        // + secondly the control values will be transferred
        // + the new control values will be divided by the new weight to make it unweighted
        for (typename TPatchType::DoubleGridFunctionContainerType::const_iterator it = DoubleGridFunctions_.begin();
                it != DoubleGridFunctions_.end(); ++it)
        {
            typename StructuredControlGrid<Dim, double>::Pointer pNewDoubleControlGrid = typename StructuredControlGrid<Dim, double>::Pointer(new StructuredControlGrid<Dim, double>(new_size));
            typename StructuredControlGrid<Dim, double>::Pointer pDoubleControlGrid = iga::dynamic_pointer_cast<StructuredControlGrid<Dim, double> >((*it)->pControlGrid());
            if (pDoubleControlGrid == nullptr)
                KRATOS_ERROR << "The cast to StructuredControlGrid is failed.";
            this->ComputeBsplinesDegreeElevation<Dim, double>(*pDoubleControlGrid, *pFESpace, order_increment, *pNewDoubleControlGrid, new_knots);
            pNewDoubleControlGrid->SetName((*it)->pControlGrid()->Name());
            pNewPatch->template CreateGridFunction<double>(pNewDoubleControlGrid);
        }

        for (typename TPatchType::Array1DGridFunctionContainerType::const_iterator it = Array1DGridFunctions_.begin();
                it != Array1DGridFunctions_.end(); ++it)
        {
            if ((*it)->pControlGrid()->Name() == "CONTROL_POINT_COORDINATES") continue;
            typename StructuredControlGrid<Dim, array_1d<double, 3> >::Pointer pNewArray1DControlGrid = typename StructuredControlGrid<Dim, array_1d<double, 3> >::Pointer(new StructuredControlGrid<Dim, array_1d<double, 3> >(new_size));
            typename StructuredControlGrid<Dim, array_1d<double, 3> >::Pointer pArray1DControlGrid = iga::dynamic_pointer_cast<StructuredControlGrid<Dim, array_1d<double, 3> > >((*it)->pControlGrid());
            if (pArray1DControlGrid == nullptr)
                KRATOS_ERROR << "The cast to StructuredControlGrid is failed.";
            this->ComputeBsplinesDegreeElevation<Dim, array_1d<double, 3> >(*pArray1DControlGrid, *pFESpace, order_increment, *pNewArray1DControlGrid, new_knots);
            pNewArray1DControlGrid->SetName((*it)->pControlGrid()->Name());
            pNewPatch->template CreateGridFunction<array_1d<double, 3> >(pNewArray1DControlGrid);
        }

        for (typename TPatchType::VectorGridFunctionContainerType::const_iterator it = VectorGridFunctions_.begin();
                it != VectorGridFunctions_.end(); ++it)
        {
            typename StructuredControlGrid<Dim, Vector>::Pointer pNewVectorControlGrid = typename StructuredControlGrid<Dim, Vector>::Pointer(new StructuredControlGrid<Dim, Vector>(new_size));
            typename StructuredControlGrid<Dim, Vector>::Pointer pVectorControlGrid = iga::dynamic_pointer_cast<StructuredControlGrid<Dim, Vector> >((*it)->pControlGrid());
            if (pVectorControlGrid == nullptr)
                KRATOS_ERROR << "The cast to StructuredControlGrid is failed.";
            this->ComputeBsplinesDegreeElevation<Dim, Vector>(*pVectorControlGrid, *pFESpace, order_increment, *pNewVectorControlGrid, new_knots);
            pNewVectorControlGrid->SetName((*it)->pControlGrid()->Name());
            pNewPatch->template CreateGridFunction<Vector>(pNewVectorControlGrid);
        }

        // mark refined patch
        if (refined_patches.find(pPatch->Id()) == refined_patches.end())
        {
            refined_patches[pPatch->Id()].resize(Dim);
            for (unsigned int i = 0; i < Dim; ++i)
                refined_patches[pPatch->Id()][i] = 0;
        }
        for (unsigned int i = 0; i < Dim; ++i)
        {
            if (order_increment[i] != 0)
            {
                refined_patches[pPatch->Id()][i] = 1;
            }
        }

        for (typename TPatchType::interface_iterator it = pPatch->InterfaceBegin(); it != pPatch->InterfaceEnd(); ++it)
        {
            typename BSplinesPatchInterfaceType::Pointer pInterface = iga::dynamic_pointer_cast<BSplinesPatchInterfaceType>(*it);
            if (pInterface == nullptr)
                KRATOS_ERROR << "The cast to BSplinesPatchInterface is failed";
            typename TPatchType::Pointer pNeighbor = pInterface->pPatch2();

            if (pNeighbor->pFESpace()->Type() != BSplinesFESpaceType::StaticType())
                KRATOS_ERROR << "The FESpace of the neighbor is not BSplinesFESpace";

            // transfer the order increment to neighbors
            std::vector<std::size_t> neib_order_increment(Dim);
            for (unsigned int i = 0; i < Dim; ++i)
                neib_order_increment[i] = 0;

            if constexpr (Dim == 2)
            {
                int dir1 = ParameterDirection<2>::Get_(pInterface->Side1());
                int dir2 = ParameterDirection<2>::Get_(pInterface->Side2());

                neib_order_increment[dir2] = order_increment[dir1];
            }
            else if constexpr (Dim == 3)
            {
                std::vector<int> param_dirs_1 = ParameterDirection<3>::Get(pInterface->Side1());
                std::vector<int> param_dirs_2 = ParameterDirection<3>::Get(pInterface->Side2());

                neib_order_increment[param_dirs_2[ pInterface->LocalParameterMapping(0) ] ] = order_increment[param_dirs_1[0]];
                neib_order_increment[param_dirs_2[ pInterface->LocalParameterMapping(1) ] ] = order_increment[param_dirs_1[1]];
            }

            DegreeElevate<TPatchType>(pNeighbor, refined_patches, neib_order_increment);

            pInterface->SetPatch1(pNewPatch);
            pInterface->SetPatch2(pNeighbor);
            pInterface->pOtherInterface()->SetPatch2(pNewPatch);

            pNewPatch->AddInterface(pInterface);
        }

        // get the parent multipatch
        auto pMultiPatch = pPatch->pParentMultiPatch();

        if (pMultiPatch != nullptr)
        {
            // set the parent multipatch
            pNewPatch->pSetParentMultiPatch(pMultiPatch);

            // remove this patch from multipatch
            pMultiPatch->Patches().erase(pPatch->Id());
        }

        // swap
        #ifdef ENABLE_VERBOSE_REFINEMENT
        std::cout << __FUNCTION__ << ": " << pPatch << " is swapped with ";
        #endif
        pPatch.swap(pNewPatch);
        #ifdef ENABLE_VERBOSE_REFINEMENT
        std::cout << pPatch << std::endl;
        #endif

        if (pMultiPatch != nullptr)
        {
            // replace the corresponding patch in multipatch
            pMultiPatch->Patches().push_back(pPatch);
            pMultiPatch->Patches().Unique();
        }

        #ifdef ENABLE_VERBOSE_REFINEMENT
        std::cout << __FUNCTION__ << " completed for patch " << pPatch->Id() << std::endl;
        #endif
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
#undef ENABLE_VERBOSE_REFINEMENT

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_REFINEMENT_UTILITY_HPP_INCLUDED defined
