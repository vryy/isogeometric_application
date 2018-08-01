//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 25 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_PATCH_UTILITY_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_PATCH_UTILITY_H_INCLUDED

// System includes
#include <vector>
#include <fstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"
#include "custom_utilities/hbsplines/hbsplines_basis_function.h"
#include "custom_utilities/hbsplines/hbsplines_fespace.h"
#include "custom_utilities/grid_function.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/multipatch.h"
#include "custom_utilities/control_grid_utility.h"
#include "isogeometric_application/isogeometric_application.h"

namespace Kratos
{

template<int TDim>
struct HBSplinesPatchUtility_Helper
{
    static typename Patch<TDim>::Pointer CreatePatchFromBSplines(typename Patch<TDim>::Pointer pPatch);
};

/**
Utility to generate hierarchical B-Splines patch.
*/
class HBSplinesPatchUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HBSplinesPatchUtility);

    /// Type definition

    /// Default constructor
    HBSplinesPatchUtility() {}

    /// Destructor
    virtual ~HBSplinesPatchUtility() {}

    /// Create hbsplines patch from bsplines patch
    /// One has to ensure that the B-Splines patch is enumerated properly before transforming to hierarchical B-Splines
    template<int TDim>
    static typename Patch<TDim>::Pointer CreatePatchFromBSplines(typename Patch<TDim>::Pointer pPatch)
    {
        return HBSplinesPatchUtility_Helper<TDim>::CreatePatchFromBSplines(pPatch);
    }

    /// List the boundary basis functions on side
    template<int TDim>
    static void ListBoundaryBfs(std::ostream& rOStream, typename HBSplinesFESpace<TDim>::Pointer pFESpace, const BoundarySide& side)
    {
        typedef typename HBSplinesFESpace<TDim>::bf_iterator bf_iterator;

        rOStream << "Listing of boundary basis function on " << BoundarySideName(side) << " side of hierarchical B-Splines space:" << std::endl;
        rOStream << "<<<<<" << std::endl;

        for (bf_iterator it = pFESpace->bf_begin(); it != pFESpace->bf_end(); ++it)
        {
            if ((*it)->IsOnSide(BOUNDARY_FLAG(side)))
            {
/*                rOStream << "  bf " << (*it)->Id() << ", eq_id: " << (*it)->EquationId() << std::endl;*/
                rOStream << *(*it) << std::endl;
            }
        }

        rOStream << ">>>>>" << std::endl;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "HBSplinesPatchUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

}; // end class HBSplinesPatchUtility

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const HBSplinesPatchUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

template<>
Patch<2>::Pointer HBSplinesPatchUtility_Helper<2>::CreatePatchFromBSplines(typename Patch<2>::Pointer pPatch)
{
    // firstly read the FESpace
    if (pPatch->pFESpace()->Type() != BSplinesFESpace<2>::StaticType())
        KRATOS_THROW_ERROR(std::logic_error, "The input patch is not B-Splines patch", "")

    typename BSplinesFESpace<2>::Pointer pFESpace = boost::dynamic_pointer_cast<BSplinesFESpace<2> >(pPatch->pFESpace());
    if (pFESpace == NULL)
        KRATOS_THROW_ERROR(std::runtime_error, "The cast to BSplinesFESpace is failed.", "")

    // create the hierarchical B-Splines FESpace
    typename HBSplinesFESpace<2>::Pointer pNewFESpace = HBSplinesFESpace<2>::Create();

    for (int dim = 0; dim < 2; ++dim)
    {
        pNewFESpace->SetInfo(dim, pPatch->pFESpace()->Order(dim));
        pNewFESpace->KnotVector(dim) = pFESpace->KnotVector(dim);
    }

    typedef typename BSplinesFESpace<2>::knot_t knot_t;
    typedef typename HBSplinesFESpace<2>::cell_t cell_t;
    typedef typename Patch<2>::ControlPointType ControlPointType;

    std::size_t level = 1;
    const double area_tol = 1.0e-6; // tolerance to accept the nonzero-area cell. We should parameterize it. TODO
    std::size_t number_1 = pFESpace->KnotVector(0).size() - pFESpace->Order(0) - 1;
    std::size_t number_2 = pFESpace->KnotVector(1).size() - pFESpace->Order(1) - 1;

    std::vector<Variable<double>*> double_var_list = pPatch->ExtractVariables<Variable<double> >();
    std::vector<Variable<array_1d<double, 3> >*> array1d_var_list = pPatch->ExtractVariables<Variable<array_1d<double, 3> > >();
    std::vector<Variable<Vector>*> vector_var_list = pPatch->ExtractVariables<Variable<Vector> >();

    std::vector<std::size_t> func_indices = pFESpace->FunctionIndices();

    std::size_t id = 0;
    for(std::size_t j = 0; j < number_2; ++j)
    {
        // create and fill the local knot vector
        std::vector<knot_t> pLocalKnots2;
        for(std::size_t k = 0; k < pFESpace->Order(1) + 2; ++k)
            pLocalKnots2.push_back(pFESpace->KnotVector(1).pKnotAt(j + k));

        for(std::size_t i = 0; i < number_1; ++i)
        {
            // create and fill the local knot vector
            std::vector<knot_t> pLocalKnots1;
            for(std::size_t k = 0; k < pFESpace->Order(0) + 2; ++k)
                pLocalKnots1.push_back(pFESpace->KnotVector(0).pKnotAt(i + k));

            // create the basis function object
            std::size_t i_func = j * number_1 + i;
            std::vector<std::vector<knot_t> > pLocalKnots = {pLocalKnots1, pLocalKnots2};

            const std::size_t& func_id = func_indices[i_func];
            typename HBSplinesBasisFunction<2>::Pointer p_bf = pNewFESpace->CreateBf(++id, level, pLocalKnots);
            p_bf->SetEquationId(func_id);

            // set the boundary information
            if (i == 0) p_bf->AddBoundary(BOUNDARY_FLAG(_BLEFT_));
            else if (i == number_1-1) p_bf->AddBoundary(BOUNDARY_FLAG(_BRIGHT_));

            if (j == 0) p_bf->AddBoundary(BOUNDARY_FLAG(_BBOTTOM_));
            else if (j == number_2-1) p_bf->AddBoundary(BOUNDARY_FLAG(_BTOP_));

            // create the cells for the basis function
            for(std::size_t i1 = 0; i1 < pFESpace->Order(0) + 1; ++i1)
            {
                knot_t pLeft = pFESpace->KnotVector(0).pKnotAt(i + i1);
                knot_t pRight = pFESpace->KnotVector(0).pKnotAt(i + i1 + 1);

                for(std::size_t j1 = 0; j1 < pFESpace->Order(1) + 1; ++j1)
                {
                    knot_t pDown = pFESpace->KnotVector(1).pKnotAt(j + j1);
                    knot_t pUp = pFESpace->KnotVector(1).pKnotAt(j + j1 + 1);

                    // check if the cell domain area is nonzero
                    double area = (pRight->Value() - pLeft->Value()) * (pUp->Value() - pDown->Value());
                    if(fabs(area) > area_tol)
                    {
                        std::vector<knot_t> pKnots = {pLeft, pRight, pDown, pUp};
                        cell_t p_cell = pNewFESpace->pCellManager()->CreateCell(pKnots);
                        p_cell->SetLevel(level);
                        p_bf->AddCell(p_cell);
                        p_cell->AddBf(p_bf);
                    }
                }
            }

            // transfer the control point
            ControlPointType c = pPatch->pControlPointGridFunction()->pControlGrid()->GetData(i_func);
            p_bf->SetValue(CONTROL_POINT, c);

            // transfer other data
            for (std::size_t i = 0; i < double_var_list.size(); ++i)
            {
                GridFunction<2, double>::Pointer pGridFunction = pPatch->pGetGridFunction<Variable<double> >(*double_var_list[i]);
                const double& v = pGridFunction->pControlGrid()->GetData(i_func);
                p_bf->SetValue(*double_var_list[i], v);
            }

            for (std::size_t i = 0; i < array1d_var_list.size(); ++i)
            {
                GridFunction<2, array_1d<double, 3> >::Pointer pGridFunction = pPatch->pGetGridFunction<Variable<array_1d<double, 3> > >(*array1d_var_list[i]);
                const array_1d<double, 3>& v = pGridFunction->pControlGrid()->GetData(i_func);
                p_bf->SetValue(*array1d_var_list[i], v);
            }

            for (std::size_t i = 0; i < vector_var_list.size(); ++i)
            {
                GridFunction<2, Vector>::Pointer pGridFunction = pPatch->pGetGridFunction<Variable<Vector> >(*vector_var_list[i]);
                const Vector& v = pGridFunction->pControlGrid()->GetData(i_func);
                p_bf->SetValue(*vector_var_list[i], v);
            }
        }
    }

    Patch<2>::Pointer pNewPatch = Patch<2>::Create(pPatch->Id(), pNewFESpace);

    // set control points grid function
    typename ControlGrid<ControlPointType>::Pointer pControlPointGrid
        = ControlGridUtility::CreatePointBasedControlGrid<ControlPointType, HBSplinesFESpace<2> >(CONTROL_POINT, pNewFESpace);
    pNewPatch->CreateControlPointGridFunction(pControlPointGrid);

    // set other control grid functions
    for (std::size_t i = 0; i < double_var_list.size(); ++i)
    {
        typename ControlGrid<double>::Pointer pControlGrid
            = ControlGridUtility::CreatePointBasedControlGrid<double, HBSplinesFESpace<2> >(*double_var_list[i], pNewFESpace);
        pNewPatch->CreateGridFunction(*double_var_list[i], pControlGrid);
    }

    for (std::size_t i = 0; i < array1d_var_list.size(); ++i)
    {
        typename ControlGrid<array_1d<double, 3> >::Pointer pControlGrid
            = ControlGridUtility::CreatePointBasedControlGrid<array_1d<double, 3>, HBSplinesFESpace<2> >(*array1d_var_list[i], pNewFESpace);
        pNewPatch->CreateGridFunction(*array1d_var_list[i], pControlGrid);
    }

    for (std::size_t i = 0; i < vector_var_list.size(); ++i)
    {
        typename ControlGrid<Vector>::Pointer pControlGrid
            = ControlGridUtility::CreatePointBasedControlGrid<Vector, HBSplinesFESpace<2> >(*vector_var_list[i], pNewFESpace);
        pNewPatch->CreateGridFunction(*vector_var_list[i], pControlGrid);
    }

    return pNewPatch;
}

template<>
Patch<3>::Pointer HBSplinesPatchUtility_Helper<3>::CreatePatchFromBSplines(typename Patch<3>::Pointer pPatch)
{
    // firstly read the FESpace
    if (pPatch->pFESpace()->Type() != BSplinesFESpace<3>::StaticType())
        KRATOS_THROW_ERROR(std::logic_error, "The input patch is not B-Splines patch", "")

    typename BSplinesFESpace<3>::Pointer pFESpace = boost::dynamic_pointer_cast<BSplinesFESpace<3> >(pPatch->pFESpace());
    if (pFESpace == NULL)
        KRATOS_THROW_ERROR(std::runtime_error, "The cast to BSplinesFESpace is failed.", "")

    // create the hierarchical B-Splines FESpace
    typename HBSplinesFESpace<3>::Pointer pNewFESpace = HBSplinesFESpace<3>::Create();

    for (int dim = 0; dim < 3; ++dim)
        pNewFESpace->SetInfo(dim, pPatch->pFESpace()->Order(dim));

    typedef typename BSplinesFESpace<3>::knot_t knot_t;
    typedef typename HBSplinesFESpace<3>::cell_t cell_t;
    typedef typename Patch<3>::ControlPointType ControlPointType;

    std::size_t level = 1;
    const double area_tol = 1.0e-6; // tolerance to accept the nonzero-area cell. We should parameterize it. TODO
    std::size_t number_1 = pFESpace->KnotVector(0).size() - pFESpace->Order(0) - 1;
    std::size_t number_2 = pFESpace->KnotVector(1).size() - pFESpace->Order(1) - 1;
    std::size_t number_3 = pFESpace->KnotVector(2).size() - pFESpace->Order(2) - 1;

    std::vector<Variable<double>*> double_var_list = pPatch->ExtractVariables<Variable<double> >();
    std::vector<Variable<array_1d<double, 3> >*> array1d_var_list = pPatch->ExtractVariables<Variable<array_1d<double, 3> > >();
    std::vector<Variable<Vector>*> vector_var_list = pPatch->ExtractVariables<Variable<Vector> >();

    std::vector<std::size_t> func_indices = pFESpace->FunctionIndices();

    std::size_t id = 0;
    for(std::size_t l = 0; l < number_3; ++l)
    {
        // create and fill the local knot vector
        std::vector<knot_t> pLocalKnots3;
        for(std::size_t k = 0; k < pFESpace->KnotVector(2).size() + 2; ++k)
            pLocalKnots3.push_back(pFESpace->KnotVector(2).pKnotAt(l + k));

        for(std::size_t j = 0; j < number_2; ++j)
        {
            // create and fill the local knot vector
            std::vector<knot_t> pLocalKnots2;
            for(std::size_t k = 0; k < pFESpace->KnotVector(1).size() + 2; ++k)
                pLocalKnots2.push_back(pFESpace->KnotVector(1).pKnotAt(j + k));

            for(std::size_t i = 0; i < number_1; ++i)
            {
                // create and fill the local knot vector
                std::vector<knot_t> pLocalKnots1;
                for(std::size_t k = 0; k < pFESpace->KnotVector(0).size() + 2; ++k)
                    pLocalKnots1.push_back(pFESpace->KnotVector(0).pKnotAt(i + k));

                // create the basis function object
                std::size_t i_func = (l * number_2 + j) * number_1 + i;
                std::vector<std::vector<knot_t> > pLocalKnots = {pLocalKnots1, pLocalKnots2, pLocalKnots3};

                const std::size_t& func_id = func_indices[i_func];
                typename HBSplinesBasisFunction<3>::Pointer p_bf = pNewFESpace->CreateBf(++id, level, pLocalKnots);
                p_bf->SetEquationId(func_id);

                // set the boundary information
                if (i == 0) p_bf->AddBoundary(BOUNDARY_FLAG(_BLEFT_));
                else if (i == number_1-1) p_bf->AddBoundary(BOUNDARY_FLAG(_BRIGHT_));

                if (j == 0) p_bf->AddBoundary(BOUNDARY_FLAG(_BFRONT_));
                else if (j == number_2-1) p_bf->AddBoundary(BOUNDARY_FLAG(_BBACK_));

                if (l == 0) p_bf->AddBoundary(BOUNDARY_FLAG(_BBOTTOM_));
                else if (l == number_3-1) p_bf->AddBoundary(BOUNDARY_FLAG(_BTOP_));

                // create the cells for the basis function
                for(std::size_t i1 = 0; i1 < pFESpace->Order(0) + 1; ++i1)
                {
                    knot_t pLeft = pFESpace->KnotVector(0).pKnotAt(i + i1);
                    knot_t pRight = pFESpace->KnotVector(0).pKnotAt(i + i1 + 1);

                    for(std::size_t j1 = 0; j1 < pFESpace->Order(1) + 1; ++j1)
                    {
                        knot_t pDown = pFESpace->KnotVector(1).pKnotAt(j + j1);
                        knot_t pUp = pFESpace->KnotVector(1).pKnotAt(j + j1 + 1);

                        for(std::size_t l1 = 0; l1 < pFESpace->Order(2) + 1; ++l1)
                        {
                            knot_t pBelow = pFESpace->KnotVector(2).pKnotAt(l + l1);
                            knot_t pAbove = pFESpace->KnotVector(2).pKnotAt(l + l1 + 1);

                            // check if the cell domain area is nonzero
                            double area = (pRight->Value() - pLeft->Value()) * (pUp->Value() - pDown->Value()) * (pAbove->Value() - pBelow->Value());
                            if(fabs(area) > area_tol)
                            {
                                std::vector<knot_t> pKnots = {pLeft, pRight, pDown, pUp, pBelow, pAbove};
                                cell_t p_cell = pNewFESpace->pCellManager()->CreateCell(pKnots);
                                p_cell->SetLevel(level);
                                p_bf->AddCell(p_cell);
                                p_cell->AddBf(p_bf);
                            }
                        }
                    }
                }
            }
        }
    }

    typename Patch<3>::Pointer pNewPatch = Patch<3>::Create(pPatch->Id(), pNewFESpace);

    // set control points grid
    typename ControlGrid<ControlPointType>::Pointer pControlPointGrid
        = ControlGridUtility::CreatePointBasedControlGrid<ControlPointType, HBSplinesFESpace<3> >(CONTROL_POINT, pNewFESpace);
    pNewPatch->CreateControlPointGridFunction(pControlPointGrid);

    // set other control grid
    for (std::size_t i = 0; i < double_var_list.size(); ++i)
    {
        typename ControlGrid<double>::Pointer pControlGrid
            = ControlGridUtility::CreatePointBasedControlGrid<double, HBSplinesFESpace<3> >(*double_var_list[i], pNewFESpace);
        pNewPatch->CreateGridFunction(*double_var_list[i], pControlGrid);
    }

    for (std::size_t i = 0; i < array1d_var_list.size(); ++i)
    {
        typename ControlGrid<array_1d<double, 3> >::Pointer pControlGrid
            = ControlGridUtility::CreatePointBasedControlGrid<array_1d<double, 3>, HBSplinesFESpace<3> >(*array1d_var_list[i], pNewFESpace);
        pNewPatch->CreateGridFunction(*array1d_var_list[i], pControlGrid);
    }

    for (std::size_t i = 0; i < vector_var_list.size(); ++i)
    {
        typename ControlGrid<Vector>::Pointer pControlGrid
            = ControlGridUtility::CreatePointBasedControlGrid<Vector, HBSplinesFESpace<3> >(*vector_var_list[i], pNewFESpace);
        pNewPatch->CreateGridFunction(*vector_var_list[i], pControlGrid);
    }

    return pNewPatch;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_PATCH_UTILITY_H_INCLUDED defined

