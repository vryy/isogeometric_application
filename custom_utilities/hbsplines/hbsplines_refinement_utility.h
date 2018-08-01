//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 28 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_REFINEMENT_UTILITY_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_REFINEMENT_UTILITY_H_INCLUDED

// System includes
#include <vector>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "custom_utilities/bspline_utils.h"
#include "custom_utilities/nurbs/knot_array_1d.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/grid_function.h"
#include "custom_utilities/patch.h"
#include "custom_utilities/hbsplines/hbsplines_fespace.h"


namespace Kratos
{


template<int TDim>
struct HBSplinesRefinementUtility_Helper
{
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;
    typedef typename HBSplinesFESpace<TDim>::bf_t bf_t;

    static void Refine(typename Patch<TDim>::Pointer pPatch, const std::size_t& Id, const int& echo_level);

    static std::pair<std::vector<std::size_t>, std::vector<bf_t> > Refine(typename Patch<TDim>::Pointer pPatch,
            typename HBSplinesFESpace<TDim>::bf_t p_bf, std::set<std::size_t>& refined_patches, const int& echo_level);

    static void RefineWindow(typename Patch<TDim>::Pointer pPatch, const std::vector<std::vector<double> >& window, const int& echo_level);

    static void LinearDependencyRefine(typename Patch<TDim>::Pointer pPatch, const std::size_t& refine_cycle, const int& echo_level);
};


/**
Class accounts for linear independent refinement of a single hierarchical B-Splines patch
 */
class HBSplinesRefinementUtility
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HBSplinesRefinementUtility);

    /// Default constructor
    HBSplinesRefinementUtility() {}

    /// Destructor
    virtual ~HBSplinesRefinementUtility() {}

    /// Refine a single B-Splines basis function
    template<int TDim>
    static void Refine(typename Patch<TDim>::Pointer pPatch, const std::size_t& Id, const int& echo_level)
    {
        HBSplinesRefinementUtility_Helper<TDim>::Refine(pPatch, Id, echo_level);
    }


    /// Refine the basis functions in a region
    template<int TDim>
    static void RefineWindow(typename Patch<TDim>::Pointer pPatch, const std::vector<std::vector<double> >& window, const int& echo_level)
    {
        HBSplinesRefinementUtility_Helper<TDim>::RefineWindow(pPatch, window, echo_level);
    }


    /// Perform additional refinement to ensure linear independence
    /// In this algorithm, every bf in each level will be checked. If the support domain of a bf contained in the domain_manager of that level, this bf will be refined. According to the paper of Vuong et al, this will produce a linear independent bases.
    template<int TDim>
    static void LinearDependencyRefine(typename Patch<TDim>::Pointer pPatch, const std::size_t& refine_cycle, const int& echo_level)
    {
        HBSplinesRefinementUtility_Helper<TDim>::LinearDependencyRefine(pPatch, refine_cycle, echo_level);
    }


    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "HBSplinesRefinementUtility";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const HBSplinesRefinementUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

template<int TDim>
inline void HBSplinesRefinementUtility_Helper<TDim>::Refine(typename Patch<TDim>::Pointer pPatch, const std::size_t& Id, const int& echo_level)
{
    typedef typename HBSplinesFESpace<TDim>::bf_t bf_t;
    typedef typename HBSplinesFESpace<TDim>::bf_container_t bf_container_t;

    if (pPatch->pFESpace()->Type() != HBSplinesFESpace<TDim>::StaticType())
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "only support the hierarchical B-Splines patch")

    // extract the hierarchical B-Splines space
    typename HBSplinesFESpace<TDim>::Pointer pFESpace = boost::dynamic_pointer_cast<HBSplinesFESpace<TDim> >(pPatch->pFESpace());
    if (pFESpace == NULL)
        KRATOS_THROW_ERROR(std::runtime_error, "The cast to HBSplinesFESpace is failed.", "")

    // get the correct basis function
    bf_t p_bf;
    bool found = false;
    for(typename bf_container_t::iterator it = pFESpace->bf_begin(); it != pFESpace->bf_end(); ++it)
    {
        if((*it)->Id() == Id)
        {
            p_bf = *it;
            found = true;
        }
    }

    if(!found) return;

    // does not refine if maximum level is reached
    if(p_bf->Level() == pFESpace->MaxLevel())
    {
        std::cout << "Maximum level is reached, basis function " << p_bf->Id() << " is skipped" << std::endl;
        return;
    }

    // refine the patch
    std::set<std::size_t> refined_patches;
    Refine(pPatch, p_bf, refined_patches, echo_level);
}

template<int TDim>
std::pair<std::vector<std::size_t>, std::vector<typename HBSplinesFESpace<TDim>::bf_t> > HBSplinesRefinementUtility_Helper<TDim>::Refine(
        typename Patch<TDim>::Pointer pPatch, typename HBSplinesFESpace<TDim>::bf_t p_bf,
        std::set<std::size_t>& refined_patches, const int& echo_level)
{
    // Type definitions
    typedef typename HBSplinesFESpace<TDim>::bf_t bf_t;
    typedef typename HBSplinesFESpace<TDim>::bf_container_t bf_container_t;
    typedef typename HBSplinesFESpace<TDim>::CellType CellType;
    typedef typename HBSplinesFESpace<TDim>::cell_t cell_t;
    typedef typename HBSplinesFESpace<TDim>::cell_container_t cell_container_t;
    typedef typename Patch<TDim>::ControlPointType ControlPointType;

    if (std::find(refined_patches.begin(), refined_patches.end(), pPatch->Id()) != refined_patches.end())
    {
        std::vector<std::size_t> aux1;
        std::vector<bf_t> aux2;

        return std::make_pair(aux1, aux2);
    }
    else
    {
        refined_patches.insert(pPatch->Id());
    }

    #ifdef ENABLE_PROFILING
    double start = OpenMPUtils::GetCurrentTime();
    #endif

    // save the equation_id
    std::size_t equation_id = p_bf->EquationId();

    // extract the hierarchical B-Splines space
    typename HBSplinesFESpace<TDim>::Pointer pFESpace = boost::dynamic_pointer_cast<HBSplinesFESpace<TDim> >(pPatch->pFESpace());
    if (pFESpace == NULL)
        KRATOS_THROW_ERROR(std::runtime_error, "The cast to HBSplinesFESpace is failed.", "")

    // get the list of variables in the patch
    std::vector<Variable<double>*> double_variables = pPatch->template ExtractVariables<Variable<double> >();
    std::vector<Variable<array_1d<double, 3> >*> array_1d_variables = pPatch->template ExtractVariables<Variable<array_1d<double, 3> > >();
    std::vector<Variable<Vector>*> vector_variables = pPatch->template ExtractVariables<Variable<Vector> >();

    /* create a list of basis function in the next level representing this basis function */
    // create a list of new knots
    double cell_tol = pFESpace->pCellManager()->GetTolerance(); // tolerance to accept the nonzero-area cell. We should parameterize it.
    std::vector<std::vector<knot_t> > pnew_local_knots(TDim);
    std::vector<std::vector<double> > ins_knots(TDim);
    for(unsigned int dim = 0; dim < TDim; ++dim)
    {
        const std::vector<knot_t>& pLocalKnots = p_bf->LocalKnots(dim);

        for(std::vector<knot_t>::const_iterator it = pLocalKnots.begin(); it != pLocalKnots.end(); ++it)
        {
            pnew_local_knots[dim].push_back(*it);

            std::vector<knot_t>::const_iterator it2 = it + 1;
            if(it2 != pLocalKnots.end())
            {
                if(fabs((*it2)->Value() - (*it)->Value()) > cell_tol)
                {
                    // now we just add the middle one, but in the general we can add arbitrary values
                    // TODO find the way to generalize this or parameterize this
                    double ins_knot = 0.5 * ((*it)->Value() + (*it2)->Value());
                    knot_t p_new_knot;
                    p_new_knot = pFESpace->KnotVector(dim).pCreateUniqueKnot(ins_knot, cell_tol);
                    pnew_local_knots[dim].push_back(p_new_knot);
                    ins_knots[dim].push_back(p_new_knot->Value());
                }
            }
        }
    }

    /* compute the refinement coefficients */
    Vector RefinedCoeffs;
    std::vector<std::vector<double> > local_knots(TDim);
    for(std::size_t dim = 0; dim < TDim; ++dim)
        p_bf->LocalKnots(dim, local_knots[dim]);

    std::vector<std::vector<double> > new_knots(TDim);
    if (TDim == 2)
    {
        BSplineUtils::ComputeBsplinesKnotInsertionCoefficients2DLocal(RefinedCoeffs,
            new_knots[0], new_knots[1],
            pFESpace->Order(0), pFESpace->Order(1),
            local_knots[0], local_knots[1],
            ins_knots[0], ins_knots[1]);
    }
    else if (TDim == 3)
    {
        BSplineUtils::ComputeBsplinesKnotInsertionCoefficients3DLocal(RefinedCoeffs,
            new_knots[0], new_knots[1], new_knots[2],
            pFESpace->Order(0), pFESpace->Order(1), pFESpace->Order(2),
            local_knots[0], local_knots[1], local_knots[2],
            ins_knots[0], ins_knots[1], ins_knots[2]);
    }
    // KRATOS_WATCH(RefinedCoeffs)

    #ifdef ENABLE_PROFILING
    double time_1 = OpenMPUtils::GetCurrentTime() - start;
    start = OpenMPUtils::GetCurrentTime();
    #endif

    /* create new basis functions */
    unsigned int next_level = p_bf->Level() + 1;
    if(next_level > pFESpace->LastLevel()) pFESpace->SetLastLevel(next_level);
    std::size_t last_id = pFESpace->LastId();
    typename cell_container_t::Pointer pnew_cells;

    std::vector<std::size_t> numbers(TDim);
    std::vector<bf_t> pnew_bfs;

    // start to enumerate from the last equation id in the multipatch
    // we always assign an incremental equation_id for the new refined bfs, so that the bfs on the boundary will automaticall match
    std::size_t starting_id;
    if (pPatch->pParentMultiPatch() != NULL)
    {
        starting_id = pPatch->pParentMultiPatch()->GetLastEquationId();
    }
    else
    {
        starting_id = pPatch->pFESpace()->GetLastEquationId();
    }

    if (TDim == 2)
    {
        pnew_cells = typename cell_container_t::Pointer(new CellManager2D<CellType>());

        numbers[0] = pnew_local_knots[0].size() - pFESpace->Order(0) - 1;
        numbers[1] = pnew_local_knots[1].size() - pFESpace->Order(1) - 1;

        for(std::size_t j = 0; j < numbers[1]; ++j)
        {
            // create and fill the local knot vector
            std::vector<knot_t> pLocalKnots2;
            for(std::size_t k = 0; k < pFESpace->Order(1) + 2; ++k)
                pLocalKnots2.push_back(pnew_local_knots[1][j + k]);

            for(std::size_t i = 0; i < numbers[0]; ++i)
            {
                // create and fill the local knot vector
                std::vector<knot_t> pLocalKnots1;
                for(std::size_t k = 0; k < pFESpace->Order(0) + 2; ++k)
                    pLocalKnots1.push_back(pnew_local_knots[0][i + k]);

                // create the basis function object
                std::vector<std::vector<knot_t> > pLocalKnots = {pLocalKnots1, pLocalKnots2};
                bf_t pnew_bf = pFESpace->CreateBf(last_id+1, next_level, pLocalKnots);
                pnew_bfs.push_back(pnew_bf);
                if (pnew_bf->Id() == last_id+1) ++last_id;

                std::size_t i_func = j * numbers[0] + i;

                // set the boundary information
                if (p_bf->IsOnSide(BOUNDARY_FLAG(_BLEFT_)))
                    if (i == 0)
                        pnew_bfs[i_func]->AddBoundary(BOUNDARY_FLAG(_BLEFT_));

                if (p_bf->IsOnSide(BOUNDARY_FLAG(_BRIGHT_)))
                    if (i == numbers[0]-1)
                        pnew_bfs[i_func]->AddBoundary(BOUNDARY_FLAG(_BRIGHT_));

                if (p_bf->IsOnSide(BOUNDARY_FLAG(_BBOTTOM_)))
                    if (j == 0)
                        pnew_bfs[i_func]->AddBoundary(BOUNDARY_FLAG(_BBOTTOM_));

                if (p_bf->IsOnSide(BOUNDARY_FLAG(_BTOP_)))
                    if (j == numbers[1]-1)
                        pnew_bfs[i_func]->AddBoundary(BOUNDARY_FLAG(_BTOP_));

                // assign new equation id
                pnew_bf->SetEquationId(++starting_id);
                if((echo_level & ECHO_REFINEMENT) == ECHO_REFINEMENT)
                {
                    std::cout << "new bf " << pnew_bf->Id() << " is assigned eq_id = " << pnew_bf->EquationId() << std::endl;
                }

                // transfer the control point information
                ControlPointType oldC = p_bf->GetValue(CONTROL_POINT);
                ControlPointType newC = pnew_bf->GetValue(CONTROL_POINT);
                newC += RefinedCoeffs[i_func] * oldC;
                pnew_bf->SetValue(CONTROL_POINT, newC);

                // transfer other control values from p_bf to pnew_bf
                for (std::size_t i = 0; i < double_variables.size(); ++i)
                {
                    double old_value = p_bf->GetValue(*double_variables[i]);
                    double new_value = pnew_bf->GetValue(*double_variables[i]);
                    new_value += RefinedCoeffs[i_func] * old_value;
                    pnew_bf->SetValue(*double_variables[i], new_value);
                }

                for (std::size_t i = 0; i < array_1d_variables.size(); ++i)
                {
                    array_1d<double, 3> old_value = p_bf->GetValue(*array_1d_variables[i]);
                    array_1d<double, 3> new_value = pnew_bf->GetValue(*array_1d_variables[i]);
                    new_value += RefinedCoeffs[i_func] * old_value;
                    pnew_bf->SetValue(*array_1d_variables[i], new_value);
                }

                for (std::size_t i = 0; i < vector_variables.size(); ++i)
                {
                    Vector old_value = p_bf->GetValue(*vector_variables[i]);
                    Vector new_value = pnew_bf->GetValue(*vector_variables[i]);
                    new_value += RefinedCoeffs[i_func] * old_value;
                    pnew_bf->SetValue(*vector_variables[i], new_value);
                }

                // create the cells for the basis function
                for(std::size_t i1 = 0; i1 < pFESpace->Order(0) + 1; ++i1)
                {
                    knot_t pLeft = pnew_local_knots[0][i + i1];
                    knot_t pRight = pnew_local_knots[0][i + i1 + 1];
                    for(std::size_t j1 = 0; j1 < pFESpace->Order(1) + 1; ++j1)
                    {
                        knot_t pDown = pnew_local_knots[1][j + j1];
                        knot_t pUp = pnew_local_knots[1][j + j1 + 1];

                        // check if the cell domain area is nonzero
                        double area = (pRight->Value() - pLeft->Value()) * (pUp->Value() - pDown->Value());
                        if(sqrt(fabs(area)) > cell_tol)
                        {
                            std::vector<knot_t> pKnots = {pLeft, pRight, pDown, pUp};
                            cell_t pnew_cell = pFESpace->pCellManager()->CreateCell(pKnots);
                            pnew_cell->SetLevel(next_level);
                            pnew_bf->AddCell(pnew_cell);
                            pnew_cell->AddBf(pnew_bf);
                            pnew_cells->insert(pnew_cell);
                        }
                    }
                }
            }
        }
    }
    else if (TDim == 3)
    {
        pnew_cells = typename cell_container_t::Pointer(new CellManager3D<CellType>());

        numbers[0] = pnew_local_knots[0].size() - pFESpace->Order(0) - 1;
        numbers[1] = pnew_local_knots[1].size() - pFESpace->Order(1) - 1;
        numbers[2] = pnew_local_knots[2].size() - pFESpace->Order(2) - 1;

        for(std::size_t l = 0; l < numbers[2]; ++l)
        {
            // create and fill the local knot vector
            std::vector<knot_t> pLocalKnots3;
            for(std::size_t k = 0; k < pFESpace->Order(2) + 2; ++k)
                pLocalKnots3.push_back(pnew_local_knots[2][l + k]);

            for(std::size_t j = 0; j < numbers[1]; ++j)
            {
                // create and fill the local knot vector
                std::vector<knot_t> pLocalKnots2;
                for(std::size_t k = 0; k < pFESpace->Order(1) + 2; ++k)
                    pLocalKnots2.push_back(pnew_local_knots[1][j + k]);

                for(std::size_t i = 0; i < numbers[0]; ++i)
                {
                    // create and fill the local knot vector
                    std::vector<knot_t> pLocalKnots1;
                    for(std::size_t k = 0; k < pFESpace->Order(0) + 2; ++k)
                        pLocalKnots1.push_back(pnew_local_knots[0][i + k]);

                    // create the basis function object
                    std::vector<std::vector<knot_t> > pLocalKnots = {pLocalKnots1, pLocalKnots2, pLocalKnots3};
                    bf_t pnew_bf = pFESpace->CreateBf(last_id+1, next_level, pLocalKnots);
                    pnew_bfs.push_back(pnew_bf);
                    if (pnew_bf->Id() == last_id+1) ++last_id;

                    // update the coordinates
                    std::size_t i_func = (l * numbers[1] + j) * numbers[0] + i;

                    // set the boundary information
                    if (p_bf->IsOnSide(BOUNDARY_FLAG(_BLEFT_)))
                        if (i == 0)
                            pnew_bfs[i_func]->AddBoundary(BOUNDARY_FLAG(_BLEFT_));

                    if (p_bf->IsOnSide(BOUNDARY_FLAG(_BRIGHT_)))
                        if (i == numbers[0]-1)
                            pnew_bfs[i_func]->AddBoundary(BOUNDARY_FLAG(_BRIGHT_));

                    if (p_bf->IsOnSide(BOUNDARY_FLAG(_BFRONT_)))
                        if (j == 0)
                            pnew_bfs[i_func]->AddBoundary(BOUNDARY_FLAG(_BFRONT_));

                    if (p_bf->IsOnSide(BOUNDARY_FLAG(_BBACK_)))
                        if (j == numbers[1]-1)
                            pnew_bfs[i_func]->AddBoundary(BOUNDARY_FLAG(_BBACK_));

                    if (p_bf->IsOnSide(BOUNDARY_FLAG(_BBOTTOM_)))
                        if (l == 0)
                            pnew_bfs[i_func]->AddBoundary(BOUNDARY_FLAG(_BBOTTOM_));

                    if (p_bf->IsOnSide(BOUNDARY_FLAG(_BTOP_)))
                        if (l == numbers[2]-1)
                            pnew_bfs[i_func]->AddBoundary(BOUNDARY_FLAG(_BTOP_));

                    // assign new equation id
                    pnew_bf->SetEquationId(++starting_id);

                    // transfer the control point information
                    ControlPointType oldC = p_bf->GetValue(CONTROL_POINT);
                    ControlPointType newC = pnew_bf->GetValue(CONTROL_POINT);
                    newC += RefinedCoeffs[i_func] * oldC;
                    pnew_bf->SetValue(CONTROL_POINT, newC);

                    // transfer other control values from p_bf to pnew_bf
                    for (std::size_t i = 0; i < double_variables.size(); ++i)
                    {
                        double old_value = p_bf->GetValue(*double_variables[i]);
                        double new_value = pnew_bf->GetValue(*double_variables[i]);
                        new_value += RefinedCoeffs[i_func] * old_value;
                        pnew_bf->SetValue(*double_variables[i], new_value);
                    }

                    for (std::size_t i = 0; i < array_1d_variables.size(); ++i)
                    {
                        array_1d<double, 3> old_value = p_bf->GetValue(*array_1d_variables[i]);
                        array_1d<double, 3> new_value = pnew_bf->GetValue(*array_1d_variables[i]);
                        new_value += RefinedCoeffs[i_func] * old_value;
                        pnew_bf->SetValue(*array_1d_variables[i], new_value);
                    }

                    for (std::size_t i = 0; i < vector_variables.size(); ++i)
                    {
                        Vector old_value = p_bf->GetValue(*vector_variables[i]);
                        Vector new_value = pnew_bf->GetValue(*vector_variables[i]);
                        new_value += RefinedCoeffs[i_func] * old_value;
                        pnew_bf->SetValue(*vector_variables[i], new_value);
                    }

                    // create the cells for the basis function
                    for(std::size_t i1 = 0; i1 < pFESpace->Order(0) + 1; ++i1)
                    {
                        knot_t pLeft = pnew_local_knots[0][i + i1];
                        knot_t pRight = pnew_local_knots[0][i + i1 + 1];

                        for(std::size_t j1 = 0; j1 < pFESpace->Order(1) + 1; ++j1)
                        {
                            knot_t pDown = pnew_local_knots[1][j + j1];
                            knot_t pUp = pnew_local_knots[1][j + j1 + 1];

                            for(std::size_t l1 = 0; l1 < pFESpace->Order(2) + 1; ++l1)
                            {
                                knot_t pBelow = pnew_local_knots[2][l + l1];
                                knot_t pAbove = pnew_local_knots[2][l + l1 + 1];

                                // check if the cell domain volume is nonzero
                                double volume = (pRight->Value() - pLeft->Value()) * (pUp->Value() - pDown->Value()) * (pAbove->Value() - pBelow->Value());
                                if(pow(fabs(volume), 1.0/3) > cell_tol)
                                {
                                    std::vector<knot_t> pKnots = {pLeft, pRight, pDown, pUp, pBelow, pAbove};
                                    cell_t pnew_cell = pFESpace->pCellManager()->CreateCell(pKnots);
                                    pnew_cell->SetLevel(next_level);
                                    pnew_bf->AddCell(pnew_cell);
                                    pnew_cell->AddBf(pnew_bf);
                                    pnew_cells->insert(pnew_cell);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // update the weight information for all the grid functions (except the control point grid function)
    std::vector<double> Weights = pFESpace->GetWeights();

    typename Patch<TDim>::DoubleGridFunctionContainerType DoubleGridFunctions_ = pPatch->DoubleGridFunctions();
    for (typename Patch<TDim>::DoubleGridFunctionContainerType::iterator it = DoubleGridFunctions_.begin();
            it != DoubleGridFunctions_.end(); ++it)
    {
        typename WeightedFESpace<TDim>::Pointer pThisFESpace = boost::dynamic_pointer_cast<WeightedFESpace<TDim> >((*it)->pFESpace());
        if (pThisFESpace == NULL)
            KRATOS_THROW_ERROR(std::runtime_error, "The cast to WeightedFESpace is failed.", "")
        pThisFESpace->SetWeights(Weights);
    }

    typename Patch<TDim>::Array1DGridFunctionContainerType Array1DGridFunctions_ = pPatch->Array1DGridFunctions();
    for (typename Patch<TDim>::Array1DGridFunctionContainerType::iterator it = Array1DGridFunctions_.begin();
            it != Array1DGridFunctions_.end(); ++it)
    {
        typename WeightedFESpace<TDim>::Pointer pThisFESpace = boost::dynamic_pointer_cast<WeightedFESpace<TDim> >((*it)->pFESpace());
        if (pThisFESpace == NULL)
            KRATOS_THROW_ERROR(std::runtime_error, "The cast to WeightedFESpace is failed.", "")
        pThisFESpace->SetWeights(Weights);
    }

    typename Patch<TDim>::VectorGridFunctionContainerType VectorGridFunctions_ = pPatch->VectorGridFunctions();
    for (typename Patch<TDim>::VectorGridFunctionContainerType::iterator it = VectorGridFunctions_.begin();
            it != VectorGridFunctions_.end(); ++it)
    {
        typename WeightedFESpace<TDim>::Pointer pThisFESpace = boost::dynamic_pointer_cast<WeightedFESpace<TDim> >((*it)->pFESpace());
        if (pThisFESpace == NULL)
            KRATOS_THROW_ERROR(std::runtime_error, "The cast to WeightedFESpace is failed.", "")
        pThisFESpace->SetWeights(Weights);
    }

    #ifdef ENABLE_PROFILING
    double time_2 = OpenMPUtils::GetCurrentTime() - start;
    start = OpenMPUtils::GetCurrentTime();
    #endif

    /* remove the cells of the old basis function (remove only the cell in the current level) */
    typename cell_container_t::Pointer pcells_to_remove;
    if (TDim == 2)
        pcells_to_remove = typename cell_container_t::Pointer(new CellManager2D<CellType>());
    else if (TDim == 3)
        pcells_to_remove = typename cell_container_t::Pointer(new CellManager3D<CellType>());

    // firstly we check if the cell c of the current bf in the current level cover any sub-cells. Then the sub-cell includes all bfs of the cell c.
    for(typename HBSplinesFESpace<TDim>::BasisFunctionType::cell_iterator it_cell = p_bf->cell_begin(); it_cell != p_bf->cell_end(); ++it_cell)
    {
        if((*it_cell)->Level() == p_bf->Level())
        {
            for(typename cell_container_t::iterator it_subcell = pnew_cells->begin(); it_subcell != pnew_cells->end(); ++it_subcell)
            {
                if((*it_subcell)->template IsCovered<TDim>(*it_cell))
                {
                    for(typename CellType::bf_iterator it_bf = (*it_cell)->bf_begin(); it_bf != (*it_cell)->bf_end(); ++it_bf)
                    {
                        (*it_subcell)->AddBf(*it_bf);
                        (*it_bf)->AddCell(*it_subcell);
                    }
                }
            }

            // mark to remove the old cell
            pcells_to_remove->insert(*it_cell);
        }
    }

    // secondly, it happens that new cell c cover several existing cells. In this case cell c must be removed, and its bfs will be transferred to sub-cells.
    for(typename cell_container_t::iterator it_cell = pnew_cells->begin(); it_cell != pnew_cells->end(); ++it_cell)
    {
        std::vector<cell_t> p_cells = pFESpace->pCellManager()->GetCells(*it_cell);
        if(p_cells.size() > 0)
        {
            if((echo_level & ECHO_REFINEMENT) == ECHO_REFINEMENT)
            {
                std::cout << "cell " << (*it_cell)->Id() << " is detected to contain some smaller cells:";
                for(std::size_t i = 0; i < p_cells.size(); ++i)
                    std::cout << " " << p_cells[i]->Id();
                std::cout << std::endl;
            }

            pcells_to_remove->insert(*it_cell);
            for(std::size_t i = 0; i < p_cells.size(); ++i)
            {
                for(typename CellType::bf_iterator it_bf = (*it_cell)->bf_begin(); it_bf != (*it_cell)->bf_end(); ++it_bf)
                {
                    p_cells[i]->AddBf(*it_bf);
                    (*it_bf)->AddCell(p_cells[i]);
                }
            }
        }
    }

    /* remove the cells from the previous step */
    for(typename cell_container_t::iterator it_cell = pcells_to_remove->begin(); it_cell != pcells_to_remove->end(); ++it_cell)
    {
        pFESpace->pCellManager()->erase(*it_cell);
        for(typename bf_container_t::iterator it_bf = pFESpace->bf_begin(); it_bf != pFESpace->bf_end(); ++it_bf)
            (*it_bf)->RemoveCell(*it_cell);
    }

    /* remove the basis function from all the cells */
    for(typename cell_container_t::iterator it_cell = pFESpace->pCellManager()->begin(); it_cell != pFESpace->pCellManager()->end(); ++it_cell)
        (*it_cell)->RemoveBf(p_bf);

    /* remove the old basis function */
    pFESpace->RemoveBf(p_bf);

    // TODO check if pFESpace->pCellManager()->CollapseCells() can help to further remove the overlapping cells

    #ifdef ENABLE_PROFILING
    double time_3 = OpenMPUtils::GetCurrentTime() - start;
    start = OpenMPUtils::GetCurrentTime();
    #endif

    /* Debug part: to be removed when the code is stable
    for(typename cell_container_t::iterator it_cell = pFESpace->pCellManager()->begin(); it_cell != pFESpace->pCellManager()->end(); ++it_cell)
    {
        std::cout << "cell " << (*it_cell)->Id() << " supports:";
        for(typename CellType::bf_iterator it_bf = (*it_cell)->bf_begin(); it_bf != (*it_cell)->bf_end(); ++it_bf)
            std::cout << " " << (*it_bf)->Id();
        std::cout << std::endl;
    }

    for(typename bf_container_t::iterator it_bf = mBasisFuncs.begin(); it_bf != mBasisFuncs.end(); ++it_bf)
    {
        std::cout << "bf " << (*it_bf)->Id() << " contains cell:";
        for(DeprecatedHBBasisFunction::cell_iterator it_cell = (*it_bf)->cell_begin(); it_cell != (*it_bf)->cell_end(); ++it_cell)
            std::cout << " " << (*it_cell)->Id();
        std::cout << std::endl;
    }
    */

    pFESpace->RecordRefinementHistory(p_bf->Id());
    if((echo_level & ECHO_REFINEMENT) == ECHO_REFINEMENT)
    {
        std::cout << "Refine patch " << pPatch->Id() << ", bf " << p_bf->Id() << ", eq_id " << p_bf->EquationId() << " completed" << std::endl;
        #ifdef ENABLE_PROFILING
        std::cout << " Time to compute the refinement coefficients: " << time_1 << " s" << std::endl;
        std::cout << " Time to create new cells and new bfs: " << time_2 << " s" << std::endl;
        std::cout << " Time to clean up: " << time_3 << " s" << std::endl;
        #endif
    }

    // refine also the neighbors
    for (std::size_t i = 0; i < pPatch->NumberOfInterfaces(); ++i)
    {
        typename PatchInterface<TDim>::Pointer pInterface = pPatch->pInterface(i);

        typename Patch<TDim>::Pointer pNeighborPatch = pInterface->pPatch2();

        // extract the hierarchical B-Splines space
        typename HBSplinesFESpace<TDim>::Pointer pNeighborFESpace = boost::dynamic_pointer_cast<HBSplinesFESpace<TDim> >(pNeighborPatch->pFESpace());
        if (pNeighborFESpace == NULL)
            KRATOS_THROW_ERROR(std::runtime_error, "The cast to HBSplinesFESpace is failed.", "")

        // get the correct basis function
        bf_t p_neighbor_bf;
        bool found = false;
        for(typename bf_container_t::iterator it = pNeighborFESpace->bf_begin(); it != pNeighborFESpace->bf_end(); ++it)
        {
            if ((*it)->EquationId() == equation_id)
            {
                p_neighbor_bf = *it;
                found = true;
            }
        }

        if (found)
        {
            if((echo_level & ECHO_REFINEMENT) == ECHO_REFINEMENT)
            {
                std::cout << "Neighbor patch " << pNeighborPatch->Id() << " of patch " << pPatch->Id() << " will be refined" << std::endl;
            }

            Refine(pNeighborPatch, p_neighbor_bf, refined_patches, echo_level);
        }
    }

    return std::make_pair(numbers, pnew_bfs);
}

template<>
inline void HBSplinesRefinementUtility_Helper<2>::RefineWindow(typename Patch<2>::Pointer pPatch,
        const std::vector<std::vector<double> >& window, const int& echo_level)
{
    if (pPatch->pFESpace()->Type() != HBSplinesFESpace<2>::StaticType())
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "only support the hierarchical B-Splines patch")

    // Type definitions
    typedef typename HBSplinesFESpace<2>::bf_t bf_t;
    typedef typename HBSplinesFESpace<2>::bf_container_t bf_container_t;
    typedef typename HBSplinesFESpace<2>::CellType CellType;
    typedef typename HBSplinesFESpace<2>::cell_t cell_t;
    typedef typename HBSplinesFESpace<2>::cell_container_t cell_container_t;
    typedef typename Patch<2>::ControlPointType ControlPointType;

    // extract the hierarchical B-Splines space
    typename HBSplinesFESpace<2>::Pointer pFESpace = boost::dynamic_pointer_cast<HBSplinesFESpace<2> >(pPatch->pFESpace());
    if (pFESpace == NULL)
        KRATOS_THROW_ERROR(std::runtime_error, "The cast to HBSplinesFESpace is failed.", "")

    // search and mark all basis functions need to refine on all level (starting from the last level) which support is contained in the refining domain
    for(typename bf_container_t::iterator it_bf = pFESpace->bf_begin(); it_bf != pFESpace->bf_end(); ++it_bf)
    {
        // get the bounding box (support domain of the basis function)
        std::vector<double> bounding_box = (*it_bf)->GetBoundingBox();

        // check if the bounding box lie in the refined domain
        // Remarks: this can be changed by a refinement indicator (i.e from error estimator)
        if(    bounding_box[0] >= window[0][0] && bounding_box[1] <= window[0][1]
            && bounding_box[2] >= window[1][0] && bounding_box[3] <= window[1][1] )
        {
            Refine(pPatch, (*it_bf)->Id(), echo_level);
            pPatch->pParentMultiPatch()->Enumerate();
        }
    }
}

template<>
inline void HBSplinesRefinementUtility_Helper<3>::RefineWindow(typename Patch<3>::Pointer pPatch,
        const std::vector<std::vector<double> >& window, const int& echo_level)
{
    if (pPatch->pFESpace()->Type() != HBSplinesFESpace<3>::StaticType())
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "only support the hierarchical B-Splines patch")

    // Type definitions
    typedef typename HBSplinesFESpace<3>::bf_t bf_t;
    typedef typename HBSplinesFESpace<3>::bf_container_t bf_container_t;
    typedef typename HBSplinesFESpace<3>::CellType CellType;
    typedef typename HBSplinesFESpace<3>::cell_t cell_t;
    typedef typename HBSplinesFESpace<3>::cell_container_t cell_container_t;
    typedef typename Patch<3>::ControlPointType ControlPointType;

    // extract the hierarchical B-Splines space
    typename HBSplinesFESpace<3>::Pointer pFESpace = boost::dynamic_pointer_cast<HBSplinesFESpace<3> >(pPatch->pFESpace());
    if (pFESpace == NULL)
        KRATOS_THROW_ERROR(std::runtime_error, "The cast to HBSplinesFESpace is failed.", "")

    // search and mark all basis functions need to refine on all level (starting from the last level) which support is contained in the refining domain
    for(typename bf_container_t::iterator it_bf = pFESpace->bf_begin(); it_bf != pFESpace->bf_end(); ++it_bf)
    {
        // get the bounding box (support domain of the basis function)
        std::vector<double> bounding_box = (*it_bf)->GetBoundingBox();

        // check if the bounding box lie in the refined domain
        // Remarks: this can be changed by a refinement indicator (i.e from error estimator)
        if(    bounding_box[0] >= window[0][0] && bounding_box[1] <= window[0][1]
            && bounding_box[2] >= window[1][0] && bounding_box[3] <= window[1][1]
            && bounding_box[4] >= window[2][0] && bounding_box[5] <= window[2][1] )
        {
            Refine(pPatch, (*it_bf)->Id(), echo_level);
            pPatch->pParentMultiPatch()->Enumerate();
        }
    }
}


template<int TDim>
inline void HBSplinesRefinementUtility_Helper<TDim>::LinearDependencyRefine(typename Patch<TDim>::Pointer pPatch, const std::size_t& refine_cycle, const int& echo_level)
{
    if (pPatch->pFESpace()->Type() != HBSplinesFESpace<TDim>::StaticType())
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "only support the hierarchical B-Splines patch")

    // Type definitions
    typedef typename HBSplinesFESpace<TDim>::bf_t bf_t;
    typedef typename HBSplinesFESpace<TDim>::bf_container_t bf_container_t;
    typedef typename HBSplinesFESpace<TDim>::CellType CellType;
    typedef typename HBSplinesFESpace<TDim>::cell_t cell_t;
    typedef typename HBSplinesFESpace<TDim>::cell_container_t cell_container_t;
    typedef typename HBSplinesFESpace<TDim>::domain_t domain_t;
    typedef typename Patch<TDim>::ControlPointType ControlPointType;

    // extract the hierarchical B-Splines space
    typename HBSplinesFESpace<TDim>::Pointer pFESpace = boost::dynamic_pointer_cast<HBSplinesFESpace<TDim> >(pPatch->pFESpace());
    if (pFESpace == NULL)
        KRATOS_THROW_ERROR(std::runtime_error, "The cast to HBSplinesFESpace is failed.", "")

    if(pFESpace->LastLevel() < 1) return;

    #ifdef ENABLE_PROFILING
    double start = OpenMPUtils::GetCurrentTime();
    #endif

    // rebuild support domain
    pFESpace->ClearSupportDomain();
    for(std::size_t level = 1; level <= pFESpace->LastLevel(); ++level)
    {
        domain_t p_domain = pFESpace->GetSupportDomain(level);

        // add the knots to the domain manager
        for(std::size_t next_level = level; next_level <= pFESpace->LastLevel(); ++next_level)
        {
            for(typename bf_container_t::iterator it_bf = pFESpace->bf_begin(); it_bf != pFESpace->bf_end(); ++it_bf)
            {
                if((*it_bf)->Level() == next_level)
                {
                    for(typename HBSplinesFESpace<TDim>::BasisFunctionType::cell_iterator it_cell = (*it_bf)->cell_begin(); it_cell != (*it_bf)->cell_end(); ++it_cell)
                    {
                        if(TDim == 2)
                        {
                            p_domain->AddXcoord((*it_cell)->LeftValue());
                            p_domain->AddXcoord((*it_cell)->RightValue());
                            p_domain->AddYcoord((*it_cell)->DownValue());
                            p_domain->AddYcoord((*it_cell)->UpValue());
                        }
                        else if(TDim == 3)
                        {
                            p_domain->AddXcoord((*it_cell)->LeftValue());
                            p_domain->AddXcoord((*it_cell)->RightValue());
                            p_domain->AddYcoord((*it_cell)->DownValue());
                            p_domain->AddYcoord((*it_cell)->UpValue());
                            p_domain->AddZcoord((*it_cell)->BelowValue());
                            p_domain->AddZcoord((*it_cell)->AboveValue());
                        }
                    }
                }
            }
        }

        // add the cells to the domain manager
        for(std::size_t next_level = level; next_level <= pFESpace->LastLevel(); ++next_level)
        {
            for(typename bf_container_t::iterator it_bf = pFESpace->bf_begin(); it_bf != pFESpace->bf_end(); ++it_bf)
            {
                if((*it_bf)->Level() == next_level)
                {
                    for(typename HBSplinesFESpace<TDim>::BasisFunctionType::cell_iterator it_cell = (*it_bf)->cell_begin(); it_cell != (*it_bf)->cell_end(); ++it_cell)
                    {
                        if(TDim == 2)
                        {
                            std::vector<double> box = {(*it_cell)->LeftValue(), (*it_cell)->RightValue(), (*it_cell)->DownValue(), (*it_cell)->UpValue()};
                            p_domain->AddCell(box);
                        }
                        else if(TDim == 3)
                        {
                            std::vector<double> box = {(*it_cell)->LeftValue(), (*it_cell)->RightValue(), (*it_cell)->DownValue(), (*it_cell)->UpValue(), (*it_cell)->BelowValue(), (*it_cell)->AboveValue()};
                            p_domain->AddCell(box);
                        }
                    }
                }
            }
        }

//            std::cout << "support domain level " << level << *p_domain << std::endl;
    }

    // refine based on the rule that if a bf has support domain contained in next level, it must be refined
    for(std::size_t level = 1; level <= pFESpace->LastLevel() - 1; ++level)
    {
        std::vector<std::size_t> refined_bfs;
        for(typename bf_container_t::iterator it_bf = pFESpace->bf_begin(); it_bf != pFESpace->bf_end(); ++it_bf)
        {
            // extract the support domain of the next level
            if((*it_bf)->Level() != level) continue;
            domain_t p_domain = pFESpace->GetSupportDomain(level + 1);

            // get the support domain of the bf
            std::vector<double> bounding_box = (*it_bf)->GetBoundingBox();

            // check if the bf support domain contained in the refined domain managed by the domain manager
            bool is_inside = p_domain->IsInside(bounding_box);

            if(is_inside)
                refined_bfs.push_back((*it_bf)->Id());
        }

        if(refined_bfs.size() > 0)
        {
            if((echo_level & ECHO_REFINEMENT) == ECHO_REFINEMENT)
            {
                std::cout << "Additional Bf";
                for(std::size_t i = 0; i < refined_bfs.size(); ++i)
                    std::cout << " " << refined_bfs[i];
                std::cout << " of level " << level << " will be refined to maintain the linear independence..." << std::endl;
            }

            for(std::size_t i = 0; i < refined_bfs.size(); ++i)
            {
                Refine(pPatch, refined_bfs[i], echo_level);
                pPatch->pParentMultiPatch()->Enumerate();
            }

            // perform another round to make sure all bfs has support domain in the domain manager of each level
            LinearDependencyRefine(pPatch, refine_cycle + 1, echo_level);
        }
    }

    #ifdef ENABLE_PROFILING
    std::cout << "LinearDependencyRefine cycle " << refine_cycle << " completed: " << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;
    #else
    std::cout << "LinearDependencyRefine cycle " << refine_cycle << " completed" << std::endl;
    #endif
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_REFINEMENT_UTILITY_H_INCLUDED defined

