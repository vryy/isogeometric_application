//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 24 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_FESPACE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_FESPACE_H_INCLUDED

// System includes
#include <vector>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "containers/pointer_vector_set.h"
#include "custom_utilities/bezier_utils.h"
#include "custom_utilities/bspline_utils.h"
#include "custom_utilities/fespace.h"
#include "custom_utilities/nurbs/knot_array_1d.h"
#include "custom_utilities/nurbs/bsplines_indexing_utility.h"
#include "custom_utilities/nurbs/cell_manager_1d.h"
#include "custom_utilities/nurbs/cell_manager_2d.h"
#include "custom_utilities/nurbs/cell_manager_3d.h"
#include "custom_utilities/nurbs/domain_manager.h"
#include "custom_utilities/nurbs/domain_manager_2d.h"
#include "custom_utilities/nurbs/domain_manager_3d.h"
#include "custom_utilities/hbsplines/hb_cell.h"
#include "custom_utilities/hbsplines/hbsplines_basis_function.h"

#define DEBUG_GEN_CELL
#define DEBUG_DESTROY

namespace Kratos
{

/**
This class represents the FESpace for a single hierarchical BSplines patch defined over parametric domain.
 */
template<int TDim>
class HBSplinesFESpace : public FESpace<TDim>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HBSplinesFESpace);

    /// Type definition
    typedef FESpace<TDim> BaseType;
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    typedef HBSplinesBasisFunction<TDim> BasisFunctionType;
    typedef typename BasisFunctionType::Pointer bf_t;
    struct bf_compare { bool operator() (const bf_t& lhs, const bf_t& rhs) const {return lhs->Id() < rhs->Id();} };
    typedef std::set<bf_t, bf_compare> bf_container_t;
    typedef typename bf_container_t::iterator bf_iterator;
    typedef typename bf_container_t::const_iterator bf_const_iterator;

    typedef HBCell<HBSplinesBasisFunction<TDim> > CellType;
    typedef CellManager<CellType> cell_container_t;
    typedef typename cell_container_t::cell_t cell_t;

    typedef DomainManager::Pointer domain_t;
    typedef std::map<std::size_t, domain_t> domain_container_t;

    typedef std::map<std::size_t, bf_t> function_map_t;

    /// Default constructor
    HBSplinesFESpace() : BaseType(), m_function_map_is_created(false), mLastLevel(1), mMaxLevel(10)
    {
        if (TDim == 1)
        {
            mpCellManager = typename cell_container_t::Pointer(new CellManager1D<CellType>());
        }
        else if (TDim == 2)
        {
            mpCellManager = typename cell_container_t::Pointer(new CellManager2D<CellType>());
        }
        else if(TDim == 3)
        {
            mpCellManager = typename cell_container_t::Pointer(new CellManager3D<CellType>());
        }
    }

    /// Destructor
    virtual ~HBSplinesFESpace()
    {
        #ifdef DEBUG_DESTROY
        std::cout << Type() << ", Addr = " << this << " is destroyed" << std::endl;
        #endif
    }

    /// Helper to create new HBSplinesFESpace pointer
    static typename HBSplinesFESpace<TDim>::Pointer Create()
    {
        return typename HBSplinesFESpace<TDim>::Pointer(new HBSplinesFESpace());
    }

    /// Add a already constructed basis function to the internal list
    void AddBf(bf_t p_bf)
    {
        mpBasisFuncs.insert(p_bf);
    }

    /// Check if the bf exists in the list; otherwise create new bf and return
    /// The Id is only used when the new bf is created. User must always check the Id of the returned function.
    bf_t CreateBf(const std::size_t& Id, const std::size_t& Level, const std::vector<std::vector<knot_t> >& rpKnots)
    {
        // search in the current list of basis functions, the one that has the same local knot vector with provided ones
        for(bf_iterator it = bf_begin(); it != bf_end(); ++it)
            if((*it)->Contain(rpKnots))
                return *it;

        // create the new bf and add the knot
        bf_t p_bf = bf_t(new BasisFunctionType(Id, Level));
        for (int dim = 0; dim < TDim; ++dim)
        {
            p_bf->SetLocalKnotVectors(dim, rpKnots[dim]);
            p_bf->SetInfo(dim, this->Order(dim));
        }
        mpBasisFuncs.insert(p_bf);
        m_function_map_is_created = false;

        return p_bf;
    }

    /// Remove the basis functions from the container
    void RemoveBf(bf_t p_bf)
    {
        mpBasisFuncs.erase(p_bf);
    }

    // Iterators for the basis functions
    bf_iterator bf_begin() {return mpBasisFuncs.begin();}
    bf_const_iterator bf_begin() const {return mpBasisFuncs.begin();}
    bf_iterator bf_end() {return mpBasisFuncs.end();}
    bf_const_iterator bf_end() const {return mpBasisFuncs.end();}

    /// Get the last id of the basis functions
    std::size_t LastId() const {bf_const_iterator it = bf_end(); --it; return (*it)->Id();}

    /// Get the last refinement level ain the hierarchical mesh
    const std::size_t& LastLevel() const {return mLastLevel;}

    /// Set the last level in the hierarchical mesh
    void SetLastLevel(const std::size_t& LastLevel) {mLastLevel = LastLevel;}

    /// Get the maximum level allowed in the hierarchical mesh
    const std::size_t& MaxLevel() const {return mMaxLevel;}

    /// Set the maximum level allowed in the hierarchical mesh
    void SetMaxLevel(const std::size_t& MaxLevel) {mMaxLevel = MaxLevel;}

    /// Set the order for the B-Splines
    void SetInfo(const std::size_t& i, const std::size_t& order) {mOrders[i] = order;}

    /// Get the order of the BSplines patch in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
    {
        if (i >= TDim) return 0;
        else return mOrders[i];
    }

    /// Get the number of basis functions defined over the BSplines
    virtual const std::size_t TotalNumber() const {return mpBasisFuncs.size();}

    /// Get the knot vector in i-direction, i=0..Dim
    /// User must be careful to use this function because it can modify the internal knot vectors
    knot_container_t& KnotVector(const std::size_t& i) {return mKnotVectors[i];}

    /// Get the knot vector in i-direction, i=0..Dim
    const knot_container_t& KnotVector(const std::size_t& i) const {return mKnotVectors[i];}

    /// Get the weights of all the basis functions
    std::vector<double> GetWeights() const
    {
        std::vector<double> weights(this->TotalNumber());
        std::size_t cnt = 0;
        for (bf_iterator it = bf_begin(); it != bf_end(); ++it)
            weights[cnt++] = (*it)->GetValue(CONTROL_POINT).W();
        return weights;
    }

    /// Get the string representing the type of the patch
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string representing the type of the patch
    static std::string StaticType()
    {
        std::stringstream ss;
        ss << "HBSplinesFESpace" << TDim << "D";
        return ss.str();
    }

    /// Validate the HBSplinesFESpace
    virtual bool Validate() const
    {
        // TODO validate more
        return BaseType::Validate();
    }

    /// Get the values of the basis function i at point xi
    /// REMARK: This function only returns the unweighted basis function value. To obtain the correct one, use WeightedFESpace
    virtual void GetValue(double& v, const std::size_t& i, const std::vector<double>& xi) const
    {
        std::size_t j = 0;
        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if (j == i)
            {
                v = (*it)->GetValue(xi);
                return;
            }
            ++j;
        }
        v = 0.0;
    }

    /// Get the values of the basis functions at point xi
    /// REMARK: This function only returns the unweighted basis function value. To obtain the correct one, use WeightedFESpace
    virtual void GetValue(std::vector<double>& values, const std::vector<double>& xi) const
    {
        if (values.size() != this->TotalNumber())
            values.resize(this->TotalNumber());
        std::size_t i = 0;
        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
            values[i++] = (*it)->GetValue(xi);
    }

    /// Get the derivative of the basis function i at point xi
    /// the output derivatives has the form of values[func_index][dim_index]
    /// REMARK: This function only returns the unweighted basis function derivatives. To obtain the correct one, use WeightedFESpace
    virtual void GetDerivative(std::vector<double>& values, const std::size_t& i, const std::vector<double>& xi) const
    {
        std::size_t j = 0;
        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if (j == i)
            {
                (*it)->GetDerivative(values, xi);
                return;
            }
            ++j;
        }
        if (values.size() != TDim)
            values.resize(TDim);
        for (int dim = 0; dim < TDim; ++dim) values[dim] = 0.0;
    }

    /// Get the derivative of the basis functions at point xi
    /// the output derivatives has the form of values[func_index][dim_index]
    /// REMARK: This function only returns the unweighted basis function derivatives. To obtain the correct one, use WeightedFESpace
    virtual void GetDerivative(std::vector<std::vector<double> >& values, const std::vector<double>& xi) const
    {
        if (values.size() != this->TotalNumber())
            values.resize(this->TotalNumber());
        std::size_t i = 0;
        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if (values[i].size() != TDim)
                values[i].resize(TDim);
            (*it)->GetDerivative(values[i], xi);
            ++i;
        }
    }

    /// Get the values and derivatives of the basis functions at point xi
    /// the output derivatives has the form of values[func_index][dim_index]
    /// REMARK: This function only returns the unweighted basis function derivatives. To obtain the correct one, use WeightedFESpace
    virtual void GetValueAndDerivative(std::vector<double>& values, std::vector<std::vector<double> >& derivatives, const std::vector<double>& xi) const
    {
        if (values.size() != this->TotalNumber())
            values.resize(this->TotalNumber());
        if (derivatives.size() != this->TotalNumber())
            derivatives.resize(this->TotalNumber());

        std::size_t i = 0;
        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            values[i] = (*it)->GetValue(xi);
            if (derivatives[i].size() != TDim)
                derivatives[i].resize(TDim);
            (*it)->GetDerivative(derivatives[i], xi);
            ++i;
        }
    }

    /// Compare between two BSplines patches in terms of parametric information
    virtual bool IsCompatible(const FESpace<TDim>& rOtherFESpace) const
    {
        if (rOtherFESpace.Type() != Type())
        {
            KRATOS_WATCH(rOtherFESpace.Type())
            KRATOS_WATCH(Type())
            std::cout << "WARNING!!! the other patch type is not " << Type() << std::endl;
            return false;
        }

        const HBSplinesFESpace<TDim>& rOtherHBSplinesFESpace = dynamic_cast<const HBSplinesFESpace<TDim>&>(rOtherFESpace);

        // compare the knot vectors and order information
        for (std::size_t i = 0; i < TDim; ++i)
        {
            if (!(this->Order(i)) == rOtherHBSplinesFESpace.Order(i))
                return false;
            if (!(this->KnotVector(i) == rOtherHBSplinesFESpace.KnotVector(i)))
                return false;
        }

        return true;
    }

    /// Get the refinement history
    const std::vector<std::size_t>& RefinementHistory() const {return mRefinementHistory;}

    /// Clear the refinement history
    void ClearRefinementHistory() {mRefinementHistory.clear();}

    /// Add the bf's id to the refinement history
    void RecordRefinementHistory(const std::size_t& Id) {mRefinementHistory.push_back(Id);}

    /// Clear the support domain container
    void ClearSupportDomain() {mSupportDomains.clear();}

    /// Get the domain manager for support domain for level k. In the case that does not exist, create the new one.
    domain_t GetSupportDomain(const std::size_t& Level)
    {
        domain_container_t::iterator it = mSupportDomains.find(Level);
        if(it != mSupportDomains.end())
            return it->second;
        else
        {
            domain_t p_domain;
            if(TDim == 2)
                p_domain = domain_t(new DomainManager2D(Level));
            else if(TDim == 3)
                p_domain = domain_t(new DomainManager3D(Level));
            mSupportDomains[Level] = p_domain;
            return p_domain;
        }
    }

    /// Reset the function indices to a given values.
    /// This is useful when assigning the id for the boundary patch.
    virtual void ResetFunctionIndices()
    {
        BaseType::mGlobalToLocal.clear();
        for (bf_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            (*it)->SetEquationId(-1);
        }
    }

    /// Reset the function indices to a given values.
    /// This is useful when assigning the id for the boundary patch.
    virtual void ResetFunctionIndices(const std::vector<std::size_t>& func_indices)
    {
        if (func_indices.size() != this->TotalNumber())
        {
            KRATOS_WATCH(this->TotalNumber())
            std::cout << "func_indices:";
            for (std::size_t i = 0; i < func_indices.size(); ++i)
                std::cout << " " << func_indices[i];
            std::cout << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, "The func_indices vector does not have the same size as total number of basis functions", "")
        }
        std::size_t cnt = 0;
        for (bf_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            (*it)->SetEquationId(func_indices[cnt]);
            BaseType::mGlobalToLocal[(*it)->EquationId()] = cnt;
            ++cnt;
        }
    }

    /// Enumerate the dofs of each grid function. The enumeration algorithm is pretty straightforward.
    /// If the dof does not have pre-existing value, which assume it is -1, it will be assigned the incremental value.
    virtual std::size_t& Enumerate(std::size_t& start)
    {
        BaseType::mGlobalToLocal.clear();
        std::size_t cnt = 0;
        for (bf_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if ((*it)->EquationId() == -1) (*it)->SetEquationId(start++);
            BaseType::mGlobalToLocal[(*it)->EquationId()] = cnt++;
        }

        return start;
    }

    /// Access the function indices (aka global ids)
    virtual std::vector<std::size_t> FunctionIndices() const
    {
        std::vector<std::size_t> func_indices(this->TotalNumber());
        std::size_t cnt = 0;
        for (bf_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            func_indices[cnt] = (*it)->EquationId();
            ++cnt;
        }
        return func_indices;
    }

    /// Update the function indices using a map. The map shall be the mapping from old index to new index.
    virtual void UpdateFunctionIndices(const std::map<std::size_t, std::size_t>& indices_map)
    {
        std::size_t cnt = 0;
        BaseType::mGlobalToLocal.clear();
        for (bf_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            std::map<std::size_t, std::size_t>::const_iterator it2 = indices_map.find((*it)->EquationId());

            if (it2 == indices_map.end())
            {
                std::cout << "WARNING!!! the indices_map does not contain " << (*it)->EquationId() << std::endl;
                continue;
            }

            (*it)->SetEquationId(it2->second);
            BaseType::mGlobalToLocal[(*it)->EquationId()] = cnt;
            ++cnt;
        }
    }

    /// Get the first equation_id in this space
    virtual std::size_t GetFirstEquationId() const
    {
        std::size_t first_id;

        for (bf_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if ((*it)->EquationId() == -1)
            {
                return -1;
            }
            else
            {
                if (it == bf_begin())
                    first_id = (*it)->EquationId();
                else
                    if ((*it)->EquationId() < first_id)
                        first_id = (*it)->EquationId();
            }
        }

        return first_id;
    }

    /// Get the last equation_id in this space
    virtual std::size_t GetLastEquationId() const
    {
        std::size_t last_id = -1;
        bool hit = false;

        for (bf_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if ((*it)->EquationId() != -1)
            {
                if (!hit)
                {
                    hit = true;
                    last_id = (*it)->EquationId();
                }
                else
                {
                    if ((*it)->EquationId() > last_id)
                        last_id = (*it)->EquationId();
                }
            }
        }

        return last_id;
    }

    /// Extract the index of the functions on the boundary
    virtual std::vector<std::size_t> ExtractBoundaryFunctionIndices(const BoundarySide& side) const
    {
        std::vector<std::size_t> func_indices;

        // firstly we organize the basis functions based on its equation_id
        std::map<std::size_t, bf_t> map_bfs;
        for (bf_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if ((*it)->IsOnSide(BOUNDARY_FLAG(side)))
                map_bfs[(*it)->EquationId()] = (*it);
        }

        // then we can extract the equation_id
        func_indices.resize(map_bfs.size());
        std::size_t cnt = 0;
        for (typename std::map<std::size_t, bf_t>::iterator it = map_bfs.begin(); it != map_bfs.end(); ++it)
        {
            func_indices[cnt++] = it->first;
        }

        return func_indices;
    }

    /// Assign the index for the functions on the boundary
    virtual void AssignBoundaryFunctionIndices(const BoundarySide& side, const std::vector<std::size_t>& func_indices)
    {
        // firstly we organize the basis functions based on its equation_id
        std::map<std::size_t, bf_t> map_bfs;
        for (bf_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if ((*it)->IsOnSide(BOUNDARY_FLAG(side)))
                map_bfs[(*it)->EquationId()] = (*it);
        }

        // then we can assign the equation_id incrementally
        std::size_t cnt = 0;
        for (typename std::map<std::size_t, bf_t>::iterator it = map_bfs.begin(); it != map_bfs.end(); ++it)
        {
            it->second->SetEquationId(func_indices[cnt++]);
        }
    }

    /// Construct the boundary FESpace based on side
    virtual typename FESpace<TDim-1>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    {
        typedef HBSplinesFESpace<TDim-1> BoundaryFESpaceType;
        typename BoundaryFESpaceType::Pointer pBFESpace = typename BoundaryFESpaceType::Pointer(new BoundaryFESpaceType());

        for (bf_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if ((*it)->IsOnSide(BOUNDARY_FLAG(side)))
            {
                typename BoundaryFESpaceType::bf_t pNewSubBf;

                if (TDim == 2)
                {
                    if ((side == _BLEFT_) || (side == _BRIGHT_))
                        pNewSubBf = (*it)->Project(1);
                    if ((side == _BTOP_) || (side == _BBOTTOM_))
                        pNewSubBf = (*it)->Project(0);
                }
                else if (TDim == 3)
                {
                    if ((side == _BFRONT_) || (side == _BBACK_))
                        pNewSubBf = (*it)->Project(0);
                    if ((side == _BLEFT_) || (side == _BRIGHT_))
                        pNewSubBf = (*it)->Project(1);
                    if ((side == _BTOP_) || (side == _BBOTTOM_))
                        pNewSubBf = (*it)->Project(2);
                }

                pBFESpace->AddBf(pNewSubBf);
            }
        }

        if (TDim == 2)
        {
            if ((side == _BLEFT_) || (side == _BRIGHT_))
            {
                pBFESpace->SetInfo(0, this->Order(1));
                pBFESpace->KnotVector(0) = this->KnotVector(1);
            }
            if ((side == _BTOP_) || (side == _BBOTTOM_))
            {
                pBFESpace->SetInfo(0, this->Order(0));
                pBFESpace->KnotVector(0) = this->KnotVector(0);
            }
        }
        else if (TDim == 3)
        {
            if ((side == _BFRONT_) || (side == _BBACK_))
            {
                pBFESpace->SetInfo(0, this->Order(1));
                pBFESpace->SetInfo(1, this->Order(2));
                pBFESpace->KnotVector(0) = this->KnotVector(1);
                pBFESpace->KnotVector(1) = this->KnotVector(2);
            }
            if ((side == _BLEFT_) || (side == _BRIGHT_))
            {
                pBFESpace->SetInfo(0, this->Order(2));
                pBFESpace->SetInfo(1, this->Order(0));
                pBFESpace->KnotVector(0) = this->KnotVector(2);
                pBFESpace->KnotVector(1) = this->KnotVector(0);
            }
            if ((side == _BTOP_) || (side == _BBOTTOM_))
            {
                pBFESpace->SetInfo(0, this->Order(0));
                pBFESpace->SetInfo(1, this->Order(1));
                pBFESpace->KnotVector(0) = this->KnotVector(0);
                pBFESpace->KnotVector(1) = this->KnotVector(1);
            }
        }

        pBFESpace->SetLastLevel(this->LastLevel());

        // construct the cells from the boundary basis functions
        typename BoundaryFESpaceType::cell_container_t::Pointer pnew_cells;
        double cell_tol = pBFESpace->pCellManager()->GetTolerance();

        if (TDim == 2)
        {
            pnew_cells = typename BoundaryFESpaceType::cell_container_t::Pointer(new CellManager1D<typename BoundaryFESpaceType::CellType>());

            for (typename BoundaryFESpaceType::bf_iterator it = pBFESpace->bf_begin(); it != pBFESpace->bf_end(); ++it)
            {
                for(std::size_t i1 = 0; i1 < pBFESpace->Order(0) + 1; ++i1)
                {
                    knot_t pLeft = (*it)->LocalKnots(0)[i1];
                    knot_t pRight = (*it)->LocalKnots(0)[i1 + 1];

                    // check if the cell domain length is nonzero
                    double length = (pRight->Value() - pLeft->Value());
                    if(fabs(length) > cell_tol)
                    {
                        std::vector<knot_t> pKnots = {pLeft, pRight};
                        typename BoundaryFESpaceType::cell_t pnew_cell = pBFESpace->pCellManager()->CreateCell(pKnots);
                        pnew_cell->SetLevel(this->LastLevel());
                        (*it)->AddCell(pnew_cell);
                        pnew_cell->AddBf(*it);
                        pnew_cells->insert(pnew_cell);
                    }
                }
            }
        }
        else if (TDim == 3)
        {
            pnew_cells = typename BoundaryFESpaceType::cell_container_t::Pointer(new CellManager2D<typename BoundaryFESpaceType::CellType>());

            for (typename BoundaryFESpaceType::bf_iterator it = pBFESpace->bf_begin(); it != pBFESpace->bf_end(); ++it)
            {
                for(std::size_t i1 = 0; i1 < pBFESpace->Order(0) + 1; ++i1)
                {
                    knot_t pLeft = (*it)->LocalKnots(0)[i1];
                    knot_t pRight = (*it)->LocalKnots(0)[i1 + 1];

                    for(std::size_t j1 = 0; j1 < pBFESpace->Order(1) + 1; ++j1)
                    {
                        knot_t pDown = (*it)->LocalKnots(1)[j1];
                        knot_t pUp = (*it)->LocalKnots(1)[j1 + 1];

                        // check if the cell domain area is nonzero
                        double area = (pRight->Value() - pLeft->Value()) * (pUp->Value() - pDown->Value());
                        if(sqrt(fabs(area)) > cell_tol)
                        {
                            std::vector<knot_t> pKnots = {pLeft, pRight, pDown, pUp};
                            typename BoundaryFESpaceType::cell_t pnew_cell = pBFESpace->pCellManager()->CreateCell(pKnots);
                            pnew_cell->SetLevel(this->LastLevel());
                            (*it)->AddCell(pnew_cell);
                            pnew_cell->AddBf(*it);
                            pnew_cells->insert(pnew_cell);
                        }
                    }
                }
            }
        }

        // collapse the overlapping cells
        pBFESpace->pCellManager()->CollapseCells();

        return pBFESpace;
    }

    /// Construct the boundary FESpace based on side and rotation
    virtual typename FESpace<TDim-1>::Pointer ConstructBoundaryFESpace(const BoundarySide& side,
        const std::map<std::size_t, std::size_t>& local_parameter_map, const std::vector<BoundaryDirection>& directions) const
    {
        return this->ConstructBoundaryFESpace(side);
    }

    /// Get the underlying cell manager
    typename cell_container_t::Pointer pCellManager() {return mpCellManager;}

    /// Get the underlying cell manager
    typename cell_container_t::ConstPointer pCellManager() const {return mpCellManager;}

    /// Update the basis functions for all cells. This function must be called before any operation on cell is required.
    void UpdateCells()
    {
        // for each cell compute the extraction operator and add to the anchor
        Vector Crow;
        for(typename cell_container_t::iterator it_cell = mpCellManager->begin(); it_cell != mpCellManager->end(); ++it_cell)
        {
            (*it_cell)->Reset();
            for(typename CellType::bf_iterator it_bf = (*it_cell)->bf_begin(); it_bf != (*it_cell)->bf_end(); ++it_bf)
            {
                (*it_bf)->ComputeExtractionOperator(Crow, *it_cell);
                (*it_cell)->AddAnchor((*it_bf)->Id(), (*it_bf)->GetValue(CONTROL_POINT).W(), Crow);
            }
        }
    }

    /// Create the cell manager for all the cells in the support domain of the HBSplinesFESpace
    virtual typename BaseType::cell_container_t::Pointer ConstructCellManager() const
    {
        // create the compatible cell manager and add to the list
        typename BaseType::cell_container_t::Pointer pCompatCellManager;

        if (TDim == 1)
            pCompatCellManager = CellManager1D<typename BaseType::cell_container_t::CellType>::Create();
        else if (TDim == 2)
            pCompatCellManager = CellManager2D<typename BaseType::cell_container_t::CellType>::Create();
        else if (TDim == 3)
            pCompatCellManager = CellManager3D<typename BaseType::cell_container_t::CellType>::Create();

        for(typename cell_container_t::iterator it_cell = mpCellManager->begin(); it_cell != mpCellManager->end(); ++it_cell)
            pCompatCellManager->insert(*it_cell);

        return pCompatCellManager;
    }

    /// Get the basis functions on side
    std::vector<bf_t> GetBoundaryBfs(const std::size_t& boundary_id) const
    {
        // firstly we organize the basis functions based on its equation_id
        std::map<std::size_t, bf_t> map_bfs;
        for (bf_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if ((*it)->IsOnSide(boundary_id))
                map_bfs[(*it)->EquationId()] = (*it);
        }

        // then we can extract the equation_id
        std::vector<bf_t> bf_list(map_bfs.size());
        std::size_t cnt = 0;
        for (typename std::map<std::size_t, bf_t>::iterator it = map_bfs.begin(); it != map_bfs.end(); ++it)
        {
            bf_list[cnt++] = it->second;
        }

        return bf_list;
    }

    /// Overload operator[], this allows to access the basis function randomly based on index
    bf_t operator[](const std::size_t& i)
    {
        typename bf_container_t::iterator it = mpBasisFuncs.begin();
        std::advance(it, i);
        return *it;
    }

    /// Overload operator(), this allows to access the basis function based on its id
    bf_t operator()(const std::size_t& Id)
    {
        // create the index map if it's not created yet
        if(!m_function_map_is_created)
            CreateFunctionsMap();

        // return the bf if its Id exist in the list
        typename function_map_t::iterator it = mFunctionsMap.find(Id);
        if(it != mFunctionsMap.end())
            return it->second;
        else
            KRATOS_THROW_ERROR(std::runtime_error, "Access index is not found:", Id)
    }

    /// Overload assignment operator
    HBSplinesFESpace<TDim>& operator=(const HBSplinesFESpace<TDim>& rOther)
    {
        // TODO copy more
        KRATOS_THROW_ERROR(std::logic_error, "The assignment oprator is not complete", "")
        BaseType::operator=(rOther);
        return *this;
    }

    /// Clone this FESpace, this is a deep copy operation
    virtual typename FESpace<TDim>::Pointer Clone() const
    {
        typename HBSplinesFESpace<TDim>::Pointer pNewFESpace = typename HBSplinesFESpace<TDim>::Pointer(new HBSplinesFESpace<TDim>());
        *pNewFESpace = *this;
        return pNewFESpace;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Type() << ", Addr = " << this << ", n = " << this->TotalNumber();
        rOStream << ", p = (";
        for (std::size_t dim = 0; dim < TDim; ++dim)
            rOStream << " " << this->Order(dim);
        rOStream << "), number of levels = " << mLastLevel << std::endl;

        rOStream << "###############Begin knot vectors################" << std::endl;
        for (int dim = 0; dim < TDim; ++dim)
        {
            rOStream << "knot vector " << dim+1 << ":";
            for (std::size_t i = 0; i < this->KnotVector(dim).size(); ++i)
                rOStream << " " << this->KnotVector(dim)[i];
            rOStream << std::endl;
        }
        rOStream << "###############End knot vectors##################" << std::endl;
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        // print the basis functions
        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            rOStream << *(*it) << std::endl;
        }

        // print the cells in each level
        for (std::size_t level = 1; level < mLastLevel+1; ++level)
        {
            rOStream << "###############Begin cells at level " << level << "################" << std::endl;
            std::size_t n = 0;
            for(typename cell_container_t::iterator it = mpCellManager->begin(); it != mpCellManager->end(); ++it)
                if((*it)->Level() == level)
                    rOStream << "(" << ++n << ") " << *(*it) << std::endl;
            rOStream << "###############End cells at level " << level << "################" << std::endl;
        }
    }

private:

    unsigned int mEchoLevel;

    boost::array<std::size_t, TDim> mOrders;

    boost::array<knot_container_t, TDim> mKnotVectors;

    std::size_t mLastLevel;
    std::size_t mMaxLevel;

    typename cell_container_t::Pointer mpCellManager;

    bf_container_t mpBasisFuncs;
    mutable function_map_t mFunctionsMap; // map from basis function id to the basis function. It's mainly used to search for the bf quickly. But it needs to be re-initialized whenever new bf is added to the set
    bool m_function_map_is_created;

    void CreateFunctionsMap()
    {
        mFunctionsMap.clear();
        for(bf_iterator it = bf_begin(); it != bf_end(); ++it)
            mFunctionsMap[(*it)->Id()] = *it;
        m_function_map_is_created = true;
    }

    domain_container_t mSupportDomains; // this domain manager manages the support of all bfs in each level

    std::vector<std::size_t> mRefinementHistory;
};

/**
 * Template specific instantiation for null-D BSplines patch to terminate the compilation
 */
template<>
class HBSplinesFESpace<0> : public FESpace<0>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HBSplinesFESpace);

    /// Type definition
    typedef FESpace<0> BaseType;
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;
    typedef HBSplinesBasisFunction<0> BasisFunctionType;
    typedef typename BasisFunctionType::Pointer bf_t;
    struct bf_compare { bool operator() (const bf_t& lhs, const bf_t& rhs) const {return lhs->Id() < rhs->Id();} };
    typedef std::set<bf_t, bf_compare> bf_container_t;
    typedef typename bf_container_t::iterator bf_iterator;
    typedef typename bf_container_t::const_iterator bf_const_iterator;

    typedef HBCell<HBSplinesBasisFunction<0> > CellType;
    typedef CellManager<CellType> cell_container_t;
    typedef typename cell_container_t::cell_t cell_t;

    /// Default constructor
    HBSplinesFESpace() : BaseType() {}

    /// Destructor
    virtual ~HBSplinesFESpace() {}

    /// Get the order of the BSplines patch in specific direction
    virtual const std::size_t Order(const std::size_t& i) const {return 0;}

    /// Get the number of basis functions defined over the BSplines HBSplinesFESpace
    virtual const std::size_t Number() const {return 0;}

    /// Get the string describing the type of the patch
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the patch
    static std::string StaticType()
    {
        return "HBSplinesFESpace0D";
    }

    /// Validate the HBSplinesFESpace before using
    virtual bool Validate() const
    {
        return BaseType::Validate();
    }

    /// Compare between two BSplines patches in terms of parametric information
    virtual bool IsCompatible(const FESpace<0>& rOtherFESpace) const
    {
        if (rOtherFESpace.Type() != Type())
        {
            KRATOS_WATCH(rOtherFESpace.Type())
            KRATOS_WATCH(Type())
            std::cout << "WARNING!!! the other FESpace type is not " << Type() << std::endl;
            return false;
        }

        return true;
    }

    /// Add a already constructed basis function to the internal list
    void AddBf(bf_t p_bf) {}

    // Iterators for the basis functions
    bf_iterator bf_begin() {return mpBasisFuncs.begin();}
    bf_const_iterator bf_begin() const {return mpBasisFuncs.begin();}
    bf_iterator bf_end() {return mpBasisFuncs.end();}
    bf_const_iterator bf_end() const {return mpBasisFuncs.end();}

    /// Set the last level in the hierarchical mesh
    void SetLastLevel(const std::size_t& LastLevel) {}

    /// Set the order for the B-Splines
    void SetInfo(const std::size_t& i, const std::size_t& order) {}

    /// Get the knot vector in i-direction, i=0..Dim
    /// User must be careful to use this function because it can modify the internal knot vectors
    knot_container_t& KnotVector(const std::size_t& i) {return mDummyKnotVector;}

    /// Get the knot vector in i-direction, i=0..Dim
    const knot_container_t& KnotVector(const std::size_t& i) const {return mDummyKnotVector;}

    /// Get the underlying cell manager
    typename cell_container_t::Pointer pCellManager() {return mpCellManager;}

    /// Get the underlying cell manager
    typename cell_container_t::ConstPointer pCellManager() const {return mpCellManager;}

    /// Update the basis functions for all cells. This function must be called before any operation on cell is required.
    void UpdateCells() {}

private:

    knot_container_t mDummyKnotVector;
    bf_container_t mpBasisFuncs;
    typename cell_container_t::Pointer mpCellManager;
};

/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const HBSplinesFESpace<TDim>& rThis)
{
    rOStream << "-------------Begin HBSplinesFESpace Info-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "-------------End HBSplinesFESpace Info-------------" << std::endl;
    return rOStream;
}

} // namespace Kratos.

#undef DEBUG_GEN_CELL    /// Get the underlying cell manager
#undef DEBUG_DESTROY

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_FESPACE_H_INCLUDED defined
