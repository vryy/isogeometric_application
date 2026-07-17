//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 24 Nov 2017 $
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
#include "custom_utilities/bezier_utils.h"
#include "custom_utilities/bspline_utils.h"
#include "custom_utilities/nurbs/bcell_manager.h"
#include "custom_utilities/nurbs/pbbsplines_fespace.h"
#include "custom_utilities/nurbs/domain_manager.h"
#include "custom_utilities/nurbs/domain_manager_2d.h"
#include "custom_utilities/nurbs/domain_manager_3d.h"
#include "custom_utilities/hbsplines/hb_cell.h"
#include "custom_utilities/hbsplines/hbsplines_basis_function.h"

#define DEBUG_GEN_CELL

namespace Kratos
{

/**
This class represents the FESpace for a single hierarchical BSplines patch defined over parametric domain.
 */
template<int TDim, typename TLocalCoordinateType = double>
class HBSplinesFESpace : public PBBSplinesFESpace<TDim, TLocalCoordinateType, HBSplinesBasisFunction<TDim>, BCellManager<TDim, typename HBSplinesBasisFunction<TDim>::CellType> >
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HBSplinesFESpace);

    /// Type definition
    typedef PBBSplinesFESpace<TDim, TLocalCoordinateType, HBSplinesBasisFunction<TDim>, BCellManager<TDim, typename HBSplinesBasisFunction<TDim>::CellType> > BaseType;
    typedef HBSplinesFESpace<TDim, TLocalCoordinateType> ThisType;
    typedef typename BaseType::FESpaceType FESpaceType;
    typedef typename BaseType::knot_container_t knot_container_t;
    typedef typename BaseType::knot_t knot_t;

    typedef typename BaseType::BasisFunctionType BasisFunctionType;
    typedef typename BaseType::bf_t bf_t;
    typedef typename BaseType::bf_container_t bf_container_t;
    typedef typename BaseType::bf_iterator bf_iterator;
    typedef typename BaseType::bf_const_iterator bf_const_iterator;

    typedef typename BaseType::CellType CellType;
    typedef typename BaseType::cell_container_t cell_container_t;
    typedef typename BaseType::cell_t cell_t;

    typedef DomainManager::Pointer domain_t;
    typedef std::map<std::size_t, domain_t> domain_container_t;

    struct bf_compare
    {
        constexpr bool operator()( const bf_t& lhs, const bf_t& rhs ) const
        {
            return *lhs < *rhs;
        }
    };

    /// Default constructor
    HBSplinesFESpace() : BaseType(), mLastLevel(1), mMaxLevel(99)
    {
        mKnotVectors.resize(mLastLevel + 1);
    }

    /// Destructor
    ~HBSplinesFESpace() override
    {
#ifdef ISOGEOMETRIC_DEBUG_DESTROY
        std::cout << Type() << ", Addr = " << this << " is destroyed" << std::endl;
#endif
    }

    /// Helper function to cast to HBSplinesFESpace pointer
    static typename ThisType::Pointer Cast(typename FESpaceType::Pointer pFESpace)
    {
#ifdef SD_APP_FORWARD_COMPATIBILITY
        return std::dynamic_pointer_cast<ThisType>(pFESpace);
#else
        return boost::dynamic_pointer_cast<ThisType>(pFESpace);
#endif
    }

    /// Helper function to cast to HBSplinesFESpace pointer
    static typename ThisType::ConstPointer Cast(typename FESpaceType::ConstPointer pFESpace)
    {
#ifdef SD_APP_FORWARD_COMPATIBILITY
        return std::dynamic_pointer_cast<const ThisType>(pFESpace);
#else
        return boost::dynamic_pointer_cast<const ThisType>(pFESpace);
#endif
    }

    /// Helper to create new HBSplinesFESpace pointer
    static typename ThisType::Pointer Create()
    {
        return typename ThisType::Pointer(new ThisType());
    }

    /// Check if the bf exists in the list; otherwise create new bf and return
    /// The Id is only used when the new bf is created. User must always check the Id of the returned function.
    bf_t CreateBf(std::size_t Id, std::size_t Level, const std::vector<std::vector<knot_t> >& rpKnots)
    {
        double Tol = this->KnotVector(0).GetResolution();

        // search in the current list of basis functions, the one that has the same local knot vector with provided ones
        for (bf_iterator it = BaseType::bf_begin(); it != BaseType::bf_end(); ++it)
        {
            if (it->CompareLocalKnots(rpKnots, Tol))
            {
                return *it.base();
            }
        }

        // create the new bf and add the knot
        bf_t p_bf = bf_t(new BasisFunctionType(Id, Level));
        for (int dim = 0; dim < TDim; ++dim)
        {
            p_bf->SetLocalKnotVectors(dim, rpKnots[dim]);
            p_bf->SetInfo(dim, this->Order(dim));
        }
        BaseType::mpBasisFuncs.insert(BaseType::mpBasisFuncs.end(), p_bf);

        return p_bf;
    }

    /// Update the basis functions for all cells. This function must be called before any operation on cell is required.
    void UpdateCells() override
    {
        this->ResetCells();

        // for each cell compute the extraction operator and add to the anchor
        for (auto it_cell = BaseType::mpCellManager->begin(); it_cell != BaseType::mpCellManager->end(); ++it_cell)
        {
            for (auto it_bf = (*it_cell)->bf_begin(); it_bf != (*it_cell)->bf_end(); ++it_bf)
            {
                const BasisFunctionType& bf = *(it_bf->lock());
                if (!bf.Is(ACTIVE)) continue;
                Vector Crow;
                bf.ComputeExtractionOperator(Crow, *it_cell);
                (*it_cell)->AddAnchor(bf.EquationId(), bf.GetData(CONTROL_POINT).W(), Crow);
            }
        }
    }

    /// Get the knot vector in i-direction, i=0..Dim, for the first level
    /// User must be careful to use this function because it can modify the internal knot vectors
    void SetKnotVector(std::size_t idir, const knot_container_t& p_knot_vector)
    {
        if (mKnotVectors.size() < mLastLevel+1)
            mKnotVectors.resize(mLastLevel+1);
        mKnotVectors[0][idir] = p_knot_vector;

        // construct knot vectors for the remaining levels
        for (std::size_t j = 0; j < mLastLevel; ++j)
            mKnotVectors[j+1][idir] = mKnotVectors[j][idir].CloneAndRefineInTheMiddle();
    }

    /// Get the knot vector in i-direction, i=0..Dim, of the first level
    const knot_container_t& KnotVector(std::size_t i) const {return mKnotVectors[0][i];}

    /// Get the knot vector in i-direction, i=0..Dim, of a specific level
    const knot_container_t& KnotVector(std::size_t lvl, std::size_t i) const {return mKnotVectors[lvl-1][i];}

    /// Get the last refinement level ain the hierarchical mesh
    std::size_t LastLevel() const {return mLastLevel;}

    /// Set the last level in the hierarchical mesh
    void SetLastLevel(std::size_t LastLevel)
    {
        mLastLevel = LastLevel;

        if (mLastLevel + 1 > mKnotVectors.size())
        {
            std::size_t current_size = mKnotVectors.size();
            mKnotVectors.resize(mLastLevel + 1);

            std::size_t n_add_level = (mLastLevel + 1) - current_size;

            for (std::size_t i = 0; i < TDim; ++i)
            {
                for (std::size_t j = current_size; j < current_size + n_add_level; ++j)
                    mKnotVectors[j][i] = mKnotVectors[j-1][i].CloneAndRefineInTheMiddle();
            }
        }
        else
        {
            mKnotVectors.resize(mLastLevel + 1);
        }
    }

    /// Get the maximum level allowed in the hierarchical mesh
    std::size_t MaxLevel() const {return mMaxLevel;}

    /// Set the maximum level allowed in the hierarchical mesh
    void SetMaxLevel(std::size_t MaxLevel) {mMaxLevel = MaxLevel;}

    /// Get the string representing the type of the patch
    std::string Type() const override
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
    bool Validate() const override
    {
        // TODO validate more
        return BaseType::Validate();
    }

    /// Compare between two BSplines patches in terms of parametric information
    bool IsCompatible(const FESpaceType& rOtherFESpace) const override
    {
        if (rOtherFESpace.Type() != Type())
        {
            KRATOS_WATCH(rOtherFESpace.Type())
            KRATOS_WATCH(Type())
            std::cout << "WARNING!!! the other patch type is not " << Type() << std::endl;
            return false;
        }

        const ThisType* pOtherHBSplinesFESpace = dynamic_cast<const ThisType*>(&rOtherFESpace);
        if (pOtherHBSplinesFESpace == nullptr)
            return false;

        // compare the knot vectors and order information
        for (std::size_t i = 0; i < TDim; ++i)
        {
            if (!(this->Order(i)) == pOtherHBSplinesFESpace->Order(i))
            {
                return false;
            }
            if (!(this->KnotVector(i) == pOtherHBSplinesFESpace->KnotVector(i)))
            {
                return false;
            }
        }

        return true;
    }

    /// Get the basis functions based on boundary flag. This allows to extract the corner bf.
    std::set<bf_t, bf_compare> ExtractBoundaryBfsByFlag(std::size_t boundary_flag) const
    {
        std::set<bf_t, bf_compare> bf_set;

        for (bf_const_iterator it_bf = this->bf_begin(); it_bf != this->bf_end(); ++it_bf)
        {
            if (it_bf->IsOnSide(boundary_flag) && it_bf->Is(ACTIVE))
            {
                bf_set.insert(*it_bf.base());
            }
        }

        return bf_set;
    }

    /// Extract the index of the functions on the boundaries
    std::vector<std::size_t> ExtractBoundaryFunctionIndicesByFlag(int boundary_flag) const override
    {
        std::set<bf_t, bf_compare> bfs = this->ExtractBoundaryBfsByFlag(boundary_flag);

        std::vector<std::size_t> func_indices;
        for (auto it = bfs.begin(); it != bfs.end(); ++it)
        {
            func_indices.push_back((*it)->EquationId());
        }

        return func_indices;
    }

    /// Extract the index of the functions on the boundary
    std::vector<std::size_t> ExtractBoundaryFunctionIndices(const BoundarySide side) const override
    {
        // collect the basis and put into an organized set
        std::set<bf_t, bf_compare> set_bfs;
        for (bf_const_iterator it = this->bf_begin(); it != this->bf_end(); ++it)
        {
            if (it->IsOnSide(BOUNDARY_FLAG(side)) && it->Is(ACTIVE))
            {
                set_bfs.insert(*it.base());
            }
        }

        // then we can extract the equation_id
        std::vector<std::size_t> func_indices;
        for (auto it = set_bfs.begin(); it != set_bfs.end(); ++it)
        {
            func_indices.push_back((*it)->EquationId());
        }

        return func_indices;
    }

    /// Assign the index for the functions on the boundary
    void AssignBoundaryFunctionIndices(const BoundarySide side, const std::vector<std::size_t>& func_indices, const bool override) override
    {
        // collect the basis and put into an organized set
        std::set<bf_t, bf_compare> set_bfs;
        for (bf_const_iterator it = this->bf_begin(); it != this->bf_end(); ++it)
        {
            if (it->IsOnSide(BOUNDARY_FLAG(side)) && it->Is(ACTIVE))
            {
                set_bfs.insert(*it.base());
            }
        }

        // then we can assign the equation_id incrementally
        std::size_t cnt = 0;
        for (auto it = set_bfs.begin(); it != set_bfs.end(); ++it, ++cnt)
        {
            if (func_indices[cnt] != -1)
            {
                if (override)
                {
                    const auto local_id = this->LocalId((*it)->EquationId()); // record the local id
                    (*it)->SetEquationId(func_indices[cnt]);
                    BaseType::mGlobalToLocal[(*it)->EquationId()] = local_id; // reassign the local id
                }
                else
                {
                    if ((*it)->EquationId() == -1)
                        (*it)->SetEquationId(func_indices[cnt]);
                }
            }
        }
    }

    /// Get the refinement history
    const std::vector<std::size_t>& RefinementHistory() const {return mRefinementHistory;}

    /// Clear the refinement history
    void ClearRefinementHistory() {mRefinementHistory.clear();}

    /// Add the bf's id to the refinement history
    void RecordRefinementHistory(std::size_t Id) {mRefinementHistory.push_back(Id);}

    /// Clear the support domain container
    void ClearSupportDomain() {mSupportDomains.clear();}

    /// Get the domain manager for support domain for level k. In the case that does not exist, create the new one.
    domain_t GetSupportDomain(std::size_t Level)
    {
        domain_container_t::iterator it = mSupportDomains.find(Level);
        if (it != mSupportDomains.end())
        {
            return it->second;
        }
        else
        {
            domain_t p_domain;
            if constexpr (TDim == 2)
            {
                p_domain = domain_t(new DomainManager2D(Level));
            }
            else if constexpr (TDim == 3)
            {
                p_domain = domain_t(new DomainManager3D(Level));
            }
            mSupportDomains[Level] = p_domain;
            for (int i = 0; i < TDim; ++i)
                p_domain->SetTolerance(i, this->KnotVector(i).GetResolution());
            return p_domain;
        }
    }

    /// Build the support domain
    void BuildSupportDomain()
    {
        for (std::size_t level = 1; level <= this->LastLevel(); ++level)
        {
            domain_t p_domain = this->GetSupportDomain(level);

            // add the knots to the domain manager
            for (std::size_t next_level = level; next_level <= this->LastLevel(); ++next_level)
            {
                for (auto it_bf = this->bf_begin(); it_bf != this->bf_end(); ++it_bf)
                {
                    if (it_bf->Level() == next_level && it_bf->Is(ACTIVE))
                    {
                        for (auto it_cell = it_bf->cell_begin(); it_cell != it_bf->cell_end(); ++it_cell)
                        {
                            if constexpr (TDim > 0)
                            {
                                p_domain->AddXcoord((*it_cell)->XiMinValue());
                                p_domain->AddXcoord((*it_cell)->XiMaxValue());
                            }

                            if constexpr (TDim > 1)
                            {
                                p_domain->AddYcoord((*it_cell)->EtaMinValue());
                                p_domain->AddYcoord((*it_cell)->EtaMaxValue());
                            }

                            if constexpr (TDim > 2)
                            {
                                p_domain->AddZcoord((*it_cell)->ZetaMinValue());
                                p_domain->AddZcoord((*it_cell)->ZetaMaxValue());
                            }
                        }
                    }
                }
            }

            // add the cells to the domain manager
            for (std::size_t next_level = level; next_level <= this->LastLevel(); ++next_level)
            {
                for (auto it_bf = this->bf_begin(); it_bf != this->bf_end(); ++it_bf)
                {
                    if (it_bf->Level() == next_level && it_bf->Is(ACTIVE))
                    {
                        for (auto it_cell = it_bf->cell_begin(); it_cell != it_bf->cell_end(); ++it_cell)
                        {
                            if constexpr (TDim == 1)
                            {
                                std::vector<TLocalCoordinateType> box = {(*it_cell)->XiMinValue(), (*it_cell)->XiMaxValue()};
                                p_domain->AddCell(box);
                            }
                            else if constexpr (TDim == 2)
                            {
                                std::vector<TLocalCoordinateType> box = {(*it_cell)->XiMinValue(), (*it_cell)->XiMaxValue(), (*it_cell)->EtaMinValue(), (*it_cell)->EtaMaxValue()};
                                p_domain->AddCell(box);
                            }
                            else if constexpr (TDim == 3)
                            {
                                std::vector<TLocalCoordinateType> box = {(*it_cell)->XiMinValue(), (*it_cell)->XiMaxValue(), (*it_cell)->EtaMinValue(), (*it_cell)->EtaMaxValue(), (*it_cell)->ZetaMinValue(), (*it_cell)->ZetaMaxValue()};
                                p_domain->AddCell(box);
                            }
                        }
                    }
                }
            }
    //            std::cout << "support domain level " << level << *p_domain << std::endl;
        }
    }

    /// Rebuild the support domain
    void RebuildSupportDomain()
    {
        this->ClearSupportDomain();
        this->BuildSupportDomain();
    }

    /// Construct the boundary FESpace based on side
    typename FESpace<TDim-1, TLocalCoordinateType>::Pointer ConstructBoundaryFESpace(const BoundarySide side) const override
    {
        typedef HBSplinesFESpace<TDim-1, TLocalCoordinateType> BoundaryFESpaceType;
        typename BoundaryFESpaceType::Pointer pBoundaryFESpace = typename BoundaryFESpaceType::Pointer(new BoundaryFESpaceType());

        std::map<std::size_t, std::size_t> ident_indices_map;

        for (auto it = BaseType::bf_begin(); it != BaseType::bf_end(); ++it)
        {
            if (it->IsOnSide(BOUNDARY_FLAG(side)) && it->Is(ACTIVE))
            {
                typename BoundaryFESpaceType::bf_t pNewSubBf;

                if constexpr (TDim == 2)
                {
                    if ((side == _BLEFT_) || (side == _BRIGHT_))
                    {
                        pNewSubBf = it->Project(1);
                    }
                    if ((side == _BTOP_) || (side == _BBOTTOM_))
                    {
                        pNewSubBf = it->Project(0);
                    }
                }
                else if constexpr (TDim == 3)
                {
                    if ((side == _BFRONT_) || (side == _BBACK_))
                    {
                        pNewSubBf = it->Project(0);
                    }
                    if ((side == _BLEFT_) || (side == _BRIGHT_))
                    {
                        pNewSubBf = it->Project(1);
                    }
                    if ((side == _BTOP_) || (side == _BBOTTOM_))
                    {
                        pNewSubBf = it->Project(2);
                    }
                }

                pNewSubBf->Set(ACTIVE, true);
                pBoundaryFESpace->AddBf(pNewSubBf);
                ident_indices_map[pNewSubBf->EquationId()] = pNewSubBf->EquationId();
            }
        }

        // update the global to local map
        pBoundaryFESpace->UpdateFunctionIndices(ident_indices_map);

        // set the B-Splines information
        if constexpr (TDim == 2)
        {
            if ((side == _BLEFT_) || (side == _BRIGHT_))
            {
                pBoundaryFESpace->SetInfo(0, this->Order(1));
                pBoundaryFESpace->SetKnotVector(0, this->KnotVector(1));
            }
            if ((side == _BTOP_) || (side == _BBOTTOM_))
            {
                pBoundaryFESpace->SetInfo(0, this->Order(0));
                pBoundaryFESpace->SetKnotVector(0, this->KnotVector(0));
            }
        }
        else if constexpr (TDim == 3)
        {
            if ((side == _BFRONT_) || (side == _BBACK_))
            {
                pBoundaryFESpace->SetInfo(0, this->Order(1));
                pBoundaryFESpace->SetInfo(1, this->Order(2));
                pBoundaryFESpace->SetKnotVector(0, this->KnotVector(1));
                pBoundaryFESpace->SetKnotVector(1, this->KnotVector(2));
            }
            if ((side == _BLEFT_) || (side == _BRIGHT_))
            {
                pBoundaryFESpace->SetInfo(0, this->Order(2));
                pBoundaryFESpace->SetInfo(1, this->Order(0));
                pBoundaryFESpace->SetKnotVector(0, this->KnotVector(2));
                pBoundaryFESpace->SetKnotVector(1, this->KnotVector(0));
            }
            if ((side == _BTOP_) || (side == _BBOTTOM_))
            {
                pBoundaryFESpace->SetInfo(0, this->Order(0));
                pBoundaryFESpace->SetInfo(1, this->Order(1));
                pBoundaryFESpace->SetKnotVector(0, this->KnotVector(0));
                pBoundaryFESpace->SetKnotVector(1, this->KnotVector(1));
            }
        }

        pBoundaryFESpace->SetLastLevel(this->LastLevel());

        // construct the cells from the boundary basis functions
        typename BoundaryFESpaceType::cell_container_t::Pointer pnew_cells;

        pnew_cells = typename BoundaryFESpaceType::cell_container_t::Pointer(new BCellManager<TDim-1, typename BoundaryFESpaceType::CellType> ());

        if constexpr (TDim == 2)
        {
            TLocalCoordinateType tol = pBoundaryFESpace->KnotVector(0).GetResolution();

            for (auto it = pBoundaryFESpace->bf_begin(); it != pBoundaryFESpace->bf_end(); ++it)
            {
                for (std::size_t i1 = 0; i1 < pBoundaryFESpace->Order(0) + 1; ++i1)
                {
                    knot_t pXiMin = it->LocalKnots(0)[i1];
                    knot_t pXiMax = it->LocalKnots(0)[i1 + 1];

                    // check if the cell domain length is nonzero
                    TLocalCoordinateType length = (pXiMax->Value() - pXiMin->Value());
                    if (std::abs(length) > tol)
                    {
                        std::vector<knot_t> pKnots = {pXiMin, pXiMax};
                        typename BoundaryFESpaceType::cell_t pnew_cell = pBoundaryFESpace->pCellManager()->CreateCell(pKnots);
                        pnew_cell->SetLevel(this->LastLevel());
                        it->AddCell(pnew_cell);
                        pnew_cell->AddBf(*it.base());
                        pnew_cells->insert(pnew_cell);
                    }
                }
            }
        }
        else if constexpr (TDim == 3)
        {
            TLocalCoordinateType tol = std::sqrt(pBoundaryFESpace->KnotVector(0).GetResolution() * pBoundaryFESpace->KnotVector(1).GetResolution());

            for (auto it = pBoundaryFESpace->bf_begin(); it != pBoundaryFESpace->bf_end(); ++it)
            {
                for (std::size_t i1 = 0; i1 < pBoundaryFESpace->Order(0) + 1; ++i1)
                {
                    knot_t pXiMin = it->LocalKnots(0)[i1];
                    knot_t pXiMax = it->LocalKnots(0)[i1 + 1];

                    for (std::size_t j1 = 0; j1 < pBoundaryFESpace->Order(1) + 1; ++j1)
                    {
                        knot_t pEtaMin = it->LocalKnots(1)[j1];
                        knot_t pEtaMax = it->LocalKnots(1)[j1 + 1];

                        // check if the cell domain area is nonzero
                        TLocalCoordinateType area = (pXiMax->Value() - pXiMin->Value()) * (pEtaMax->Value() - pEtaMin->Value());
                        if (std::sqrt(std::abs(area)) > tol)
                        {
                            std::vector<knot_t> pKnots = {pXiMin, pXiMax, pEtaMin, pEtaMax};
                            typename BoundaryFESpaceType::cell_t pnew_cell = pBoundaryFESpace->pCellManager()->CreateCell(pKnots);
                            pnew_cell->SetLevel(this->LastLevel());
                            it->AddCell(pnew_cell);
                            pnew_cell->AddBf(*it.base());
                            pnew_cells->insert(pnew_cell);
                        }
                    }
                }
            }
        }

        // collapse the overlapping cells
        pBoundaryFESpace->pCellManager()->CollapseCells();

        // update the cell support and extraction operator
        pBoundaryFESpace->UpdateCells();

        // re-add the supporting cells
        for (auto it_cell = pBoundaryFESpace->pCellManager()->begin(); it_cell != pBoundaryFESpace->pCellManager()->end(); ++it_cell)
        {
            for (typename BoundaryFESpaceType::CellType::bf_iterator it_bf = (*it_cell)->bf_begin(); it_bf != (*it_cell)->bf_end(); ++it_bf)
            {
                (*it_bf).lock()->AddCell(*it_cell);
            }
        }

        return pBoundaryFESpace;
    }

    /// Construct the boundary FESpace based on side and rotation
    typename FESpace<TDim-1, TLocalCoordinateType>::Pointer ConstructBoundaryFESpace(const BoundarySide side,
            const std::map<std::size_t, std::size_t>& local_parameter_map,
            const std::vector<BoundaryDirection>& directions) const override
    {
        return this->ConstructBoundaryFESpace(side);
    }

    /// Overload assignment operator
    ThisType& operator=(const ThisType& rOther)
    {
        // TODO copy more
        KRATOS_ERROR << "The assignment oprator is not yet implemented";
        BaseType::operator=(rOther);
        return *this;
    }

    /// Clone this FESpace, this is a deep copy operation
    typename FESpaceType::Pointer Clone() const override
    {
        typename ThisType::Pointer pNewFESpace = typename ThisType::Pointer(new ThisType());
        *pNewFESpace = *this;
        return pNewFESpace;
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        BaseType::PrintInfo(rOStream);
        rOStream << "Number of levels = " << mLastLevel << std::endl;

        rOStream << "###############Begin knot vectors################" << std::endl;
        for (std::size_t lvl = 1; lvl <= this->LastLevel(); ++lvl)
        {
            rOStream << "Level " << lvl << std::endl;
            for (int dim = 0; dim < TDim; ++dim)
            {
                rOStream << " knot vector " << dim + 1 << ":";
                for (std::size_t i = 0; i < this->KnotVector(lvl, dim).size(); ++i)
                {
                    rOStream << " " << this->KnotVector(lvl, dim)[i];
                }
                rOStream << std::endl;
            }
        }
        rOStream << "###############End knot vectors##################" << std::endl;
    }

    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "###################" << std::endl;

        // print the cells in each level
        for (std::size_t level = 1; level < mLastLevel + 1; ++level)
        {
            rOStream << "###############Begin cells at level " << level << "################" << std::endl;
            std::size_t n = 0;
            for (typename cell_container_t::iterator it = BaseType::mpCellManager->begin(); it != BaseType::mpCellManager->end(); ++it)
            {
                if ((*it)->Level() == level)
                {
                    rOStream << "(" << ++n << ") " << *(*it) << std::endl;
                }
            }
            rOStream << "###############End cells at level " << level << "################" << std::endl;
        }
    }

private:

    std::size_t mLastLevel;
    std::size_t mMaxLevel;

    std::vector<boost::array<knot_container_t, TDim> > mKnotVectors;

    domain_container_t mSupportDomains; // this domain manager manages the support of all bfs in each level

    std::vector<std::size_t> mRefinementHistory;
};

/**
 * Template specific instantiation for null-D BSplines patch to terminate the compilation
 */
template<typename TLocalCoordinateType>
class HBSplinesFESpace<0, TLocalCoordinateType> : public FESpace<0, TLocalCoordinateType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HBSplinesFESpace);

    /// Type definition
    typedef FESpace<0, TLocalCoordinateType> BaseType;
    typedef KnotArray1D<TLocalCoordinateType> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;
    typedef HBSplinesBasisFunction<0> BasisFunctionType;
    typedef typename BasisFunctionType::Pointer bf_t;
    struct bf_compare { bool operator() (const bf_t& lhs, const bf_t& rhs) const {return lhs->Id() < rhs->Id();} };
    typedef std::set<bf_t, bf_compare> bf_container_t;
    typedef typename bf_container_t::iterator bf_iterator;
    typedef typename bf_container_t::const_iterator bf_const_iterator;

    typedef HBCell<HBSplinesBasisFunction<0> > CellType;
    typedef BaseBCellManager<CellType> cell_container_t;
    typedef typename cell_container_t::cell_t cell_t;

    /// Default constructor
    HBSplinesFESpace() : BaseType() {}

    /// Destructor
    ~HBSplinesFESpace() override {}

    /// Get the order of the HBSplines patch in specific direction
    std::size_t Order(std::size_t i) const override {return 0;}

    /// Get the number of basis functions defined over the HBSplinesFESpace
    std::size_t TotalNumber() const override {return 0;}

    /// Get the string describing the type of the patch
    std::string Type() const override
    {
        return StaticType();
    }

    /// Get the string describing the type of the patch
    static std::string StaticType()
    {
        return "HBSplinesFESpace0D";
    }

    /// Validate the HBSplinesFESpace before using
    bool Validate() const override
    {
        return BaseType::Validate();
    }

    /// Compare between two HBSplines patches in terms of parametric information
    bool IsCompatible(const BaseType& rOtherFESpace) const override
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
    void SetLastLevel(std::size_t LastLevel) {}

    /// Set the order for the B-Splines
    void SetInfo(std::size_t i, std::size_t order) {}

    /// Get the knot vector in i-direction, i=0..Dim
    /// User must be careful to use this function because it can modify the internal knot vectors
    knot_container_t& KnotVector(std::size_t i) {return mDummyKnotVector;}

    /// Get the knot vector in i-direction, i=0..Dim
    const knot_container_t& KnotVector(std::size_t i) const {return mDummyKnotVector;}

    /// Get the underlying cell manager
    typename cell_container_t::Pointer pCellManager() {return mpCellManager;}

    /// Get the underlying cell manager
    typename cell_container_t::ConstPointer pCellManager() const {return mpCellManager;}

    /// Update the function indices using a map. The map shall be the mapping from old index to new index.
    void UpdateFunctionIndices(const std::map<std::size_t, std::size_t>& indices_map) override {}

    /// Update the basis functions for all cells. This function must be called before any operation on cell is required.
    void UpdateCells() {}

private:

    knot_container_t mDummyKnotVector;
    bf_container_t mpBasisFuncs;
    typename cell_container_t::Pointer mpCellManager;
};

/// output stream function
template<int TDim, typename TLocalCoordinateType>
inline std::ostream& operator <<(std::ostream& rOStream, const HBSplinesFESpace<TDim, TLocalCoordinateType>& rThis)
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

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_FESPACE_H_INCLUDED defined
