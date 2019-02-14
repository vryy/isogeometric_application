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
template<int TDim>
class HBSplinesFESpace : public PBBSplinesFESpace<TDim, HBSplinesBasisFunction<TDim>, BCellManager<TDim, typename HBSplinesBasisFunction<TDim>::CellType> >
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(HBSplinesFESpace);

    /// Type definition
    typedef PBBSplinesFESpace<TDim, HBSplinesBasisFunction<TDim>, BCellManager<TDim, typename HBSplinesBasisFunction<TDim>::CellType> > BaseType;
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

    typedef typename BaseType::function_map_t function_map_t;

    /// Default constructor
    HBSplinesFESpace() : BaseType(), mLastLevel(1), mMaxLevel(10)
    {}

    /// Destructor
    virtual ~HBSplinesFESpace()
    {
        #ifdef ISOGEOMETRIC_DEBUG_DESTROY
        std::cout << Type() << ", Addr = " << this << " is destroyed" << std::endl;
        #endif
    }

    /// Helper to create new HBSplinesFESpace pointer
    static typename HBSplinesFESpace<TDim>::Pointer Create()
    {
        return typename HBSplinesFESpace<TDim>::Pointer(new HBSplinesFESpace());
    }

    /// Check if the bf exists in the list; otherwise create new bf and return
    /// The Id is only used when the new bf is created. User must always check the Id of the returned function.
    bf_t CreateBf(const std::size_t& Id, const std::size_t& Level, const std::vector<std::vector<knot_t> >& rpKnots)
    {
        // search in the current list of basis functions, the one that has the same local knot vector with provided ones
        for(bf_iterator it = BaseType::bf_begin(); it != BaseType::bf_end(); ++it)
            if((*it)->Contain(rpKnots))
                return *it;

        // create the new bf and add the knot
        bf_t p_bf = bf_t(new BasisFunctionType(Id, Level));
        for (int dim = 0; dim < TDim; ++dim)
        {
            p_bf->SetLocalKnotVectors(dim, rpKnots[dim]);
            p_bf->SetInfo(dim, this->Order(dim));
        }
        BaseType::mpBasisFuncs.insert(p_bf);
        BaseType::m_function_map_is_created = false;

        return p_bf;
    }

    /// EtaMaxdate the basis functions for all cells. This function must be called before any operation on cell is required.
    virtual void UpdateCells()
    {
        this->ResetCells();

        // for each cell compute the extraction operator and add to the anchor
        Vector Crow;
        for(typename cell_container_t::iterator it_cell = BaseType::mpCellManager->begin(); it_cell != BaseType::mpCellManager->end(); ++it_cell)
        {
            for(typename CellType::bf_iterator it_bf = (*it_cell)->bf_begin(); it_bf != (*it_cell)->bf_end(); ++it_bf)
            {
                BasisFunctionType& bf = *(it_bf->lock());
                bf.ComputeExtractionOperator(Crow, *it_cell);
                (*it_cell)->AddAnchor(bf.EquationId(), bf.GetValue(CONTROL_POINT).W(), Crow);
            }
        }
    }

    /// Get the knot vector in i-direction, i=0..Dim
    /// User must be careful to use this function because it can modify the internal knot vectors
    knot_container_t& KnotVector(const std::size_t& i) {return mKnotVectors[i];}

    /// Get the knot vector in i-direction, i=0..Dim
    const knot_container_t& KnotVector(const std::size_t& i) const {return mKnotVectors[i];}

    /// Get the last refinement level ain the hierarchical mesh
    const std::size_t& LastLevel() const {return mLastLevel;}

    /// Set the last level in the hierarchical mesh
    void SetLastLevel(const std::size_t& LastLevel) {mLastLevel = LastLevel;}

    /// Get the maximum level allowed in the hierarchical mesh
    const std::size_t& MaxLevel() const {return mMaxLevel;}

    /// Set the maximum level allowed in the hierarchical mesh
    void SetMaxLevel(const std::size_t& MaxLevel) {mMaxLevel = MaxLevel;}

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

    /// Construct the boundary FESpace based on side
    virtual typename FESpace<TDim-1>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    {
        typedef HBSplinesFESpace<TDim-1> BoundaryFESpaceType;
        typename BoundaryFESpaceType::Pointer pBFESpace = typename BoundaryFESpaceType::Pointer(new BoundaryFESpaceType());

        std::map<std::size_t, std::size_t> ident_indices_map;

        for (bf_iterator it = BaseType::bf_begin(); it != BaseType::bf_end(); ++it)
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
                ident_indices_map[pNewSubBf->EquationId()] = pNewSubBf->EquationId();
            }
        }

        // update the global to local map
        pBFESpace->UpdateFunctionIndices(ident_indices_map);

        // set the B-Splines information
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

        pnew_cells = typename BoundaryFESpaceType::cell_container_t::Pointer(new BCellManager<TDim-1, typename BoundaryFESpaceType::CellType>());

        if (TDim == 2)
        {
            for (typename BoundaryFESpaceType::bf_iterator it = pBFESpace->bf_begin(); it != pBFESpace->bf_end(); ++it)
            {
                for(std::size_t i1 = 0; i1 < pBFESpace->Order(0) + 1; ++i1)
                {
                    knot_t pXiMin = (*it)->LocalKnots(0)[i1];
                    knot_t pXiMax = (*it)->LocalKnots(0)[i1 + 1];

                    // check if the cell domain length is nonzero
                    double length = (pXiMax->Value() - pXiMin->Value());
                    if(fabs(length) > cell_tol)
                    {
                        std::vector<knot_t> pKnots = {pXiMin, pXiMax};
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
            for (typename BoundaryFESpaceType::bf_iterator it = pBFESpace->bf_begin(); it != pBFESpace->bf_end(); ++it)
            {
                for(std::size_t i1 = 0; i1 < pBFESpace->Order(0) + 1; ++i1)
                {
                    knot_t pXiMin = (*it)->LocalKnots(0)[i1];
                    knot_t pXiMax = (*it)->LocalKnots(0)[i1 + 1];

                    for(std::size_t j1 = 0; j1 < pBFESpace->Order(1) + 1; ++j1)
                    {
                        knot_t pEtaMin = (*it)->LocalKnots(1)[j1];
                        knot_t pEtaMax = (*it)->LocalKnots(1)[j1 + 1];

                        // check if the cell domain area is nonzero
                        double area = (pXiMax->Value() - pXiMin->Value()) * (pEtaMax->Value() - pEtaMin->Value());
                        if(sqrt(fabs(area)) > cell_tol)
                        {
                            std::vector<knot_t> pKnots = {pXiMin, pXiMax, pEtaMin, pEtaMax};
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

        // update the cell support and extraction operator
        pBFESpace->UpdateCells();

        // re-add the supporting cells
        for(typename BoundaryFESpaceType::cell_container_t::iterator it_cell = pBFESpace->pCellManager()->begin(); it_cell != pBFESpace->pCellManager()->end(); ++it_cell)
        {
            for(typename BoundaryFESpaceType::CellType::bf_iterator it_bf = (*it_cell)->bf_begin(); it_bf != (*it_cell)->bf_end(); ++it_bf)
            {
                (*it_bf).lock()->AddCell(*it_cell);
            }
        }

        return pBFESpace;
    }

    /// Construct the boundary FESpace based on side and rotation
    virtual typename FESpace<TDim-1>::Pointer ConstructBoundaryFESpace(const BoundarySide& side,
        const std::map<std::size_t, std::size_t>& local_parameter_map, const std::vector<BoundaryDirection>& directions) const
    {
        return this->ConstructBoundaryFESpace(side);
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
        BaseType::PrintInfo(rOStream);
        rOStream << "Number of levels = " << mLastLevel << std::endl;

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
        BaseType::PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "###################" << std::endl;

        // print the cells in each level
        for (std::size_t level = 1; level < mLastLevel+1; ++level)
        {
            rOStream << "###############Begin cells at level " << level << "################" << std::endl;
            std::size_t n = 0;
            for(typename cell_container_t::iterator it = BaseType::mpCellManager->begin(); it != BaseType::mpCellManager->end(); ++it)
                if((*it)->Level() == level)
                    rOStream << "(" << ++n << ") " << *(*it) << std::endl;
            rOStream << "###############End cells at level " << level << "################" << std::endl;
        }
    }

private:

    std::size_t mLastLevel;
    std::size_t mMaxLevel;

    boost::array<knot_container_t, TDim> mKnotVectors;

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
    typedef BaseBCellManager<CellType> cell_container_t;
    typedef typename cell_container_t::cell_t cell_t;

    /// Default constructor
    HBSplinesFESpace() : BaseType() {}

    /// Destructor
    virtual ~HBSplinesFESpace() {}

    /// Get the order of the BSplines patch in specific direction
    virtual std::size_t Order(const std::size_t& i) const {return 0;}

    /// Get the number of basis functions defined over the BSplines HBSplinesFESpace
    virtual std::size_t Number() const {return 0;}

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

    /// Update the function indices using a map. The map shall be the mapping from old index to new index.
    virtual void UpdateFunctionIndices(const std::map<std::size_t, std::size_t>& indices_map) {}

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

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_FESPACE_H_INCLUDED defined
