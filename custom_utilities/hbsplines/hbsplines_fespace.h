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
#include "custom_utilities/fespace.h"
#include "custom_utilities/nurbs/knot_array_1d.h"
#include "custom_utilities/nurbs/bsplines_indexing_utility.h"
#include "custom_utilities/nurbs/cell_manager_2d.h"
#include "custom_utilities/nurbs/cell_manager_3d.h"
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
        if (TDim == 2)
        {
            mpCellManager = typename cell_container_t::Pointer(new CellManager2D<CellType>());
        }
        else if(TDim == 3)
        {
            mpCellManager = typename cell_container_t::Pointer(new CellManager3D<CellType>());
        }
    }

    /// Destructor
    virtual ~HBSplinesFESpace() {}

    /// Helper to create new HBSplinesFESpace pointer
    static HBSplinesFESpace<TDim>::Pointer Create()
    {
        return HBSplinesFESpace<TDim>::Pointer(new HBSplinesFESpace());
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
        {
            weights[cnt++] = (*it)->GetValue(CONTROL_POINT).W();
        }
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
    virtual double GetValue(const std::size_t& i, const std::vector<double>& xi) const
    {
        std::vector<double> weights = this->GetWeights();

        double sum = 0.0;
        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
            sum += weights[i]*(*it)->GetValue(xi);

        return weights[i]/sum;
    }

    /// Get the values of the basis functions at point xi
    virtual std::vector<double> GetValue(const std::vector<double>& xi) const
    {
        std::vector<double> values(this->TotalNumber());
        std::vector<double> weights = this->GetWeights();

        double sum = 0.0;
        std::size_t i = 0;
        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            values[i] = (*it)->GetValue(xi);
            sum += weights[i]*values[i];
            ++i;
        }

        for (i = 0; i < values.size(); ++i)
            values[i] *= weights[i]/sum;

        return values;
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

    /// Enumerate the dofs of each grid function. The enumeration algorithm is pretty straightforward.
    /// If the dof does not have pre-existing value, which assume it is -1, it will be assigned the incremental value.
    virtual std::size_t& Enumerate(std::size_t& start)
    {
        // enumerate all basis functions
        BaseType::mGlobalToLocal.clear();
        BaseType::mFunctionsIds.resize(this->TotalNumber());
        std::size_t cnt = 0;
        for (bf_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            (*it)->SetId(start++);
            BaseType::mFunctionsIds[cnt] = (*it)->Id();
            BaseType::mGlobalToLocal[BaseType::mFunctionsIds[cnt]] = cnt;
            ++cnt;
        }

        return start;
    }

    /// Extract the index of the functions on the boundary
    virtual std::vector<std::size_t> ExtractBoundaryFunctionIndices(const BoundarySide& side) const
    {
        std::vector<std::size_t> func_indices;

        // TODO

        return func_indices;
    }

    /// Assign the index for the functions on the boundary
    virtual void AssignBoundaryFunctionIndices(const BoundarySide& side, const std::vector<std::size_t>& func_indices)
    {
        // TODO
    }

    /// Construct the boundary patch based on side
    virtual typename FESpace<TDim-1>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    {
        typename HBSplinesFESpace<TDim-1>::Pointer pBFESpace = typename HBSplinesFESpace<TDim-1>::Pointer(new HBSplinesFESpace<TDim-1>());

        // TODO construct the basis function and cells on the boundary

        // transfer the function indices
        std::vector<std::size_t> b_func_indices = ExtractBoundaryFunctionIndices(side);
        pBFESpace->ResetFunctionIndices(b_func_indices);

        return pBFESpace;
    }

    /// Get the underlying cell manager
    typename cell_container_t::Pointer pCellManager() {return mpCellManager;}

    /// Get the underlying cell manager
    typename cell_container_t::ConstPointer pCellManager() const {return mpCellManager;}

    /// Create the cell manager for all the cells in the support domain of the HBSplinesFESpace
    virtual typename BaseType::cell_container_t::Pointer ConstructCellManager() const
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

        // create the compatible cell manager and add to the list
        typename BaseType::cell_container_t::Pointer pCompatCellManager;
        if (TDim == 2)
            pCompatCellManager = CellManager2D<typename BaseType::cell_container_t::CellType>::Create();
        else if (TDim == 3)
            pCompatCellManager = CellManager3D<typename BaseType::cell_container_t::CellType>::Create();
        for(typename cell_container_t::iterator it_cell = mpCellManager->begin(); it_cell != mpCellManager->end(); ++it_cell)
        {
            pCompatCellManager->insert(*it_cell);
        }

        return pCompatCellManager;
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
        // TODO compare more
        BaseType::operator=(rOther);
        return *this;
    }

    /// Clone this FESpace, this is a deep copy operation
    virtual typename FESpace<TDim>::Pointer Clone() const
    {
        typename HBSplinesFESpace<TDim>::Pointer pNewFESpace = typename HBSplinesFESpace<TDim>::Pointer(new HBSplinesFESpace<TDim>());
        // TODO copy more
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

#undef DEBUG_GEN_CELL

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_FESPACE_H_INCLUDED defined
