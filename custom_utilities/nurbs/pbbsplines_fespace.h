//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 24 Nov 2017 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_PBBSPLINES_FESPACE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_PBBSPLINES_FESPACE_H_INCLUDED

// System includes
#include <vector>

// External includes
#include <boost/array.hpp>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "containers/pointer_vector_set.h"
#include "custom_utilities/iga_define.h"
#include "custom_utilities/fespace.h"
#include "isogeometric_application_variables.h"

#define DEBUG_GEN_CELL

namespace Kratos
{

/**
 * Abstract class represents the FESpace for a single point-based B-Splines patch.
 */
template<int TDim, typename TLocalCoordinateType, typename TBasisFunctionType, typename TCellManagerType>
class PBBSplinesFESpace : public FESpace<TDim, TLocalCoordinateType>, public IsogeometricEcho
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PBBSplinesFESpace);

    /// Type definition
    typedef FESpace<TDim, TLocalCoordinateType> BaseType;
    typedef PBBSplinesFESpace<TDim, TLocalCoordinateType, TBasisFunctionType, TCellManagerType> ThisType;
    typedef TBasisFunctionType BasisFunctionType;
    typedef typename BasisFunctionType::Pointer bf_t;
    struct bf_compare { bool operator() (const bf_t& lhs, const bf_t& rhs) const {return lhs->Id() < rhs->Id();} };
    typedef PointerVectorSet<BasisFunctionType, IndexedObject> bf_container_t;
    typedef typename bf_container_t::iterator bf_iterator;
    typedef typename bf_container_t::const_iterator bf_const_iterator;

    typedef typename TBasisFunctionType::CellType CellType;
    typedef TCellManagerType cell_container_t;
    typedef typename cell_container_t::cell_t cell_t;
    typedef typename CellType::knot_container_t knot_container_t;
    typedef typename CellType::knot_t knot_t;

    typedef std::map<std::size_t, bf_t> function_map_t;

    /// Default constructor
    PBBSplinesFESpace() : BaseType()
    {
        mpCellManager = typename cell_container_t::Pointer(new TCellManagerType());
    }

    /// Destructor
    ~PBBSplinesFESpace() override
    {
#ifdef ISOGEOMETRIC_DEBUG_DESTROY
        std::cout << Type() << ", Addr = " << this << " is destroyed" << std::endl;
#endif
    }

    /// Helper to create new PBBSplinesFESpace pointer
    static typename ThisType::Pointer Create()
    {
        return typename ThisType::Pointer(new ThisType());
    }

    /// Add a already constructed basis function to the internal list
    void AddBf(bf_t p_bf)
    {
        mpBasisFuncs.insert(p_bf);
    }

    /// Check if the bf exists in the list; otherwise create new bf and return
    /// The Id is only used when the new bf is created. User must always check the Id of the returned function.
    bf_t CreateBf(std::size_t Id, const std::vector<std::vector<knot_t> >& rpKnots)
    {
        // search in the current list of basis functions, the one that has the same local knot vector with provided ones
        for (bf_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if (it->CompareLocalKnots(rpKnots))
            {
                return *it.base();
            }
        }

        // create the new bf and add the knot
        bf_t p_bf = bf_t(new BasisFunctionType(Id));
        for (int dim = 0; dim < TDim; ++dim)
        {
            p_bf->SetLocalKnotVectors(dim, rpKnots[dim]);
            p_bf->SetInfo(dim, this->Order(dim));
        }
        mpBasisFuncs.insert(mpBasisFuncs.end(), p_bf);

        return p_bf;
    }

    /// Remove the basis functions from the container
    void RemoveBf(bf_t p_bf)
    {
        mpBasisFuncs.erase(p_bf->Id());
    }

    // Iterators for the basis functions
    bf_iterator bf_begin() {return mpBasisFuncs.begin();}
    bf_const_iterator bf_begin() const {return mpBasisFuncs.begin();}
    bf_iterator bf_end() {return mpBasisFuncs.end();}
    bf_const_iterator bf_end() const {return mpBasisFuncs.end();}

    /// Get the last id of the basis functions
    std::size_t LastId() const {bf_const_iterator it = bf_end(); --it; return it->Id();}

    /// Set the order for the point-based B-Splines
    void SetInfo(std::size_t i, std::size_t order) {mOrders[i] = order;}

    /// Get the order of the BSplines patch in specific direction
    std::size_t Order(std::size_t i) const override
    {
        if (i >= TDim) { return 0; }
        else { return mOrders[i]; }
    }

    /// Get the number of basis functions defined over the BSplines
    std::size_t TotalNumber() const override {return mpBasisFuncs.size();}

    /// Get the lower and upper bound of the parametric space in a specific direction
    std::vector<double> ParametricBounds(std::size_t di) const override
    {
        std::vector<double> bound = {1.0e99, -1.0e99};

        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            std::vector<double> bb = it->GetBoundingBox();
            if (bound[0] > bb[2 * di]) { bound[0] = bb[2 * di]; }
            if (bound[1] < bb[2 * di + 1]) { bound[1] = bb[2 * di + 1]; }
        }

        return bound;
    }

    /// Get the weights of all the basis functions
    std::vector<double> GetWeights() const
    {
        std::vector<double> weights(this->TotalNumber());
        std::size_t cnt = 0;
        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it, ++cnt)
        {
            weights[cnt] = it->GetData(CONTROL_POINT).W();
        }
        return weights;
    }

    /// Get the string representing the type of the patch
    std::string Type() const override
    {
        return StaticType();
    }

    /// Get the string representing the type of the patch
    static std::string StaticType()
    {
        std::stringstream ss;
        ss << "PBBSplinesFESpace" << TDim << "D";
        return ss.str();
    }

    /// Validate the PBBSplinesFESpace
    bool Validate() const override
    {
        // TODO validate more
        return BaseType::Validate();
    }

    /// Get the values of the basis function i at point xi
    /// REMARK: This function only returns the unweighted basis function value. To obtain the correct one, use WeightedFESpace
    void GetValue(double& v, std::size_t i, const std::vector<double>& xi) const override
    {
        std::size_t j = 0;
        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if (j == i)
            {
                v = it->GetValueAt(xi);
                return;
            }
            ++j;
        }
        v = 0.0;
    }

    /// Get the values of the basis functions at point xi
    /// REMARK: This function only returns the unweighted basis function value. To obtain the correct one, use WeightedFESpace
    void GetValues(std::vector<double>& values, const std::vector<double>& xi) const override
    {
        if (values.size() != this->TotalNumber())
        {
            values.resize(this->TotalNumber());
        }
        std::size_t i = 0;
        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            values[i++] = it->GetValueAt(xi);
        }
    }

    /// Get the derivative of the basis function i at point xi
    /// the output derivatives has the form of values[func_index][dim_index]
    /// REMARK: This function only returns the unweighted basis function derivatives. To obtain the correct one, use WeightedFESpace
    void GetDerivative(std::vector<double>& values, std::size_t i, const std::vector<double>& xi) const override
    {
        std::size_t j = 0;
        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if (j == i)
            {
                it->GetDerivativeAt(values, xi);
                return;
            }
            ++j;
        }
        if (values.size() != TDim)
        {
            values.resize(TDim);
        }
        for (int dim = 0; dim < TDim; ++dim) { values[dim] = 0.0; }
    }

    /// Get the derivative of the basis functions at point xi
    /// the output derivatives has the form of values[func_index][dim_index]
    /// REMARK: This function only returns the unweighted basis function derivatives. To obtain the correct one, use WeightedFESpace
    void GetDerivatives(std::vector<std::vector<double> >& values, const std::vector<double>& xi) const override
    {
        if (values.size() != this->TotalNumber())
        {
            values.resize(this->TotalNumber());
        }
        std::size_t i = 0;
        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if (values[i].size() != TDim)
            {
                values[i].resize(TDim);
            }
            it->GetDerivativeAt(values[i], xi);
            ++i;
        }
    }

    /// Get the values and derivatives of the basis functions at point xi
    /// the output derivatives has the form of values[func_index][dim_index]
    /// REMARK: This function only returns the unweighted basis function derivatives. To obtain the correct one, use WeightedFESpace
    void GetValuesAndDerivatives(std::vector<double>& values, std::vector<std::vector<double> >& derivatives,
                                 const std::vector<double>& xi) const override
    {
        if (values.size() != this->TotalNumber())
        {
            values.resize(this->TotalNumber());
        }
        if (derivatives.size() != this->TotalNumber())
        {
            derivatives.resize(this->TotalNumber());
        }

        std::size_t i = 0;
        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            values[i] = it->GetValueAt(xi);
            if (derivatives[i].size() != TDim)
            {
                derivatives[i].resize(TDim);
            }
            it->GetDerivativeAt(derivatives[i], xi);
            ++i;
        }
    }

    /// Check if a point lies inside the parametric domain of the BSplinesFESpace
    bool IsInside(const std::vector<double>& xi) const override
    {
        // TODO
        KRATOS_ERROR << "Not yet implemented";
    }

    /// Compare between two BSplines patches in terms of parametric information
    bool IsCompatible(const BaseType& rOtherFESpace) const override
    {
        if (rOtherFESpace.Type() != Type())
        {
            KRATOS_WATCH(rOtherFESpace.Type())
            KRATOS_WATCH(Type())
            std::cout << "WARNING!!! the other patch type is not " << Type() << std::endl;
            return false;
        }

        const ThisType* pOtherPBBSplinesFESpace = dynamic_cast<const ThisType*>(&rOtherFESpace);
        if (pOtherPBBSplinesFESpace == nullptr)
            return false;

        // compare the knot vectors and order information
        for (std::size_t i = 0; i < TDim; ++i)
        {
            if (!(this->Order(i)) == pOtherPBBSplinesFESpace->Order(i))
            {
                return false;
            }
        }

        return true;
    }

    /// Reset the function indices to -1.
    /// This is useful when assigning the id for the boundary patch.
    void ResetFunctionIndices() override
    {
        BaseType::mGlobalToLocal.clear();
        for (bf_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            it->SetEquationId(-1);
        }
    }

    /// Reset the function indices to a given values.
    /// This is useful when assigning the id for the boundary patch.
    void ResetFunctionIndices(const std::vector<std::size_t>& func_indices) override
    {
        if (func_indices.size() != this->TotalNumber())
        {
            // KRATOS_WATCH(this->TotalNumber())
            std::cout << "func_indices:";
            for (std::size_t i = 0; i < func_indices.size(); ++i)
            {
                std::cout << " " << func_indices[i];
            }
            std::cout << std::endl;
            KRATOS_ERROR << "The func_indices vector does not have the same size as total number of basis functions";
        }
        std::size_t cnt = 0;
        BaseType::mGlobalToLocal.clear();
        for (bf_iterator it = bf_begin(); it != bf_end(); ++it, ++cnt)
        {
            it->SetEquationId(func_indices[cnt]);
            BaseType::mGlobalToLocal[it->EquationId()] = cnt;
        }
    }

    /// Enumerate the dofs of each grid function. This function is used to initialize the equation id for all basis functions.
    /// If the dof does not have pre-existing value, which assumes -1, it will be assigned the incremental value.
    /// This function can be used after refinement, providing that all the -1 are overrided, and the start value must be the latest one on the multipatch.
    std::size_t& Enumerate(std::size_t& start) override
    {
        BaseType::mGlobalToLocal.clear();
        std::size_t cnt = 0;
        for (bf_iterator it = bf_begin(); it != bf_end(); ++it, ++cnt, ++start)
        {
            if (it->EquationId() == -1) { it->SetEquationId(start); }
            BaseType::mGlobalToLocal[it->EquationId()] = cnt;
        }

        return start;
    }

    /// Access the function indices (aka global ids)
    std::vector<std::size_t> FunctionIndices() const override
    {
        std::vector<std::size_t> func_indices(this->TotalNumber());
        std::size_t cnt = 0;
        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it, ++cnt)
        {
            func_indices[cnt] = it->EquationId();
        }
        return func_indices;
    }

    /// Update the function indices using a map. The map shall be the mapping from old index to new index.
    void UpdateFunctionIndices(const std::map<std::size_t, std::size_t>& indices_map) override
    {
        std::size_t cnt = 0;
        BaseType::mGlobalToLocal.clear();
        for (bf_iterator it = bf_begin(); it != bf_end(); ++it, ++cnt)
        {
            std::map<std::size_t, std::size_t>::const_iterator it2 = indices_map.find(it->EquationId());

            if (it2 == indices_map.end())
            {
                std::cout << "WARNING!!! the indices_map does not contain " << it->EquationId() << std::endl;
                continue;
            }

            it->SetEquationId(it2->second);
            BaseType::mGlobalToLocal[it->EquationId()] = cnt;
        }
    }

    /// Get the first equation_id in this space
    std::size_t GetFirstEquationId() const override
    {
        std::size_t first_id;

        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if (it->EquationId() == -1)
            {
                return -1;
            }
            else
            {
                if (it == bf_begin())
                {
                    first_id = it->EquationId();
                }
                else if (it->EquationId() < first_id)
                {
                    first_id = it->EquationId();
                }
            }
        }

        return first_id;
    }

    /// Get the last equation_id in this space
    std::size_t GetLastEquationId() const override
    {
        std::size_t last_id = -1;
        bool hit = false;

        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if (it->EquationId() != -1)
            {
                if (!hit)
                {
                    hit = true;
                    last_id = it->EquationId();
                }
                else
                {
                    if (it->EquationId() > last_id)
                    {
                        last_id = it->EquationId();
                    }
                }
            }
        }

        return last_id;
    }

    /// Get the basis functions based on boundary flag. This allows to extract the corner bf.
    std::vector<bf_t> ExtractBoundaryBfsByFlag(std::size_t boundary_id) const
    {
        std::vector<bf_t> bf_list;

        for (bf_const_iterator it_bf = bf_begin(); it_bf != bf_end(); ++it_bf)
        {
            if (it_bf->IsOnSide(boundary_id))
            {
                bf_list.push_back(*it_bf.base());
            }
        }

        return bf_list;
    }

    /// Extract the index of the functions on the boundaries
    std::vector<std::size_t> ExtractBoundaryFunctionIndicesByFlag(int boundary_id) const override
    {
        std::vector<bf_t> bfs = this->ExtractBoundaryBfsByFlag(boundary_id);

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
        std::vector<std::size_t> func_indices;

        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if (it->IsOnSide(BOUNDARY_FLAG(side)))
            {
                func_indices.push_back(it->EquationId());
            }
        }

        return func_indices;
    }

    /// Assign the index for the functions on the boundary
    void AssignBoundaryFunctionIndices(const BoundarySide side, const std::vector<std::size_t>& func_indices, const bool override) override
    {
        std::vector<bf_t> bf_list;
        for (bf_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            if (it->IsOnSide(BOUNDARY_FLAG(side)))
            {
                bf_list.push_back(*it.base());
            }
        }

        std::size_t cnt = 0;
        for (auto it = bf_list.begin(); it != bf_list.end(); ++it, ++cnt)
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

    /// Get the underlying cell manager
    typename cell_container_t::Pointer pCellManager() {return mpCellManager;}

    /// Get the underlying cell manager
    typename cell_container_t::ConstPointer pCellManager() const {return mpCellManager;}

    /// Clean the internal data of all the cells
    void ResetCells()
    {
        for (typename cell_container_t::iterator it_cell = mpCellManager->begin(); it_cell != mpCellManager->end(); ++it_cell)
        {
            (*it_cell)->Reset();
        }
    }

    /// Update the basis functions for all cells. This function must be called before any operation on cell is required.
    virtual void UpdateCells()
    {
        this->ResetCells();

        // for each cell compute the extraction operator and add to the anchor
        Vector Crow;
        for (bf_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            for (typename BasisFunctionType::cell_iterator it_cell = it->cell_begin(); it_cell != it->cell_end(); ++it_cell)
            {
                it->ComputeExtractionOperator(Crow, *it_cell);
                (*it_cell)->AddAnchor(it->EquationId(), it->GetData(CONTROL_POINT).W(), Crow);
            }
        }
    }

    /// Create the cell manager for all the cells in the support domain of the PBBSplinesFESpace
    typename BaseType::cell_container_t::Pointer ConstructCellManager() const override
    {
        return mpCellManager;
    }

    /// Overload operator[], this allows to access the basis function randomly based on index
    bf_t operator[](std::size_t i)
    {
        typename bf_container_t::iterator it = mpBasisFuncs.begin();
        std::advance(it, i);
        return *it.base();
    }

    /// Overload operator(), this allows to access the basis function based on its id
    bf_t operator()(std::size_t Id)
    {
        return mpBasisFuncs(Id);
    }

    /// Check if the functional space has the function with equation id
    bool HasBfByEquationId(std::size_t EquationId) const
    {
        return (BaseType::mGlobalToLocal.find(EquationId) != BaseType::mGlobalToLocal.end());
    }

    /// Check if the functional space has the function with Id
    bool HasBfById(std::size_t Id) const
    {
        return (mpBasisFuncs.find(Id) != mpBasisFuncs.end());
    }

    /// Get the basis function by equation id
    bf_t pGetBfByEquationId(std::size_t EquationId)
    {
        return this->operator[](BaseType::mGlobalToLocal[EquationId]);
    }

    /// Overload assignment operator
    ThisType& operator=(const ThisType& rOther)
    {
        // TODO copy more
        KRATOS_ERROR << "The assignment operator is not fully implemented";
        BaseType::operator=(rOther);
        return *this;
    }

    /// Clone this FESpace, this is a deep copy operation
    typename BaseType::Pointer Clone() const override
    {
        typename ThisType::Pointer pNewFESpace = typename ThisType::Pointer(new ThisType());
        *pNewFESpace = *this;
        return pNewFESpace;
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Type() << ", Addr = " << this << ", n = " << this->TotalNumber();
        rOStream << ", p = (";
        for (std::size_t dim = 0; dim < TDim; ++dim)
        {
            rOStream << " " << this->Order(dim);
        }
        rOStream << ")";
    }

    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "###################" << std::endl;

        // print the basis functions
        rOStream << "Basis functions:" << std::endl;
        for (bf_const_iterator it = bf_begin(); it != bf_end(); ++it)
        {
            rOStream << " ++ " << (*it) << std::endl;
        }

        // print the cells
        rOStream << "Cells:" << std::endl;
        for (typename cell_container_t::iterator it = mpCellManager->begin(); it != mpCellManager->end(); ++it)
        {
            rOStream << " ++ " << *(*it) << std::endl;
        }
    }

protected:

    boost::array<std::size_t, TDim> mOrders;

    typename cell_container_t::Pointer mpCellManager;

    bf_container_t mpBasisFuncs;

private:

    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
        rSerializer.save( "Orders", mOrders );
        // TODO
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
        rSerializer.load( "Orders", mOrders );
        // TODO
    }

    ///@}
};

/// output stream function
template<int TDim, typename TLocalCoordinateType, typename TBasisFunctionType, typename TCellManagerType>
inline std::ostream& operator <<(std::ostream& rOStream, const PBBSplinesFESpace<TDim, TLocalCoordinateType, TBasisFunctionType, TCellManagerType>& rThis)
{
    rOStream << "-------------Begin PBBSplinesFESpace Info-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "-------------End PBBSplinesFESpace Info-------------" << std::endl;
    return rOStream;
}

} // namespace Kratos.

#undef DEBUG_GEN_CELL

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_PBBSPLINES_FESPACE_H_INCLUDED defined
