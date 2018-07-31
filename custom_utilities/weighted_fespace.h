//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2 Dec 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_WEIGHTED_WeightedFESpace_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_WEIGHTED_WeightedFESpace_H_INCLUDED

// System includes
#include <vector>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "custom_utilities/fespace.h"


namespace Kratos
{

/**
 * A weighted FESpace add the weighted information to the FESpace.
 */
template<int TDim>
class WeightedFESpace : public FESpace<TDim>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(WeightedFESpace);

    /// Type definition
    typedef FESpace<TDim> BaseType;
    typedef typename BaseType::cell_container_t cell_container_t;

    /// Default constructor
    WeightedFESpace(typename BaseType::Pointer pFESpace, const std::vector<double>& weights)
    : BaseType(), mpFESpace(pFESpace)
    {
        mWeights.resize(weights.size());
        std::copy(weights.begin(), weights.end(), mWeights.begin());
    }

    /// Destructor
    virtual ~WeightedFESpace()
    {
    }

    /// Helper to create new WeightedFESpace pointer
    static typename WeightedFESpace<TDim>::Pointer Create(typename BaseType::Pointer pFESpace, const std::vector<double>& weights)
    {
        return typename WeightedFESpace<TDim>::Pointer(new WeightedFESpace(pFESpace, weights));
    }

    /// Get the number of basis functions defined over the WeightedFESpace
    virtual const std::size_t TotalNumber() const
    {
        return mpFESpace->TotalNumber();
    }

    /// Get the order of the WeightedFESpace in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
    {
        return mpFESpace->Order(i);
    }

    /// Set the weight vector
    void SetWeights(const std::vector<double>& weights)
    {
        if (mWeights.size() != weights.size())
            mWeights.resize(weights.size());
        std::copy(weights.begin(), weights.end(), mWeights.begin());
    }

    /// Get the weight vector
    const std::vector<double>& Weights() const {return mWeights;}

    /// Get the string representing the type of the WeightedFESpace
    virtual std::string Type() const
    {
        std::stringstream ss;
        ss << StaticType() << "_over_" << mpFESpace->Type();
        return ss.str();
    }

    /// Get the string representing the type of the WeightedFESpace
    static std::string StaticType()
    {
        std::stringstream ss;
        ss << "WeightedFESpace" << TDim << "D";
        return ss.str();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Get the values of the basis function i at point xi
    virtual void GetValue(double& v, const std::size_t& i, const std::vector<double>& xi) const
    {
        std::vector<double> values;
        mpFESpace->GetValue(values, xi);
        double sum_value = 0.0;
        for (std::size_t j = 0; j < values.size(); ++j)
            sum_value += mWeights[j] * values[j];
        v = mWeights[i]*values[i] / sum_value;
    }

    /// Get the values of the basis functions at point xi
    virtual void GetValue(std::vector<double>& new_values, const std::vector<double>& xi) const
    {
        std::vector<double> values;
        mpFESpace->GetValue(values, xi);

        if (new_values.size() != values.size())
            new_values.resize(values.size());

        double sum_value = 0.0;
        for (std::size_t i = 0; i < values.size(); ++i)
            sum_value += mWeights[i] * values[i];
        for (std::size_t i = 0; i < new_values.size(); ++i)
            new_values[i] = mWeights[i]*values[i] / sum_value;
    }

    /// Get the derivatives of the basis function i at point xi
    virtual void GetDerivative(std::vector<double>& new_dvalues, const std::size_t& i, const std::vector<double>& xi) const
    {
        std::vector<double> values;
        std::vector<std::vector<double> > dvalues;
        mpFESpace->GetValueAndDerivative(values, dvalues, xi);

        double sum_value = 0.0;
        std::vector<double> dsum_value(TDim);
        std::fill(dsum_value.begin(), dsum_value.end(), 0.0);
        for (std::size_t j = 0; j < values.size(); ++j)
        {
            sum_value += mWeights[j] * values[j];
            for (int dim = 0; dim < TDim; ++dim)
                dsum_value[dim] += mWeights[j] * dvalues[j][dim];
        }

        if (new_dvalues.size() != TDim)
            new_dvalues.resize(TDim);

        for (int dim = 0; dim < TDim; ++dim)
        {
            new_dvalues[dim] = mWeights[i] * (dvalues[i][dim]/sum_value - values[i]*dsum_value[dim]/pow(sum_value, 2));
        }
    }

    /// Get the derivatives of the basis functions at point xi
    /// the output values has the form of values[func_index][dim_index]
    virtual void GetDerivative(std::vector<std::vector<double> >& new_dvalues, const std::vector<double>& xi) const
    {
        std::vector<double> values;
        std::vector<std::vector<double> > dvalues;
        mpFESpace->GetValueAndDerivative(values, dvalues, xi);

        double sum_value = 0.0;
        std::vector<double> dsum_value(TDim);
        std::fill(dsum_value.begin(), dsum_value.end(), 0.0);
        for (std::size_t i = 0; i < values.size(); ++i)
        {
            sum_value += mWeights[i] * values[i];
            for (int dim = 0; dim < TDim; ++dim)
                dsum_value[dim] += mWeights[i] * dvalues[i][dim];
        }

        // KRATOS_WATCH(sum_value)
        // KRATOS_WATCH(dsum_value[0])

        if (new_dvalues.size() != dvalues.size())
            new_dvalues.resize(dvalues.size());
        for (std::size_t i = 0; i < new_dvalues.size(); ++i)
            if (new_dvalues[i].size() != TDim)
                new_dvalues[i].resize(TDim);

        for (std::size_t i = 0; i < new_dvalues.size(); ++i)
        {
            for (int dim = 0; dim < TDim; ++dim)
            {
                new_dvalues[i][dim] = mWeights[i] * (dvalues[i][dim]/sum_value - values[i]*dsum_value[dim]/pow(sum_value, 2));
            }
        }

        // for (std::size_t i = 0; i < new_dvalues.size(); ++i)
        //     KRATOS_WATCH(new_dvalues[i][0])
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Reset all the dof numbers for each grid function to -1.
    virtual void ResetFunctionIndices()
    {
        mpFESpace->ResetFunctionIndices();
    }

    /// Reset the function indices to a given values.
    /// This is useful when assigning the id for the boundary patch.
    virtual void ResetFunctionIndices(const std::vector<std::size_t>& func_indices)
    {
        mpFESpace->ResetFunctionIndices(func_indices);
    }

    /// Enumerate the dofs of each grid function. The enumeration algorithm is pretty straightforward.
    /// If the dof does not have pre-existing value, which assume it is -1, it will be assigned the incremental value.
    virtual std::size_t& Enumerate(std::size_t& start)
    {
        return mpFESpace->Enumerate(start);
    }

    /// Access the function indices (aka global ids)
    virtual std::vector<std::size_t> FunctionIndices() const
    {
        return mpFESpace->FunctionIndices();
    }

    /// Update the function indices using a map. The map shall be the mapping from old index to new index.
    virtual void UpdateFunctionIndices(const std::map<std::size_t, std::size_t>& indices_map)
    {
        mpFESpace->UpdateFunctionIndices(indices_map);
    }

    /// Get the first equation_id in this space
    virtual std::size_t GetFirstEquationId() const
    {
        return mpFESpace->GetFirstEquationId();
    }

    /// Get the last equation_id in this space
    virtual std::size_t GetLastEquationId() const
    {
        return mpFESpace->GetLastEquationId();
    }

    /// Check the compatibility between boundaries of two WeightedFESpacees
    virtual bool CheckBoundaryCompatibility(const FESpace<TDim>& rFESpace1, const BoundarySide& side1,
            const FESpace<TDim>& rFESpace2, const BoundarySide& side2) const
    {
        return rFESpace1 == rFESpace2;
    }

    /// Validate the WeightedFESpace before using
    virtual bool Validate() const
    {
        if (mWeights.size() != this->TotalNumber())
        {
            KRATOS_THROW_ERROR(std::logic_error, "The weight information is incorrect", "")
            return false;
        }
        return mpFESpace->Validate();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Extract the index of the functions on the boundary
    virtual std::vector<std::size_t> ExtractBoundaryFunctionIndices(const BoundarySide& side) const
    {
        return mpFESpace->ExtractBoundaryFunctionIndices(side);
    }

    /// Assign the index for the functions on the boundary
    virtual void AssignBoundaryFunctionIndices(const BoundarySide& side, const std::vector<std::size_t>& func_indices)
    {
        mpFESpace->AssignBoundaryFunctionIndices(side, func_indices);
    }

    /// Construct the boundary WeightedFESpace based on side
    virtual typename FESpace<TDim-1>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    {
        typename FESpace<TDim-1>::Pointer pBFESpace = mpFESpace->ConstructBoundaryFESpace(side);
        // TODO extract/compute the weights on the boundary
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not completed")
        std::vector<double> boundary_weights;

        return typename WeightedFESpace<TDim-1>::Pointer(new WeightedFESpace<TDim-1>(pBFESpace, boundary_weights));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compare the two weighted FESpaces in terms of its parametric information.
    virtual bool IsCompatible(const FESpace<TDim>& rOtherFESpace) const
    {
        if (rOtherFESpace.Type() == this->Type())
        {
            const WeightedFESpace<TDim>& rOtherWeightedFESpace = dynamic_cast<const WeightedFESpace<TDim>&>(rOtherFESpace);
            if (this->Weights().size() != rOtherWeightedFESpace.Weights().size())
            {
                return false;
            }
            else
            {
                for (std::size_t i = 0; i < this->Weights().size(); ++i)
                    if (this->Weights()[i] != rOtherWeightedFESpace.Weights()[i])
                        return false;
                return rOtherFESpace.IsCompatible(static_cast<const FESpace<TDim>&>(*this));
            }
        }
        else
        {
            for (std::size_t i = 0; i < this->Weights().size(); ++i)
                if (this->Weights()[i] != 1.0)
                    return false;
            return rOtherFESpace.IsCompatible(static_cast<const FESpace<TDim>&>(*this));
        }
        return false;
    }

    /// Overload comparison operator
    virtual bool operator==(const WeightedFESpace<TDim>& rOther) const
    {
        return this->IsCompatible(rOther);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Create the cell manager for all the cells in the support domain of the WeightedFESpace
    virtual typename cell_container_t::Pointer ConstructCellManager() const
    {
        typename cell_container_t::Pointer pCellManager = mpFESpace->ConstructCellManager();
        // TODO add weight information to the cells
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not completed")
        return pCellManager;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Overload assignment operator
    WeightedFESpace<TDim>& operator=(const WeightedFESpace<TDim>& rOther)
    {
        BaseType::operator=(rOther);
        this->mpFESpace = rOther.mpFESpace;
        this->mWeights = rOther.mWeights;
        return *this;
    }

    /// Clone this WeightedFESpace, this is a deep copy operation
    virtual typename FESpace<TDim>::Pointer Clone() const
    {
        typename WeightedFESpace<TDim>::Pointer pNewWeightedFESpace = typename WeightedFESpace<TDim>::Pointer(new WeightedFESpace<TDim>(mpFESpace, mWeights));
        *pNewWeightedFESpace = *this;
        return pNewWeightedFESpace;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        BaseType::PrintInfo(rOStream);
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        BaseType::PrintData(rOStream);
        rOStream << " Weights:";
        for (std::size_t i = 0; i < mWeights.size(); ++i)
            rOStream << " " << mWeights[i];
    }

private:

    typename BaseType::Pointer mpFESpace;
    std::vector<double> mWeights;

    /// Serializer
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer,  BaseType );
        rSerializer.save( "mWeights", mWeights );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer,  BaseType );
        rSerializer.load( "mWeights", mWeights );
    }
};

/**
 * Template specific instantiation for null-D WeightedFESpace to terminate the compilation
 */
template<>
class WeightedFESpace<0> : public FESpace<0>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(WeightedFESpace);

    // Type definitions
    typedef CellManager<Cell> cell_container_t;

    /// Default constructor
    WeightedFESpace(FESpace<0>::Pointer pFESpace, const std::vector<double>& weights) : mFunctionId(-1) {}

    /// Destructor
    virtual ~WeightedFESpace() {}

    /// Get the number of basis functions defined over the WeightedFESpace
    virtual const std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the WeightedFESpace in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
    {
        return 0;
    }

    /// Get the string describing the type of the WeightedFESpace
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the WeightedFESpace
    static std::string StaticType()
    {
        return "FESpace0D";
    }

    /// Overload comparison operator
    virtual bool operator==(const WeightedFESpace<0>& rOther) const
    {
        return true;
    }

    /// Reset all the dof numbers for each grid function to -1
    void ResetFunctionIndices()
    {
        mFunctionId = -1;
    }

    /// Reset the function indices to a given values
    void ResetFunctionIndices(const std::vector<std::size_t>& func_indices)
    {
        assert(func_indices.size() == 1);
        mFunctionId = func_indices[0];
    }

    /// Get the vector of function indices
    std::vector<std::size_t> FunctionIndices() const {return std::vector<std::size_t>{mFunctionId};}

    /// Check the compatibility between boundaries of two FESpacees
    virtual bool CheckBoundaryCompatibility(const WeightedFESpace<0>& rFESpace1, const BoundarySide& side1,
            const WeightedFESpace<0>& rFESpace2, const BoundarySide& side2) const
    {
        return true;
    }

    /// Validate the WeightedFESpace before using
    virtual bool Validate() const
    {
        return true;
    }

    // /// Construct the boundary WeightedFESpace based on side
    // virtual typename WeightedFESpace<-1>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    // {
    //     return NULL;
    // }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "WeightedFESpace<0>";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

protected:

    std::size_t mFunctionId;

private:

    /// Serializer
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
    }

    virtual void load(Serializer& rSerializer)
    {
    }
};


/**
 * Template specific instantiation for -1-D WeightedFESpace to terminate the compilation
 */
template<>
class WeightedFESpace<-1> : public FESpace<-1>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(WeightedFESpace);

    /// Default constructor
    WeightedFESpace(FESpace<-1>::Pointer pFESpace, const std::vector<double>& weights) {}

    /// Destructor
    virtual ~WeightedFESpace() {}

    /// Get the number of basis functions defined over the WeightedFESpace
    virtual const std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the WeightedFESpace in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
    {
        return 0;
    }

    /// Get the string describing the type of the WeightedFESpace
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the WeightedFESpace
    static std::string StaticType()
    {
        return "WeightedFESpace<-1>D";
    }

    /// Overload comparison operator
    virtual bool operator==(const WeightedFESpace<-1>& rOther) const
    {
        return true;
    }

    /// Check the compatibility between boundaries of two FESpacees
    virtual bool CheckBoundaryCompatibility(const WeightedFESpace<-1>& rFESpace1, const BoundarySide& side1,
            const WeightedFESpace<-1>& rFESpace2, const BoundarySide& side2) const
    {
        return true;
    }

    /// Validate the WeightedFESpace before using
    virtual bool Validate() const
    {
        return true;
    }

    // /// Construct the boundary WeightedFESpace based on side
    // virtual typename WeightedFESpace<-2>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    // {
    //     return NULL;
    // }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "WeightedFESpace<-1>";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};


/**
 * Template specific instantiation for -2-D WeightedFESpace to terminate the compilation
 */
template<>
class WeightedFESpace<-2> : public FESpace<-2>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(WeightedFESpace);

    /// Default constructor
    WeightedFESpace(FESpace<-2>::Pointer pFESpace, const std::vector<double>& weights) {}

    /// Destructor
    virtual ~WeightedFESpace() {}

    /// Get the number of basis functions defined over the WeightedFESpace
    virtual const std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the WeightedFESpace in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
    {
        return 0;
    }

    /// Get the string describing the type of the WeightedFESpace
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the WeightedFESpace
    static std::string StaticType()
    {
        return "WeightedFESpace<-2>D";
    }

    /// Overload comparison operator
    virtual bool operator==(const WeightedFESpace<-2>& rOther) const
    {
        return true;
    }

    /// Check the compatibility between boundaries of two FESpacees
    virtual bool CheckBoundaryCompatibility(const WeightedFESpace<-2>& rFESpace1, const BoundarySide& side1,
            const WeightedFESpace<-2>& rFESpace2, const BoundarySide& side2) const
    {
        return true;
    }

    /// Validate the WeightedFESpace before using
    virtual bool Validate() const
    {
        return true;
    }

    // /// Construct the boundary WeightedFESpace based on side
    // virtual typename WeightedFESpace<-3>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    // {
    //     return NULL;
    // }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "WeightedFESpace<-2>";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};

/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const WeightedFESpace<TDim>& rThis)
{
    rOStream << "-------------Begin WeightedFESpaceInfo-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "-------------End WeightedFESpaceInfo-------------";
    return rOStream;
}

} // namespace Kratos.

#undef DEBUG_DESTROY

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_WEIGHTED_FESPACE_H_INCLUDED defined

