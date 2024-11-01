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
#ifdef ISOGEOMETRIC_DEBUG_DESTROY
        std::cout << this->Type() << ", Addr = " << this << " is destroyed" << std::endl;
#endif
    }

    /// Helper to create new WeightedFESpace pointer
    static typename WeightedFESpace<TDim>::Pointer Create(typename BaseType::Pointer pFESpace, const std::vector<double>& weights)
    {
        return typename WeightedFESpace<TDim>::Pointer(new WeightedFESpace(pFESpace, weights));
    }

    /// Get the number of basis functions defined over the WeightedFESpace
    std::size_t TotalNumber() const override
    {
        return mpFESpace->TotalNumber();
    }

    /// Get the order of the WeightedFESpace in specific direction
    std::size_t Order(std::size_t i) const override
    {
        return mpFESpace->Order(i);
    }

    /// Set the weight vector
    void SetWeights(const std::vector<double>& weights)
    {
        if (mWeights.size() != weights.size())
        {
            mWeights.resize(weights.size());
        }
        std::copy(weights.begin(), weights.end(), mWeights.begin());
    }

    /// Get the weight vector
    const std::vector<double>& Weights() const {return mWeights;}

    /// Get the string representing the type of the WeightedFESpace
    std::string Type() const override
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
    void GetValue(double& v, std::size_t i, const std::vector<double>& xi) const override
    {
        std::vector<double> values;
        mpFESpace->GetValues(values, xi);

        double sum_value = 0.0;
        for (std::size_t j = 0; j < values.size(); ++j)
        {
            sum_value += mWeights[j] * values[j];
        }

        if (sum_value == 0.0)
        {
            v = 0.0;
        }
        else
        {
            v = mWeights[i] * values[i] / sum_value;
        }
    }

    /// Get the values of the basis functions at point xi
    void GetValues(std::vector<double>& new_values, const std::vector<double>& xi) const override
    {
        std::vector<double> values;
        mpFESpace->GetValues(values, xi);

        if (new_values.size() != values.size())
        {
            new_values.resize(values.size());
        }

        double sum_value = 0.0;
        for (std::size_t i = 0; i < values.size(); ++i)
        {
            sum_value += mWeights[i] * values[i];
        }

        if (sum_value == 0.0)
        {
            std::fill(new_values.begin(), new_values.end(), 0.0);
        }
        else
            for (std::size_t i = 0; i < new_values.size(); ++i)
            {
                new_values[i] = mWeights[i] * values[i] / sum_value;
            }
    }

    /// Get the derivatives of the basis function i at point xi
    void GetDerivative(std::vector<double>& new_dvalues, std::size_t i, const std::vector<double>& xi) const override
    {
        std::vector<double> values;
        std::vector<std::vector<double> > dvalues;
        mpFESpace->GetValuesAndDerivatives(values, dvalues, xi);

        if (new_dvalues.size() != TDim)
        {
            new_dvalues.resize(TDim);
        }

        double sum_value = 0.0;
        std::vector<double> dsum_value(TDim);
        std::fill(dsum_value.begin(), dsum_value.end(), 0.0);
        for (std::size_t j = 0; j < values.size(); ++j)
        {
            sum_value += mWeights[j] * values[j];
            for (int dim = 0; dim < TDim; ++dim)
            {
                dsum_value[dim] += mWeights[j] * dvalues[j][dim];
            }
        }

        if (sum_value == 0.0)
        {
            std::fill(new_dvalues.begin(), new_dvalues.end(), 0.0);
        }
        else
            for (int dim = 0; dim < TDim; ++dim)
            {
                new_dvalues[dim] = mWeights[i] * (dvalues[i][dim] / sum_value - values[i] * dsum_value[dim] / pow(sum_value, 2));
            }
    }

    /// Get the derivatives of the basis functions at point xi
    /// the output values has the form of values[func_index][dim_index]
    void GetDerivatives(std::vector<std::vector<double> >& new_dvalues, const std::vector<double>& xi) const override
    {
        std::vector<double> values;
        std::vector<std::vector<double> > dvalues;
        mpFESpace->GetValuesAndDerivatives(values, dvalues, xi);

        if (new_dvalues.size() != dvalues.size())
        {
            new_dvalues.resize(dvalues.size());
        }
        for (std::size_t i = 0; i < new_dvalues.size(); ++i)
        {
            if (new_dvalues[i].size() != TDim)
            {
                new_dvalues[i].resize(TDim);
            }
        }

        double sum_value = 0.0;
        std::vector<double> dsum_value(TDim);
        std::fill(dsum_value.begin(), dsum_value.end(), 0.0);
        for (std::size_t i = 0; i < values.size(); ++i)
        {
            sum_value += mWeights[i] * values[i];
            for (int dim = 0; dim < TDim; ++dim)
            {
                dsum_value[dim] += mWeights[i] * dvalues[i][dim];
            }
        }

        if (sum_value == 0.0)
        {
            for (std::size_t i = 0; i < new_dvalues.size(); ++i)
            {
                std::fill(new_dvalues[i].begin(), new_dvalues[i].end(), 0.0);
            }
        }
        else
        {
            for (std::size_t i = 0; i < new_dvalues.size(); ++i)
            {
                for (int dim = 0; dim < TDim; ++dim)
                {
                    new_dvalues[i][dim] = mWeights[i] * (dvalues[i][dim] / sum_value - values[i] * dsum_value[dim] / pow(sum_value, 2));
                }
            }
        }
    }

    /// Get the derivatives of the basis functions at point xi
    /// the output values has the form of values[func_index][dim_index]
    void GetDerivatives(const unsigned int nd, std::vector<std::vector<std::vector<double> > >& new_dvalues, const std::vector<double>& xi) const override
    {
        std::vector<double> values;
        std::vector<std::vector<std::vector<double> > > dvalues;
        mpFESpace->GetValuesAndDerivatives(nd, values, dvalues, xi);

        if (new_dvalues.size() != nd)
            new_dvalues.resize(nd);

        for (unsigned int i = 0; i < nd; ++i)
        {
            if (new_dvalues[i].size() != dvalues[i].size())
            {
                new_dvalues[i].resize(dvalues[i].size());
            }
            for (std::size_t j = 0; j < new_dvalues[i].size(); ++j)
            {
                if (new_dvalues[i][j].size() != dvalues[i][j].size())
                {
                    new_dvalues[i][j].resize(dvalues[i][j].size());
                }
            }
        }

        double sum_value = 0.0;
        for (std::size_t i = 0; i < values.size(); ++i)
        {
            sum_value += mWeights[i] * values[i];
        }

        std::vector<double> dsum_value;
        if (nd > 0)
        {
            dsum_value.resize(TDim);
            std::fill(dsum_value.begin(), dsum_value.end(), 0.0);
            for (std::size_t i = 0; i < values.size(); ++i)
            {
                for (int dim = 0; dim < TDim; ++dim)
                {
                    dsum_value[dim] += mWeights[i] * dvalues[0][i][dim];
                }
            }
        }

        std::vector<double> d2sum_value;
        if (nd > 0)
        {
            const unsigned int nders = TDim*(TDim+1)/2;
            d2sum_value.resize(nders);
            std::fill(d2sum_value.begin(), d2sum_value.end(), 0.0);
            for (std::size_t i = 0; i < values.size(); ++i)
            {
                for (int dim = 0; dim < nders; ++dim)
                {
                    d2sum_value[dim] += mWeights[i] * dvalues[1][i][dim];
                }
            }
        }

        if (sum_value == 0.0)
        {
            for (unsigned int i = 0; i < nd; ++i)
            {
                for (std::size_t j = 0; j < new_dvalues[i].size(); ++j)
                {
                    std::fill(new_dvalues[i][j].begin(), new_dvalues[i][j].end(), 0.0);
                }
            }
        }
        else
        {
            if (nd > 0)
            {
                for (std::size_t i = 0; i < new_dvalues[0].size(); ++i)
                {
                    for (int dim = 0; dim < TDim; ++dim)
                    {
                        new_dvalues[0][i][dim] = mWeights[i] * (
                                dvalues[0][i][dim] / sum_value
                              - values[i] * dsum_value[dim] / pow(sum_value, 2)
                            );
                    }
                }
            }

            if (nd > 1)
            {
                for (std::size_t i = 0; i < new_dvalues[1].size(); ++i)
                {
                    for (int dim = 0; dim < TDim; ++dim)
                    {
                        new_dvalues[1][i][dim] = mWeights[i] * (
                                dvalues[1][i][dim] / sum_value
                              - 2 * dvalues[0][i][dim] * dsum_value[dim] / pow(sum_value, 2)
                              - values[i] * d2sum_value[dim] / pow(sum_value, 2)
                              + 2 * values[i] * pow(dsum_value[dim], 2) / pow(sum_value, 3)
                            );
                    }

                    if constexpr (TDim == 2)
                    {
                        new_dvalues[1][i][2] = mWeights[i] * (
                                dvalues[1][i][2] / sum_value
                              - dvalues[0][i][0] * dsum_value[1] / pow(sum_value, 2)
                              - dvalues[0][i][1] * dsum_value[0] / pow(sum_value, 2)
                              - values[i] * d2sum_value[2] / pow(sum_value, 2)
                              + 2 * values[i] * dsum_value[0] * dsum_value[1] / pow(sum_value, 3)
                            );
                    }
                    else if constexpr (TDim == 3)
                    {
                        new_dvalues[1][i][3] = mWeights[i] * (
                                dvalues[1][i][3] / sum_value
                              - dvalues[0][i][0] * dsum_value[1] / pow(sum_value, 2)
                              - dvalues[0][i][1] * dsum_value[0] / pow(sum_value, 2)
                              - values[i] * d2sum_value[3] / pow(sum_value, 2)
                              + 2 * values[i] * dsum_value[0] * dsum_value[1] / pow(sum_value, 3)
                            );

                        new_dvalues[1][i][4] = mWeights[i] * (
                                dvalues[1][i][4] / sum_value
                              - dvalues[0][i][1] * dsum_value[2] / pow(sum_value, 2)
                              - dvalues[0][i][2] * dsum_value[1] / pow(sum_value, 2)
                              - values[i] * d2sum_value[4] / pow(sum_value, 2)
                              + 2 * values[i] * dsum_value[1] * dsum_value[2] / pow(sum_value, 3)
                            );

                        new_dvalues[1][i][5] = mWeights[i] * (
                                dvalues[1][i][5] / sum_value
                              - dvalues[0][i][0] * dsum_value[2] / pow(sum_value, 2)
                              - dvalues[0][i][2] * dsum_value[0] / pow(sum_value, 2)
                              - values[i] * d2sum_value[5] / pow(sum_value, 2)
                              + 2 * values[i] * dsum_value[0] * dsum_value[2] / pow(sum_value, 3)
                            );
                    }
                }
            }

            if (nd > 2)
            {
                KRATOS_ERROR << "Higher order " << nd << " > 2 is not yet supported";
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Check if a point lies inside the parametric domain of the BSplinesFESpace
    bool IsInside(const std::vector<double>& xi) const override
    {
        return mpFESpace->IsInside(xi);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Reset all the dof numbers for each grid function to -1.
    void ResetFunctionIndices() override
    {
        mpFESpace->ResetFunctionIndices();
    }

    /// Reset the function indices to a given values.
    /// This is useful when assigning the id for the boundary patch.
    void ResetFunctionIndices(const std::vector<std::size_t>& func_indices) override
    {
        mpFESpace->ResetFunctionIndices(func_indices);
    }

    /// Enumerate the dofs of each grid function. The enumeration algorithm is pretty straightforward.
    /// If the dof does not have pre-existing value, which assume it is -1, it will be assigned the incremental value.
    std::size_t& Enumerate(std::size_t& start) override
    {
        return mpFESpace->Enumerate(start);
    }

    /// Access the function indices (aka global ids)
    std::vector<std::size_t> FunctionIndices() const override
    {
        return mpFESpace->FunctionIndices();
    }

    /// Update the function indices using a map. The map shall be the mapping from old index to new index.
    void UpdateFunctionIndices(const std::map<std::size_t, std::size_t>& indices_map) override
    {
        mpFESpace->UpdateFunctionIndices(indices_map);
    }

    /// Get the first equation_id in this space
    std::size_t GetFirstEquationId() const override
    {
        return mpFESpace->GetFirstEquationId();
    }

    /// Get the last equation_id in this space
    std::size_t GetLastEquationId() const override
    {
        return mpFESpace->GetLastEquationId();
    }

    /// Check the compatibility between boundaries of two WeightedFESpacees
    bool CheckBoundaryCompatibility(const FESpace<TDim>& rFESpace1, const BoundarySide& side1,
                                    const FESpace<TDim>& rFESpace2, const BoundarySide& side2) const override
    {
        return rFESpace1 == rFESpace2;
    }

    /// Validate the WeightedFESpace before using
    bool Validate() const override
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
    std::vector<std::size_t> ExtractBoundaryFunctionIndices(const BoundarySide& side) const override
    {
        return mpFESpace->ExtractBoundaryFunctionIndices(side);
    }

    /// Assign the index for the functions on the boundary
    void AssignBoundaryFunctionIndices(const BoundarySide& side, const std::vector<std::size_t>& func_indices) override
    {
        mpFESpace->AssignBoundaryFunctionIndices(side, func_indices);
    }

    /// Construct the boundary WeightedFESpace based on side
    typename FESpace < TDim - 1 >::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const override
    {
        typename FESpace < TDim - 1 >::Pointer pBFESpace = mpFESpace->ConstructBoundaryFESpace(side);
        // TODO extract/compute the weights on the boundary
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not completed")
        std::vector<double> boundary_weights;

        return typename WeightedFESpace < TDim - 1 >::Pointer(new WeightedFESpace < TDim - 1 > (pBFESpace, boundary_weights));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compare the two weighted FESpaces in terms of its parametric information.
    bool IsCompatible(const FESpace<TDim>& rOtherFESpace) const override
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
                    {
                        return false;
                    }
                return rOtherFESpace.IsCompatible(static_cast<const FESpace<TDim>&>(*this));
            }
        }
        else
        {
            for (std::size_t i = 0; i < this->Weights().size(); ++i)
                if (this->Weights()[i] != 1.0)
                {
                    return false;
                }
            return rOtherFESpace.IsCompatible(static_cast<const FESpace<TDim>&>(*this));
        }
        return false;
    }

    /// Overload comparison operator
    bool operator==(const FESpace<TDim>& rOther) const override
    {
        return this->IsCompatible(rOther);
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
    typename FESpace<TDim>::Pointer Clone() const override
    {
        typename WeightedFESpace<TDim>::Pointer pNewWeightedFESpace = typename WeightedFESpace<TDim>::Pointer(new WeightedFESpace<TDim>(mpFESpace, mWeights));
        *pNewWeightedFESpace = *this;
        return pNewWeightedFESpace;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        BaseType::PrintInfo(rOStream);
    }

    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
        rOStream << " Weights:";
        for (std::size_t i = 0; i < mWeights.size(); ++i)
        {
            rOStream << " " << mWeights[i];
        }
    }

private:

    typename BaseType::Pointer mpFESpace;
    std::vector<double> mWeights;

    /// Serializer
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer,  BaseType );
        rSerializer.save( "mWeights", mWeights );
    }

    void load(Serializer& rSerializer) override
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

    /// Default constructor
    WeightedFESpace(FESpace<0>::Pointer pFESpace, const std::vector<double>& weights) : mFunctionId(-1) {}

    /// Destructor
    virtual ~WeightedFESpace() {}

    /// Get the number of basis functions defined over the WeightedFESpace
    std::size_t TotalNumber() const override
    {
        return 0;
    }

    /// Get the order of the WeightedFESpace in specific direction
    std::size_t Order(std::size_t i) const override
    {
        return 0;
    }

    /// Get the string describing the type of the WeightedFESpace
    std::string Type() const override
    {
        return StaticType();
    }

    /// Get the string describing the type of the WeightedFESpace
    static std::string StaticType()
    {
        return "FESpace0D";
    }

    /// Overload comparison operator
    bool operator==(const FESpace<0>& rOther) const override
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
    std::vector<std::size_t> FunctionIndices() const {return std::vector<std::size_t> {mFunctionId};}

    /// Check the compatibility between boundaries of two FESpacees
    bool CheckBoundaryCompatibility(const FESpace<0>& rFESpace1, const BoundarySide& side1,
                                    const FESpace<0>& rFESpace2, const BoundarySide& side2) const override
    {
        return true;
    }

    /// Validate the WeightedFESpace before using
    bool Validate() const override
    {
        return true;
    }

    // /// Construct the boundary WeightedFESpace based on side
    // virtual typename WeightedFESpace<-1>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    // {
    //     return NULL;
    // }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "WeightedFESpace<0>";
    }

    void PrintData(std::ostream& rOStream) const override
    {
    }

protected:

    std::size_t mFunctionId;

private:

    /// Serializer
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
    }

    void load(Serializer& rSerializer) override
    {
    }
};

/**
 * Template specific instantiation for -1-D WeightedFESpace to terminate the compilation
 */
template<>
class WeightedFESpace < -1 > : public FESpace < -1 >
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(WeightedFESpace);

    /// Default constructor
    WeightedFESpace(FESpace < -1 >::Pointer pFESpace, const std::vector<double>& weights) {}

    /// Destructor
    virtual ~WeightedFESpace() {}

    /// Get the number of basis functions defined over the WeightedFESpace
    std::size_t TotalNumber() const override
    {
        return 0;
    }

    /// Get the order of the WeightedFESpace in specific direction
    std::size_t Order(std::size_t i) const override
    {
        return 0;
    }

    /// Get the string describing the type of the WeightedFESpace
    std::string Type() const override
    {
        return StaticType();
    }

    /// Get the string describing the type of the WeightedFESpace
    static std::string StaticType()
    {
        return "WeightedFESpace<-1>D";
    }

    /// Overload comparison operator
    bool operator==(const FESpace < -1 > & rOther) const override
    {
        return true;
    }

    /// Check the compatibility between boundaries of two FESpacees
    bool CheckBoundaryCompatibility(const FESpace < -1 > & rFESpace1, const BoundarySide& side1,
                                    const FESpace < -1 > & rFESpace2, const BoundarySide& side2) const override
    {
        return true;
    }

    /// Validate the WeightedFESpace before using
    bool Validate() const override
    {
        return true;
    }

    // /// Construct the boundary WeightedFESpace based on side
    // virtual typename WeightedFESpace<-2>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    // {
    //     return NULL;
    // }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "WeightedFESpace<-1>";
    }

    void PrintData(std::ostream& rOStream) const override
    {
    }
};

/**
 * Template specific instantiation for -2-D WeightedFESpace to terminate the compilation
 */
template<>
class WeightedFESpace < -2 > : public FESpace < -2 >
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(WeightedFESpace);

    /// Default constructor
    WeightedFESpace(FESpace < -2 >::Pointer pFESpace, const std::vector<double>& weights) {}

    /// Destructor
    virtual ~WeightedFESpace() {}

    /// Get the number of basis functions defined over the WeightedFESpace
    std::size_t TotalNumber() const override
    {
        return 0;
    }

    /// Get the order of the WeightedFESpace in specific direction
    std::size_t Order(std::size_t i) const override
    {
        return 0;
    }

    /// Get the string describing the type of the WeightedFESpace
    std::string Type() const override
    {
        return StaticType();
    }

    /// Get the string describing the type of the WeightedFESpace
    static std::string StaticType()
    {
        return "WeightedFESpace<-2>D";
    }

    /// Overload comparison operator
    bool operator==(const FESpace < -2 > & rOther) const override
    {
        return true;
    }

    /// Check the compatibility between boundaries of two FESpacees
    bool CheckBoundaryCompatibility(const FESpace < -2 > & rFESpace1, const BoundarySide& side1,
                                    const FESpace < -2 > & rFESpace2, const BoundarySide& side2) const override
    {
        return true;
    }

    /// Validate the WeightedFESpace before using
    bool Validate() const override
    {
        return true;
    }

    // /// Construct the boundary WeightedFESpace based on side
    // virtual typename WeightedFESpace<-3>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    // {
    //     return NULL;
    // }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "WeightedFESpace<-2>";
    }

    void PrintData(std::ostream& rOStream) const override
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

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_WEIGHTED_FESPACE_H_INCLUDED defined
