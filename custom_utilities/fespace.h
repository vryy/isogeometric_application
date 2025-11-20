//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Nov 2017 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_FESPACE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_FESPACE_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "custom_utilities/iga_define.h"
#include "custom_utilities/cell_container.h"

namespace Kratos
{

/**
 * An FESpace is a collection of shape function defined over the parametric domain.
 * An isogeometric FESpace can be a NURBS FESpace, a hierarchical NURBS FESpace, or a T-Splines FESpace.
 */
template<int TDim, typename TLocalCoordinateType = double>
class FESpace
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(FESpace);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const FESpace> ConstPointer;
#endif

    /// Type definition
    typedef CellContainer cell_container_t;
    typedef TLocalCoordinateType LocalCoordinateType;

    typedef FESpace<TDim, TLocalCoordinateType> FESpaceType;
    typedef FESpace<TDim-1, TLocalCoordinateType> BoundaryFESpaceType;

    /// Default constructor
    FESpace() {}

    /// Copy constructor
    FESpace(const FESpace& rOther)
    : mGlobalToLocal(rOther.mGlobalToLocal)
    {}

    /// Destructor
    virtual ~FESpace()
    {
#ifdef ISOGEOMETRIC_DEBUG_DESTROY
        this->PrintInfo(std::cout); std::cout << " is destroyed" << std::endl;
#endif
    }

    /// Helper to create new BSplinesFESpace pointer
    static typename FESpaceType::Pointer Create()
    {
        return typename FESpaceType::Pointer(new FESpace());
    }

    /// Get the dimension of the FESpace
    static constexpr int Dim()
    {
        return TDim;
    }

    /// Get the number of basis functions defined over the FESpace
    virtual std::size_t TotalNumber() const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Get the order of the FESpace in specific direction
    virtual std::size_t Order(std::size_t di) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Get the lower and upper bound of the parametric space in a specific direction
    virtual std::vector<LocalCoordinateType> ParametricBounds(std::size_t di) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Get the string representing the type of the FESpace
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string representing the type of the FESpace
    static std::string StaticType()
    {
        std::stringstream ss;
        ss << "FESpace" << TDim << "D";
        return ss.str();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Get the value of the basis function i at point xi
    virtual void GetValue(double& v, std::size_t i, const std::vector<LocalCoordinateType>& xi) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Get the value of the basis function i at point xi
    double GetValue(std::size_t i, const std::vector<LocalCoordinateType>& xi) const
    {
        double v;
        this->GetValue(v, i, xi);
        return v;
    }

    /// Get the values of the basis functions at point xi
    virtual void GetValues(std::vector<double>& values, const std::vector<LocalCoordinateType>& xi) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Get the values of the basis functions at point xi
    std::vector<double> GetValues(const std::vector<LocalCoordinateType>& xi) const
    {
        std::vector<double> values;
        this->GetValues(values, xi);
        return values;
    }

    ///////////////

    /// Get the derivative of the basis function i at point xi
    virtual void GetDerivative(std::vector<double>& values, std::size_t i, const std::vector<LocalCoordinateType>& xi) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Get the derivatives of the basis function i at point xi
    std::vector<double> GetDerivatives(std::size_t i, const std::vector<LocalCoordinateType>& xi) const
    {
        std::vector<double> values;
        this->GetDerivatives(values, i, xi);
        return values;
    }

    /// Get the derivatives of the basis functions at point xi
    /// the output values has the form of values[func_index][dim_index]
    virtual void GetDerivatives(std::vector<std::vector<double> >& values, const std::vector<LocalCoordinateType>& xi) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Get the derivatives of the basis functions at point xi
    /// the return values has the form of values[func_index][dim_index]
    std::vector<std::vector<double> > GetDerivatives(const std::vector<LocalCoordinateType>& xi) const
    {
        std::vector<std::vector<double> > values;
        this->GetDerivatives(values);
        return values;
    }

    /// Get the derivatives upto nd of the basis functions at point xi
    /// the output values has the form of values[func_index][dim_index]
    virtual void GetDerivatives(const unsigned int nd, std::vector<std::vector<std::vector<double> > >& values, const std::vector<LocalCoordinateType>& xi) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    ///////////////

    /// Get the values and derivatives of the basis functions at point xi
    /// the output derivatives has the form of values[func_index][dim_index]
    virtual void GetValuesAndDerivatives(std::vector<double>& values,
            std::vector<std::vector<double> >& derivatives,
            const std::vector<LocalCoordinateType>& xi) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Get the values and derivatives upto nd (>=0) of the basis functions at point xi
    /// the output derivatives has the form of values[der_index][func_index][param_index]
    /// For the first derivatives, param_index is 0,1,..,dim
    /// For the second derivatives, param_index is 00,11,..,01,02,..,11,12,..
    /// For the higher derivatives, param_index is [d1][d2][d3] (d1<=d2<=d3)
    virtual void GetValuesAndDerivatives(const unsigned int nd, std::vector<double>& values,
            std::vector<std::vector<std::vector<double> > >& derivatives,
            const std::vector<LocalCoordinateType>& xi) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Check if a point lies inside the parametric domain of the FESpace
    virtual bool IsInside(const std::vector<LocalCoordinateType>& xi) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Reset all the dof numbers for each grid function to -1.
    virtual void ResetFunctionIndices()
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Reset the function indices to a given values.
    /// This is useful when assigning the id for the boundary patch.
    virtual void ResetFunctionIndices(const std::vector<std::size_t>& func_indices)
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Enumerate the dofs of each grid function. The enumeration algorithm is pretty straightforward.
    /// If the dof does not have pre-existing value, which assume it is -1, it will be assigned the incremental value.
    virtual std::size_t& Enumerate(std::size_t& start)
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Access the function indices (aka global ids)
    virtual std::vector<std::size_t> FunctionIndices() const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Update the function indices using a map. The map shall be the mapping from old index to new index.
    virtual void UpdateFunctionIndices(const std::map<std::size_t, std::size_t>& indices_map)
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Get the first equation_id in this space
    virtual std::size_t GetFirstEquationId() const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Get the last equation_id in this space
    virtual std::size_t GetLastEquationId() const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Return the local id of a given global id
    std::size_t LocalId(std::size_t global_id) const
    {
        std::map<std::size_t, std::size_t>::const_iterator it = mGlobalToLocal.find(global_id);

        if (it == mGlobalToLocal.end())
        {
            KRATOS_WATCH(TDim)
            std::cout << "mGlobalToLocal:";
            for (std::map<std::size_t, std::size_t>::const_iterator it2 = mGlobalToLocal.begin(); it2 != mGlobalToLocal.end(); ++it2)
            {
                std::cout << " " << it2->first << "->" << it2->second;
            }
            std::cout << std::endl;
            KRATOS_ERROR << "The global id " << global_id << " does not exist in mGlobalToLocal map";
        }

        return it->second;
    }

    /// Return the local ids of given global ids
    std::vector<std::size_t> LocalId(const std::vector<std::size_t>& global_ids) const
    {
        std::vector<std::size_t> local_ids(global_ids.size());
        for (std::size_t i = 0; i < global_ids.size(); ++i)
        {
            local_ids[i] = this->LocalId(global_ids[i]);
        }
        return local_ids;
    }

    /// Check the compatibility between boundaries of two FESpacees
    virtual bool CheckBoundaryCompatibility(const FESpaceType& rFESpace1, const BoundarySide& side1,
                                            const FESpaceType& rFESpace2, const BoundarySide& side2) const
    {
        typename BoundaryFESpaceType::Pointer pBFESpace1 = rFESpace1.ConstructBoundaryFESpace(side1);
        typename BoundaryFESpaceType::Pointer pBFESpace2 = rFESpace1.ConstructBoundaryFESpace(side2);

        return (*pBFESpace1) == (*pBFESpace2);
    }

    /// Validate the FESpace before using
    virtual bool Validate() const
    {
        return true;
    }

    /// Reverse the evaluation in i-direction
    virtual void Reverse(std::size_t idir)
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Transpose (swap) i- and j-direction
    virtual void Transpose(std::size_t idir, std::size_t jdir)
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Extract the index of the functions based on the boundary flag. It allows to extract the corner one.
    virtual std::vector<std::size_t> ExtractBoundaryFunctionIndicesByFlag(int boundary_id) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Extract the index of the functions on the boundary
    virtual std::vector<std::size_t> ExtractBoundaryFunctionIndices(std::vector<std::size_t>& size_info, const BoundarySide& side) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Extract the index of the functions on the boundary
    virtual std::vector<std::size_t> ExtractBoundaryFunctionIndices(const BoundarySide& side) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Extract the index of the functions on the boundary down to some level
    virtual std::vector<std::size_t> ExtractBoundaryFunctionIndices(const BoundarySide& side, std::size_t level) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Assign the index for the functions on the boundary
    /// If the override flag is on, all the boundary function indices will be overwritten if the target value is not -1
    /// If the override flag is off, only the boundary function index of not -1 is assigned
    virtual void AssignBoundaryFunctionIndices(const BoundarySide& side, const std::vector<std::size_t>& func_indices, const bool override)
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Construct the boundary FESpace based on side
    virtual typename BoundaryFESpaceType::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Construct the boundary FESpace based on side and local relative configuration
    virtual typename BoundaryFESpaceType::Pointer ConstructBoundaryFESpace(const BoundarySide& side,
            const std::map<std::size_t, std::size_t>& local_parameter_map, const std::vector<BoundaryDirection>& directions) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Construct the sliced FESpace
    virtual typename BoundaryFESpaceType::Pointer ConstructSlicedFESpace(int idir, double xi) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compare the two FESpace's in terms of its parametric information.
    virtual bool IsCompatible(const FESpaceType& rOtherFESpace) const
    {
        return false;
    }

    /// Comparison operator
    virtual bool operator==(const FESpaceType& rOther) const
    {
        const std::vector<std::size_t> my_func_indices = this->FunctionIndices();
        const std::vector<std::size_t> other_func_indices = rOther.FunctionIndices();
        return this->IsCompatible(rOther) && (my_func_indices == other_func_indices);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Fast function to get the opposite boundary side
    static BoundarySide OppositeBoundarySide(const BoundarySide& side)
    {
        if (side == _BLEFT_) { return _BRIGHT_; }
        else if (side == _BRIGHT_) { return _BLEFT_; }
        else if (side == _BTOP_) { return _BBOTTOM_; }
        else if (side == _BBOTTOM_) { return _BTOP_; }
        else if (side == _BFRONT_) { return _BBACK_; }
        else if (side == _BBACK_) { return _BFRONT_; }
        else
            KRATOS_ERROR << "Invalid boundary side " << side;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Create the cell manager for all the cells in the support domain of the FESpace
    virtual cell_container_t::Pointer ConstructCellManager() const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Overload assignment operator
    FESpaceType& operator=(const FESpaceType& rOther)
    {
        mGlobalToLocal = rOther.mGlobalToLocal;
        return *this;
    }

    /// Clone this FESpace, this is a deep copy operation
    virtual typename FESpaceType::Pointer Clone() const
    {
        typename FESpaceType::Pointer pNewFESpace = typename FESpaceType::Pointer(new FESpaceType());
        *pNewFESpace = *this;
        return pNewFESpace;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Type() << ", Add = " << this;
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << " Function Indices:";
        std::vector<std::size_t> func_indices = this->FunctionIndices();
        for (std::size_t i = 0; i < func_indices.size(); ++i)
        {
            rOStream << " " << func_indices[i];
        }
        rOStream << std::endl;
        rOStream << " GlobalToLocal:";
        for (std::map<std::size_t, std::size_t>::const_iterator it2 = mGlobalToLocal.begin(); it2 != mGlobalToLocal.end(); ++it2)
        {
            rOStream << " " << it2->first << "->" << it2->second;
        }
    }

protected:

    std::map<std::size_t, std::size_t> mGlobalToLocal;

private:

    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        rSerializer.save( "mGlobalToLocal", mGlobalToLocal );
    }

    virtual void load(Serializer& rSerializer)
    {
        rSerializer.load( "mGlobalToLocal", mGlobalToLocal );
    }
    ///@}
};

/**
 * Template specific instantiation for null-D FESpace to terminate the compilation
 */
template<typename TLocalCoordinateType>
class FESpace<0, TLocalCoordinateType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(FESpace);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const FESpace> ConstPointer;
#endif

    /// Type definition
    typedef CellContainer cell_container_t;

    typedef FESpace<0, TLocalCoordinateType> FESpaceType;

    /// Default constructor
    FESpace() : mFunctionId(-1) {}

    /// Destructor
    virtual ~FESpace() {}

    /// Helper to create new BSplinesFESpace pointer
    static typename FESpaceType::Pointer Create()
    {
        return typename FESpaceType::Pointer(new FESpace());
    }

    /// Get the number of basis functions defined over the FESpace
    virtual std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the FESpace in specific direction
    virtual std::size_t Order(std::size_t i) const
    {
        return 0;
    }

    /// Get the string describing the type of the FESpace
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the FESpace
    static std::string StaticType()
    {
        return "FESpace0D";
    }

    /// Overload comparison operator
    virtual bool operator==(const FESpaceType& rOther) const
    {
        return true;
    }

    /// Return the local id of a given global id
    std::size_t LocalId(std::size_t global_id) const
    {
        return 0;
    }

    /// Reset all the dof numbers for each grid function to -1
    virtual void ResetFunctionIndices()
    {
        mFunctionId = -1;
    }

    /// Reset the function indices to a given values
    virtual void ResetFunctionIndices(const std::vector<std::size_t>& func_indices)
    {
        assert(func_indices.size() == 1);
        mFunctionId = func_indices[0];
    }

    /// Update the function indices using a map. The map shall be the mapping from old index to new index.
    virtual void UpdateFunctionIndices(const std::map<std::size_t, std::size_t>& indices_map)
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Create the cell manager for all the cells in the support domain of the FESpace
    virtual cell_container_t::Pointer ConstructCellManager() const
    {
        KRATOS_ERROR << "ConstructCellManager does not support FESpace<0>";
    }

    /// Get the vector of function indices
    virtual std::vector<std::size_t> FunctionIndices() const {return std::vector<std::size_t> {mFunctionId};}

    /// Check the compatibility between boundaries of two FESpacees
    virtual bool CheckBoundaryCompatibility(const FESpaceType& rFESpace1, const BoundarySide& side1,
                                            const FESpaceType& rFESpace2, const BoundarySide& side2) const
    {
        return true;
    }

    /// Validate the FESpace before using
    virtual bool Validate() const
    {
        return true;
    }

    // /// Construct the boundary FESpace based on side
    // virtual typename FESpace<-1>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    // {
    //     return NULL;
    // }

    /// Compare the two FESpacees in terms of its parametric information.
    virtual bool IsCompatible(const FESpaceType& rOtherFESpace) const
    {
        return true;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "FESpace<0>";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << " Function Indices: " << mFunctionId
                 << std::endl;
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
 * Template specific instantiation for -1-D FESpace to terminate the compilation
 */
template<typename TLocalCoordinateType>
class FESpace<-1, TLocalCoordinateType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(FESpace);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const FESpace> ConstPointer;
#endif

    typedef FESpace<-1, TLocalCoordinateType> FESpaceType;

    /// Default constructor
    FESpace() {}

    /// Destructor
    virtual ~FESpace() {}

    /// Get the number of basis functions defined over the FESpace
    virtual std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the FESpace in specific direction
    virtual std::size_t Order(std::size_t i) const
    {
        return 0;
    }

    /// Get the string describing the type of the FESpace
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the FESpace
    static std::string StaticType()
    {
        return "FESpace<-1>D";
    }

    /// Overload comparison operator
    virtual bool operator==(const FESpaceType & rOther) const
    {
        return true;
    }

    /// Check the compatibility between boundaries of two FESpacees
    virtual bool CheckBoundaryCompatibility(const FESpaceType& rFESpace1, const BoundarySide& side1,
                                            const FESpaceType& rFESpace2, const BoundarySide& side2) const
    {
        return true;
    }

    /// Validate the FESpace before using
    virtual bool Validate() const
    {
        return true;
    }

    // /// Construct the boundary FESpace based on side
    // virtual typename FESpace<-2>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    // {
    //     return NULL;
    // }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "FESpace<-1>";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};

/**
 * Template specific instantiation for -2-D FESpace to terminate the compilation
 */
template<typename TLocalCoordinateType>
class FESpace<-2, TLocalCoordinateType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(FESpace);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const FESpace> ConstPointer;
#endif

    typedef FESpace<-2, TLocalCoordinateType> FESpaceType;

    /// Default constructor
    FESpace() {}

    /// Destructor
    virtual ~FESpace() {}

    /// Get the number of basis functions defined over the FESpace
    virtual std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the FESpace in specific direction
    virtual std::size_t Order(std::size_t i) const
    {
        return 0;
    }

    /// Get the string describing the type of the FESpace
    virtual std::string Type() const
    {
        return StaticType();
    }

    /// Get the string describing the type of the FESpace
    static std::string StaticType()
    {
        return "FESpace<-2>D";
    }

    /// Overload comparison operator
    virtual bool operator==(const FESpaceType& rOther) const
    {
        return true;
    }

    /// Check the compatibility between boundaries of two FESpacees
    virtual bool CheckBoundaryCompatibility(const FESpaceType& rFESpace1, const BoundarySide& side1,
                                            const FESpaceType& rFESpace2, const BoundarySide& side2) const
    {
        return true;
    }

    /// Validate the FESpace before using
    virtual bool Validate() const
    {
        return true;
    }

    // /// Construct the boundary FESpace based on side
    // virtual typename FESpace<-3>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    // {
    //     return NULL;
    // }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "FESpace<-2>";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};

/// output stream function
template<int TDim, typename TLocalCoordinateType>
inline std::ostream& operator <<(std::ostream& rOStream, const FESpace<TDim, TLocalCoordinateType>& rThis)
{
    rOStream << "-------------Begin FESpaceInfo-------------" << std::endl;
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "-------------End FESpaceInfo-------------";
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_FESPACE_H_INCLUDED defined
