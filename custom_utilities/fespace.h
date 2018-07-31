//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_FESPACE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_FESPACE_H_INCLUDED

// System includes
#include <vector>

// External includes
#include <boost/enable_shared_from_this.hpp>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "custom_utilities/iga_define.h"
#include "custom_utilities/nurbs/cell.h"
#include "custom_utilities/nurbs/cell_manager.h"


namespace Kratos
{

/**
An FESpace is a collection of shape function defined over the parametric domain. An isogeometric FESpace can be a NURBS FESpace, a hierarchical NURBS FESpace, or a T-Splines FESpace.
 */
template<int TDim>
class FESpace
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(FESpace);

    /// Type definition
    typedef CellManager<Cell> cell_container_t;

    /// Default constructor
    FESpace() {}

    /// Destructor
    virtual ~FESpace()
    {
    }

    /// Helper to create new BSplinesFESpace pointer
    static typename FESpace<TDim>::Pointer Create()
    {
        return typename FESpace<TDim>::Pointer(new FESpace());
    }

    /// Get the number of basis functions defined over the FESpace
    virtual const std::size_t TotalNumber() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Get the order of the FESpace in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
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

    /// Get the values of the basis function i at point xi
    virtual void GetValue(double& v, const std::size_t& i, const std::vector<double>& xi) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Get the values of the basis function i at point xi
    double GetValue(const std::size_t& i, const std::vector<double>& xi) const
    {
        double v;
        this->GetValue(v, i, xi);
        return v;
    }

    /// Get the values of the basis functions at point xi
    virtual void GetValue(std::vector<double>& values, const std::vector<double>& xi) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Get the values of the basis functions at point xi
    std::vector<double> GetValue(const std::vector<double>& xi) const
    {
        std::vector<double> values;
        this->GetValue(values, xi);
        return values;
    }

    ///////////////

    /// Get the derivatives of the basis function i at point xi
    virtual void GetDerivative(std::vector<double>& values, const std::size_t& i, const std::vector<double>& xi) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Get the derivatives of the basis function i at point xi
    std::vector<double> GetDerivative(const std::size_t& i, const std::vector<double>& xi) const
    {
        std::vector<double> values;
        this->GetDerivative(values, i, xi);
        return values;
    }

    /// Get the derivatives of the basis functions at point xi
    /// the output values has the form of values[func_index][dim_index]
    virtual void GetDerivative(std::vector<std::vector<double> >& values, const std::vector<double>& xi) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Get the derivatives of the basis functions at point xi
    /// the return values has the form of values[func_index][dim_index]
    std::vector<std::vector<double> > GetDerivative(const std::vector<double>& xi) const
    {
        std::vector<std::vector<double> > values;
        this->GetDerivative(values);
        return values;
    }

    ///////////////

    /// Get the values and derivatives of the basis functions at point xi
    /// the output derivatives has the form of values[func_index][dim_index]
    virtual void GetValueAndDerivative(std::vector<double>& values, std::vector<std::vector<double> >& derivatives, const std::vector<double>& xi) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Reset all the dof numbers for each grid function to -1.
    virtual void ResetFunctionIndices()
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Reset the function indices to a given values.
    /// This is useful when assigning the id for the boundary patch.
    virtual void ResetFunctionIndices(const std::vector<std::size_t>& func_indices)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Enumerate the dofs of each grid function. The enumeration algorithm is pretty straightforward.
    /// If the dof does not have pre-existing value, which assume it is -1, it will be assigned the incremental value.
    virtual std::size_t& Enumerate(std::size_t& start)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Access the function indices (aka global ids)
    virtual std::vector<std::size_t> FunctionIndices() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Update the function indices using a map. The map shall be the mapping from old index to new index.
    virtual void UpdateFunctionIndices(const std::map<std::size_t, std::size_t>& indices_map)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Get the first equation_id in this space
    virtual std::size_t GetFirstEquationId() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Get the last equation_id in this space
    virtual std::size_t GetLastEquationId() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Return the local id of a given global id
    std::size_t LocalId(const std::size_t& global_id) const
    {
        std::map<std::size_t, std::size_t>::const_iterator it = mGlobalToLocal.find(global_id);

        if (it == mGlobalToLocal.end())
        {
            KRATOS_WATCH(TDim)
            KRATOS_WATCH(global_id)
            std::cout << "mGlobalToLocal:";
            for(std::map<std::size_t, std::size_t>::const_iterator it2 = mGlobalToLocal.begin(); it2 != mGlobalToLocal.end(); ++it2)
                std::cout << " " << it2->first << "->" << it2->second;
            std::cout << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, "The global id does not exist in global_to_local map", "")
        }

        return it->second;
    }

    /// Return the local ids of given global ids
    std::vector<std::size_t> LocalId(const std::vector<std::size_t>& global_ids) const
    {
        std::vector<std::size_t> local_ids(global_ids.size());
        for (std::size_t i = 0; i < global_ids.size(); ++i)
            local_ids[i] = this->LocalId(global_ids[i]);
        return local_ids;
    }

    /// Check the compatibility between boundaries of two FESpacees
    virtual bool CheckBoundaryCompatibility(const FESpace<TDim>& rFESpace1, const BoundarySide& side1,
            const FESpace<TDim>& rFESpace2, const BoundarySide& side2) const
    {
        typename FESpace<TDim-1>::Pointer pBFESpace1 = rFESpace1.ConstructBoundaryFESpace(side1);
        typename FESpace<TDim-1>::Pointer pBFESpace2 = rFESpace1.ConstructBoundaryFESpace(side2);

        return (*pBFESpace1) == (*pBFESpace2);
    }

    /// Validate the FESpace before using
    virtual bool Validate() const
    {
        return true;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Extract the index of the functions on the boundary
    virtual std::vector<std::size_t> ExtractBoundaryFunctionIndices(const BoundarySide& side) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Extract the index of the functions on the boundary down to some level
    virtual std::vector<std::size_t> ExtractBoundaryFunctionIndices(const BoundarySide& side, const std::size_t& level) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Assign the index for the functions on the boundary
    virtual void AssignBoundaryFunctionIndices(const BoundarySide& side, const std::vector<std::size_t>& func_indices)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Construct the boundary FESpace based on side
    virtual typename FESpace<TDim-1>::Pointer ConstructBoundaryFESpace(const BoundarySide& side) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Construct the boundary FESpace based on side and local relative configuration between two patches
    virtual typename FESpace<TDim-1>::Pointer ConstructBoundaryFESpace(const BoundarySide& side,
        const std::map<std::size_t, std::size_t>& local_parameter_map, const std::vector<BoundaryDirection>& directions) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Compare the two FESpacees in terms of its parametric information.
    virtual bool IsCompatible(const FESpace<TDim>& rOtherFESpace) const
    {
        return false;
    }

    /// Overload comparison operator
    virtual bool operator==(const FESpace<TDim>& rOther) const
    {
        return this->IsCompatible(rOther);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Fast function to get the opposite boundary side
    static BoundarySide OppositeBoundarySide(const BoundarySide& side)
    {
        if (side == _BLEFT_) return _BRIGHT_;
        else if (side == _BRIGHT_) return _BLEFT_;
        else if (side == _BTOP_) return _BBOTTOM_;
        else if (side == _BBOTTOM_) return _BTOP_;
        else if (side == _BFRONT_) return _BBACK_;
        else if (side == _BBACK_) return _BFRONT_;
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid boundary side", side)
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Create the cell manager for all the cells in the support domain of the FESpace
    virtual typename cell_container_t::Pointer ConstructCellManager() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Overload assignment operator
    FESpace<TDim>& operator=(const FESpace<TDim>& rOther)
    {
        mGlobalToLocal = rOther.mGlobalToLocal;
        return *this;
    }

    /// Clone this FESpace, this is a deep copy operation
    virtual typename FESpace<TDim>::Pointer Clone() const
    {
        typename FESpace<TDim>::Pointer pNewFESpace = typename FESpace<TDim>::Pointer(new FESpace<TDim>());
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
            rOStream << " " << func_indices[i];
    }

protected:

    std::map<std::size_t, std::size_t> mGlobalToLocal;

private:

    /// Serializer
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
//        rSerializer.save( "mFunctionsIds", mFunctionsIds ); // can't save std::vector
//        rSerializer.save( "mGlobalToLocal", mGlobalToLocal ); // can't save std::map
    }

    virtual void load(Serializer& rSerializer)
    {
//        rSerializer.load( "mFunctionsIds", mFunctionsIds );
//        rSerializer.load( "mGlobalToLocal", mGlobalToLocal );
    }
};



/**
 * Template specific instantiation for null-D FESpace to terminate the compilation
 */
template<>
class FESpace<0>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(FESpace);

    // Type definitions
    typedef CellManager<Cell> cell_container_t;

    /// Default constructor
    FESpace() : mFunctionId(-1) {}

    /// Destructor
    virtual ~FESpace() {}

    /// Get the number of basis functions defined over the FESpace
    virtual const std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the FESpace in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
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
    virtual bool operator==(const FESpace<0>& rOther) const
    {
        return true;
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

    /// Get the vector of function indices
    virtual std::vector<std::size_t> FunctionIndices() const {return std::vector<std::size_t>{mFunctionId};}

    /// Check the compatibility between boundaries of two FESpacees
    virtual bool CheckBoundaryCompatibility(const FESpace<0>& rFESpace1, const BoundarySide& side1,
            const FESpace<0>& rFESpace2, const BoundarySide& side2) const
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

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "FESpace<0>";
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
 * Template specific instantiation for -1-D FESpace to terminate the compilation
 */
template<>
class FESpace<-1>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(FESpace);

    /// Default constructor
    FESpace() {}

    /// Destructor
    virtual ~FESpace() {}

    /// Get the number of basis functions defined over the FESpace
    virtual const std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the FESpace in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
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
    virtual bool operator==(const FESpace<-1>& rOther) const
    {
        return true;
    }

    /// Check the compatibility between boundaries of two FESpacees
    virtual bool CheckBoundaryCompatibility(const FESpace<-1>& rFESpace1, const BoundarySide& side1,
            const FESpace<-1>& rFESpace2, const BoundarySide& side2) const
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
template<>
class FESpace<-2>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(FESpace);

    /// Default constructor
    FESpace() {}

    /// Destructor
    virtual ~FESpace() {}

    /// Get the number of basis functions defined over the FESpace
    virtual const std::size_t TotalNumber() const
    {
        return 0;
    }

    /// Get the order of the FESpace in specific direction
    virtual const std::size_t Order(const std::size_t& i) const
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
    virtual bool operator==(const FESpace<-2>& rOther) const
    {
        return true;
    }

    /// Check the compatibility between boundaries of two FESpacees
    virtual bool CheckBoundaryCompatibility(const FESpace<-2>& rFESpace1, const BoundarySide& side1,
            const FESpace<-2>& rFESpace2, const BoundarySide& side2) const
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
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const FESpace<TDim>& rThis)
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

#undef DEBUG_DESTROY

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_FESPACE_H_INCLUDED defined

