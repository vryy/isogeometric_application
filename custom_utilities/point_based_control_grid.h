//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 27 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_POINT_BASED_CONTROL_GRID_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_POINT_BASED_CONTROL_GRID_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/control_grid.h"

namespace Kratos
{

/**
The point-based control grid allows to access the control values of point-based Splines. It relies on the FESpace to provide the basis function necessary for value extraction.
It is designed to be the control grid for point-based Splines, e.g. hierarchical B-Splines, T-Splines, ...
 */
template<typename TVariableType, class TFESpaceType>
class PointBasedControlGrid : public ControlGrid<typename TVariableType::Type>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PointBasedControlGrid);

    /// Type definition
    typedef ControlGrid<typename TVariableType::Type> BaseType;
    typedef typename BaseType::DataType DataType; // which is the same as TVariableType::Type
    typedef TVariableType VariableType;
    typedef TFESpaceType FESpaceType;

    /// Constructor with Variable and FESpace
    PointBasedControlGrid(const VariableType& rVariable, typename FESpaceType::Pointer pFESpace)
    : BaseType(rVariable.Name()), mrVariable(rVariable), mpFESpace(pFESpace) {}

    /// Destructor
    virtual ~PointBasedControlGrid() {}

    /// Create a new control grid pointer
    static typename PointBasedControlGrid::Pointer Create(const VariableType& rVariable, typename FESpaceType::Pointer pFESpace)
    {
        return typename PointBasedControlGrid::Pointer(new PointBasedControlGrid(rVariable, pFESpace));
    }

    /// Clone this grid function
    virtual typename BaseType::Pointer Clone() const
    {
        return PointBasedControlGrid::Pointer(new PointBasedControlGrid(mrVariable, mpFESpace));
    }

    /// Access the underlying FESpace
    typename FESpaceType::Pointer pFESpace() {return mpFESpace;}

    /// Access the underlying FESpace
    typename FESpaceType::ConstPointer pFESpace() const {return mpFESpace;}

    /// Get the size of underlying data
    virtual std::size_t Size() const
    {
        return mpFESpace->TotalNumber();
    }

    /// Get the size of underlying data
    virtual std::size_t size() const
    {
        return mpFESpace->TotalNumber();
    }

    /// Get the data at specific point
    virtual DataType GetData(const std::size_t& i) const
    {
        // TODO Get and Set data in the sequential manner can be expensive if the underlying FESPace uses std::set to store the basis functions.
        // It is suggested to implement the iterator for get and set the values.
        return (*mpFESpace)[i]->GetValue(mrVariable) / (*mpFESpace)[i]->Weight();
    }

    /// Set the data at specific point
    /// Be careful with this method. You can destroy the coherency of internal data.
    virtual void SetData(const std::size_t& i, const DataType& value)
    {
        // TODO see comment in GetData
        (*mpFESpace)[i]->SetValue(mrVariable, value * (*mpFESpace)[i]->Weight());
    }

    // overload operator []
    virtual DataType& operator[] (const std::size_t& i)
    {
        // TODO see comment in GetData
        return (*mpFESpace)[i]->GetValue(mrVariable);
    }

    // overload operator []
    virtual const DataType& operator[] (const std::size_t& i) const
    {
        // TODO see comment in GetData
        return (*mpFESpace)[i]->GetValue(mrVariable);
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Point-Based Control Grid " << BaseType::Name() << "[" << Size() << "]";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        // print out the control values
        // TODO we shall use the iterator here for more efficiency, especially for hierarchical B-Splines
        for (std::size_t i = 0; i < this->size(); ++i)
            rOStream << this->GetData(i) << std::endl;
    }

private:

    const VariableType& mrVariable;
    typename FESpaceType::Pointer mpFESpace;
};


/// Partial template specialization for ControlPoint variable
template<class TFESpaceType>
class PointBasedControlGrid<Variable<ControlPoint<double> >, TFESpaceType> : public ControlGrid<ControlPoint<double> >
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PointBasedControlGrid);

    /// Type definition
    typedef ControlGrid<ControlPoint<double> > BaseType;
    typedef typename BaseType::DataType DataType; // which is the same as TVariableType::Type
    typedef Variable<ControlPoint<double> > VariableType;
    typedef TFESpaceType FESpaceType;

    /// Constructor with Variable and FESpace
    PointBasedControlGrid(const VariableType& rVariable, typename FESpaceType::Pointer pFESpace)
    : BaseType(rVariable.Name()), mrVariable(rVariable), mpFESpace(pFESpace) {}

    /// Destructor
    virtual ~PointBasedControlGrid() {}

    /// Create a new control grid pointer
    static typename PointBasedControlGrid::Pointer Create(const VariableType& rVariable, typename FESpaceType::Pointer pFESpace)
    {
        return typename PointBasedControlGrid::Pointer(new PointBasedControlGrid(rVariable, pFESpace));
    }

    /// Clone this grid function
    virtual typename BaseType::Pointer Clone() const
    {
        return PointBasedControlGrid::Pointer(new PointBasedControlGrid(mrVariable, mpFESpace));
    }

    /// Access the underlying FESpace
    typename FESpaceType::Pointer pFESpace() {return mpFESpace;}

    /// Access the underlying FESpace
    typename FESpaceType::ConstPointer pFESpace() const {return mpFESpace;}

    /// Get the size of underlying data
    virtual std::size_t Size() const
    {
        return mpFESpace->TotalNumber();
    }

    /// Get the size of underlying data
    virtual std::size_t size() const
    {
        return mpFESpace->TotalNumber();
    }

    /// Get the data at specific point
    virtual DataType GetData(const std::size_t& i) const
    {
        // TODO Get and Set data in the sequential manner can be expensive if the underlying FESPace uses set to store the basis functions. It is suggested to implement the iterator for get and set the values.
        return (*mpFESpace)[i]->GetValue(mrVariable);
    }

    // overload operator []
    virtual DataType& operator[] (const std::size_t& i)
    {
        // TODO see comment in GetData
        return (*mpFESpace)[i]->GetValue(mrVariable);
    }

    // overload operator []
    virtual const DataType& operator[] (const std::size_t& i) const
    {
        // TODO see comment in GetData
        return (*mpFESpace)[i]->GetValue(mrVariable);
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Point-Based Control Grid " << BaseType::Name() << "[" << Size() << "]";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        // print out the control values
        // TODO we shall use the iterator here for more efficiency, especially for hierarchical B-Splines
        for (std::size_t i = 0; i < this->size(); ++i)
            rOStream << this->GetData(i) << std::endl;
    }

private:

    const VariableType& mrVariable;
    typename FESpaceType::Pointer mpFESpace;
};

/// output stream function
template<class TVariableType, class TFESpaceType>
inline std::ostream& operator <<(std::ostream& rOStream, const PointBasedControlGrid<TVariableType, TFESpaceType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_POINT_BASED_CONTROL_GRID_H_INCLUDED defined
