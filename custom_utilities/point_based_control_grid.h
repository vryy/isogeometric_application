//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 27 Nov 2017 $
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
    typename BaseType::Pointer Clone() const override
    {
        return PointBasedControlGrid::Pointer(new PointBasedControlGrid(mrVariable, mpFESpace));
    }

    /// Access the underlying FESpace
    typename FESpaceType::Pointer pFESpace() {return mpFESpace;}

    /// Access the underlying FESpace
    typename FESpaceType::ConstPointer pFESpace() const {return mpFESpace;}

    /// Get the size of underlying data
    std::size_t Size() const override
    {
        return mpFESpace->TotalNumber();
    }

    /// Get the size of underlying data
    std::size_t size() const override
    {
        return mpFESpace->TotalNumber();
    }

    /// Get the data at specific point
    /// It is noted that the return value is unweighted one
    DataType GetData(std::size_t i) const override
    {
        // TODO Get and Set data in the sequential manner can be expensive if the underlying FESPace uses std::set to store the basis functions.
        // It is suggested to implement the iterator for get and set the values.
        return (*mpFESpace)[i]->GetValue(mrVariable) / (*mpFESpace)[i]->Weight();
    }

    /// Set the data at specific point
    /// It is noted that the setting value is unweighted one
    void SetData(std::size_t i, const DataType& value) override
    {
        // TODO see comment in GetData
        (*mpFESpace)[i]->SetValue(mrVariable, value * (*mpFESpace)[i]->Weight());
    }

    // overload operator []
    DataType& operator[] (std::size_t i) override
    {
        // TODO see comment in GetData
        return (*mpFESpace)[i]->GetValue(mrVariable);
    }

    // overload operator []
    const DataType& operator[] (std::size_t i) const override
    {
        // TODO see comment in GetData
        return (*mpFESpace)[i]->GetValue(mrVariable);
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Point-Based Control Grid " << BaseType::Name() << "[" << Size() << "]";
    }

    void PrintData(std::ostream& rOStream) const override
    {
        // print out the control values
        // TODO we shall use the iterator here for more efficiency, especially for hierarchical B-Splines
        for (std::size_t i = 0; i < this->size(); ++i)
        {
            rOStream << this->GetData(i) << std::endl;
        }
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
    typename BaseType::Pointer Clone() const override
    {
        return PointBasedControlGrid::Pointer(new PointBasedControlGrid(mrVariable, mpFESpace));
    }

    /// Access the underlying FESpace
    typename FESpaceType::Pointer pFESpace() {return mpFESpace;}

    /// Access the underlying FESpace
    typename FESpaceType::ConstPointer pFESpace() const {return mpFESpace;}

    /// Get the size of underlying data
    std::size_t Size() const override
    {
        return mpFESpace->TotalNumber();
    }

    /// Get the size of underlying data
    std::size_t size() const override
    {
        return mpFESpace->TotalNumber();
    }

    /// Get the data at specific point
    DataType GetData(std::size_t i) const override
    {
        // TODO Get and Set data in the sequential manner can be expensive if the underlying FESPace uses set to store the basis functions. It is suggested to implement the iterator for get and set the values.
        return (*mpFESpace)[i]->GetValue(mrVariable);
    }

    // overload operator []
    DataType& operator[] (std::size_t i) override
    {
        // TODO see comment in GetData
        return (*mpFESpace)[i]->GetValue(mrVariable);
    }

    // overload operator []
    const DataType& operator[] (std::size_t i) const override
    {
        // TODO see comment in GetData
        return (*mpFESpace)[i]->GetValue(mrVariable);
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Point-Based Control Grid " << BaseType::Name() << "[" << Size() << "]";
    }

    void PrintData(std::ostream& rOStream) const override
    {
        // print out the control values
        // TODO we shall use the iterator here for more efficiency, especially for hierarchical B-Splines
        for (std::size_t i = 0; i < this->size(); ++i)
        {
            rOStream << this->GetData(i) << std::endl;
        }
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
