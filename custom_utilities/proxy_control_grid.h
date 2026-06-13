//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 13 Jun 2026 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_PROXY_CONTROL_GRID_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_PROXY_CONTROL_GRID_H_INCLUDED

// System includes

// External includes

// Project includes
#include "control_grid.h"

namespace Kratos
{

/**
 * A proxy control grid provides access to the unweighted control value through the control grid
 */
template<typename TControlValueType>
class ProxyControlGrid : public ControlGrid<typename TControlValueType::DataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ProxyControlGrid);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const ProxyControlGrid> ConstPointer;
#endif

    /// Type definition
    typedef ControlGrid<TControlValueType> ProxyType;
    typedef ControlGrid<typename TControlValueType::DataType> BaseType;
    typedef typename BaseType::DataType DataType;

    /// Default constructor
    ProxyControlGrid(typename ProxyType::Pointer pControlGrid)
    : BaseType(), mpControlGrid(pControlGrid)
    {}

    /// Destructor
    ~ProxyControlGrid() override {}

    /// Create a new control grid pointer
    static typename BaseType::Pointer Create(typename ProxyType::Pointer pControlGrid)
    {
        return typename BaseType::Pointer(new ProxyControlGrid(pControlGrid));
    }

    /// Create the clone
    typename BaseType::Pointer Clone() const
    {
        typename ProxyType::Pointer pNewControlGrid = mpControlGrid->Clone();
        if (pNewControlGrid == nullptr)
            KRATOS_ERROR << "Error cloning " << this->Type();
        return Create(pNewControlGrid);
    }

    /// Get the size of underlying data
    std::size_t Size() const override
    {
        return mpControlGrid->Size();
    }

    /// Get the data at specific point
    DataType GetData(std::size_t i) const override
    {
        return mpControlGrid->GetData(i).V();
    }

    /// Get the type of the control grid
    std::string Type() const override
    {
        return "ProxyControlGrid_" + mpControlGrid->Type();
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Proxy Control Grid of ";
        mpControlGrid->PrintInfo(rOStream);
    }

    void PrintData(std::ostream& rOStream) const override
    {
        mpControlGrid->PrintData(rOStream);
    }

private:

    typename ProxyType::Pointer mpControlGrid;

    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        // TODO
    }

    virtual void load(Serializer& rSerializer)
    {
        // TOOD
    }
    ///@}
};

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_PROXY_CONTROL_GRID_H_INCLUDED defined
