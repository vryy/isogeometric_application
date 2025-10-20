//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 Nov 2017 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_GRID_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_GRID_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"

namespace Kratos
{

/**
This class is an abstract container to keep the control values.
TODO implement iterator to iterate through control values.
 */
template<typename TDataType>
class ControlGrid
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ControlGrid);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const ControlGrid> ConstPointer;
#endif

    /// Type definition
    typedef TDataType DataType;

    /// Default constructor
    ControlGrid() : mName("UNKNOWN") {}

    /// Constructor with name
    ControlGrid(const std::string& Name) : mName(Name) {}

    /// Destructor
    virtual ~ControlGrid() {}

    /// Create a new control grid pointer
    static typename ControlGrid::Pointer Create() {return typename ControlGrid::Pointer(new ControlGrid());}

    /// Overload assignment operator
    ControlGrid<TDataType>& operator=(const ControlGrid<TDataType>& rOther)
    {
        this->mName = rOther.mName;
        return *this;
    }

    /// Create the clone
    virtual typename ControlGrid<TDataType>::Pointer Clone() const
    {
        typename ControlGrid<TDataType>::Pointer pNewControlGrid = typename ControlGrid<TDataType>::Pointer(new ControlGrid<TDataType>());
        *pNewControlGrid = *this;
        return pNewControlGrid;
    }

    /// Set the name
    void SetName(const std::string& Name) {mName = Name;}

    /// Get the name
    const std::string& Name() const {return mName;}

    /// Get the size of underlying data
    virtual std::size_t Size() const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Get the size of underlying data
    virtual std::size_t size() const
    {
        return this->Size();
    }

    /// Get the data at specific point
    virtual TDataType GetData(std::size_t i) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Set the data at specific point
    /// Be careful with this method. You can destroy the coherency of internal data.
    virtual void SetData(std::size_t i, const TDataType& value)
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// overload operator []
    virtual TDataType& operator[] (std::size_t i)
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// overload operator []
    virtual const TDataType& operator[] (std::size_t i) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    // /// extract the sub grid, given local indices of the control values
    // void ControlGrid::Pointer ExtractSubGrid(const std::vector<std::size_t>& local_ids) const
    // {
    //     ControlGrid::Pointer pSubControlGrid = Create();
    //     for (std::size_t i = 0; i < local_ids.size(); ++i)
    //     {

    //     }
    //     return pSubControlGrid;
    // }

    /// Copy the data the other grid. The size of two grids must be equal.
    virtual void CopyFrom(const ControlGrid<TDataType>& rOther)
    {
        if (rOther.Size() != this->Size())
        {
            KRATOS_ERROR << "The size of the input grid is incompatible";
        }

        for (std::size_t i = 0; i < this->size(); ++i)
        {
            this->SetData(i, rOther.GetData(i));
        }
    }

    /// Copy the data the other grid. The size of two grids must be equal.
    virtual void CopyFrom(const typename ControlGrid<TDataType>::Pointer pOther)
    {
        this->CopyFrom(*pOther);
    }

    /// Copy the data the other grid. In the case that the source has different size, the grid is resized.
    virtual void ResizeAndCopyFrom(ControlGrid<TDataType>& rOther)
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Copy the data the other grid. In the case that the source has different size, the grid is resized.
    virtual void ResizeAndCopyFrom(const typename ControlGrid<TDataType>::Pointer pOther)
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Get the type of the control grid
    virtual std::string Type() const
    {
        return "ControlGrid";
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Control Grid " << Name() << "[" << Size() << "]";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:

    std::string mName;

    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        // std::cout << "Serialization - calling ControlGrid " << this->Type() << " " << __FUNCTION__ << std::endl;
        rSerializer.save( "Name", mName );
    }

    virtual void load(Serializer& rSerializer)
    {
        // std::cout << "Serialization - calling ControlGrid " << this->Type() << " " << __FUNCTION__ << std::endl;
        rSerializer.load( "Name", mName );
    }
    ///@}
};

/// output stream function
template<typename TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const ControlGrid<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_GRID_H_INCLUDED defined
