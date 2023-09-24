//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 30 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_UNSTRUCTURED_CONTROL_GRID_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_UNSTRUCTURED_CONTROL_GRID_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/control_grid.h"

namespace Kratos
{

/**
This class is an ordinary container to keep the control values.
 */
template<typename TDataType>
class UnstructuredControlGrid : public ControlGrid<TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(UnstructuredControlGrid);

    /// Type definition
    typedef ControlGrid<TDataType> BaseType;
    typedef std::vector<TDataType> DataContainerType;

    /// Default constructor
    UnstructuredControlGrid(std::size_t size) : BaseType() {mData.resize(size);}

    /// Constructor with name
    UnstructuredControlGrid(const std::string& Name, std::size_t size) : BaseType(Name) {mData.resize(size);}

    /// Destructor
    virtual ~UnstructuredControlGrid() {}

    /// Create a new control grid pointer
    static typename UnstructuredControlGrid::Pointer Create(std::size_t size)
    {
        return typename UnstructuredControlGrid::Pointer(new UnstructuredControlGrid(size));
    }

    /// Clone this grid function
    typename BaseType::Pointer Clone() const override
    {
        UnstructuredControlGrid::Pointer pNewControlGrid = Create(size());
        for (std::size_t i = 0; i < size(); ++i)
            pNewControlGrid->SetData(i, GetData(i));
        return pNewControlGrid;
    }

    /// Get the size of underlying data
    std::size_t Size() const override {return mData.size();}

    /// Get the size of underlying data
    std::size_t size() const override {return mData.size();}

    /// Resize the underlying container
    void resize(std::size_t new_size) {mData.resize(new_size);}

    /// Get the data at specific point
    TDataType GetData(std::size_t i) const override {return mData[i];}

    /// Set the data at specific point
    /// Be careful with this method. You can destroy the coherency of internal data.
    void SetData(std::size_t i, const TDataType& value) override {mData[i] = value;}

    /// overload operator []
    TDataType& operator[] (std::size_t i) override {return mData[i];}

    /// overload operator []
    const TDataType& operator[] (std::size_t i) const override {return mData[i];}

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Unstructured Control Grid " << BaseType::Name() << "[" << Size() << "]";
    }

    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "Data:";
        for (std::size_t i = 0; i < mData.size(); ++i)
            rOStream << " " << mData[i];
    }

private:

    DataContainerType mData;
};

/// output stream function
template<typename TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const UnstructuredControlGrid<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_UNSTRUCTURED_CONTROL_GRID_H_INCLUDED defined
