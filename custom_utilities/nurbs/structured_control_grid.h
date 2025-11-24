//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 Nov 2017 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_REGULAR_CONTROL_GRID_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_REGULAR_CONTROL_GRID_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/iga_define.h"
#include "custom_utilities/control_grid.h"
#include "custom_utilities/nurbs/bsplines_indexing_utility.h"

namespace Kratos
{

/**
Base class for control value container by a regular grid
*/
template<typename TDataType>
class BaseStructuredControlGrid : public ControlGrid<TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BaseStructuredControlGrid);

    /// Type definition
    typedef ControlGrid<TDataType> BaseType;
    typedef std::vector<TDataType> DataContainerType;

    /// Default constructor
    BaseStructuredControlGrid() : BaseType() {}

    /// Constructor with name
    BaseStructuredControlGrid(const std::string& Name) : BaseType(Name) {}

    /// Copy constructor
    BaseStructuredControlGrid(const BaseStructuredControlGrid& rOther)
    : BaseType(rOther), mData(rOther.mData)
    {}

    /// Destructor
    ~BaseStructuredControlGrid() override {}

    /************************************/
    /********* INHERIT UPSTREAM *********/
    /************************************/

    /// Get the size of underlying data
    std::size_t Size() const override {return mData.size();}

    /// Get the size of underlying data
    std::size_t size() const override {return mData.size();}

    /// Get the data at specific point
    TDataType GetData(std::size_t i) const override {return mData[i];}

    /// Set the data at specific point
    /// Be careful with this method. You can destroy the coherency of internal data.
    void SetData(std::size_t i, const TDataType& value) override {mData[i] = value;}

    /// overload operator []
    TDataType& operator[] (std::size_t i) override {return mData[i];}

    /// overload operator []
    const TDataType& operator[] (std::size_t i) const override {return mData[i];}

    /// Overload assignment operator
    BaseStructuredControlGrid<TDataType>& operator=(const BaseStructuredControlGrid<TDataType>& rOther)
    {
        BaseType::operator=(rOther);
        this->mData = rOther.mData;
        return *this;
    }

    /// Create the clone
    typename ControlGrid<TDataType>::Pointer Clone() const override
    {
        typename BaseStructuredControlGrid<TDataType>::Pointer pNewControlGrid = typename BaseStructuredControlGrid<TDataType>::Pointer(new BaseStructuredControlGrid<TDataType>());
        *pNewControlGrid = *this;
        return pNewControlGrid;
    }

    /************************************/
    /****** EXCLUSIVE SUBROUTINES *******/
    /************************************/

    /// Get the global index of the grid value providing the index vector
    virtual std::size_t GetIndex(const std::vector<std::size_t>& index) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// resize the underlying container
    void Resize(std::size_t new_size) {mData.resize(new_size);}

    /// resize the underlying container
    void resize(std::size_t new_size) {mData.resize(new_size);}

    /// Access the underlying data
    DataContainerType& Data() {return mData;}

    /// Access the underlying data
    const DataContainerType& Data() const {return mData;}

    /************************************/
    /******** SUCCEED DOWNSTREAM ********/
    /************************************/

    /// Copy the data the other grid. The size of two grids must be equal.
    void CopyFrom(const ControlGrid<TDataType>& rOther) override
    {
        BaseType::CopyFrom(rOther);
    }

    /// Copy the data the other grid. The size of two grids must be equal.
    void CopyFrom(const typename ControlGrid<TDataType>::Pointer pOther) override
    {
        BaseType::CopyFrom(pOther);
    }

    /// Copy the data the other grid. In the case that the source has different size, the grid is resized.
    void ResizeAndCopyFrom(ControlGrid<TDataType>& rOther) override
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Copy the data the other grid. In the case that the source has different size, the grid is resized.
    void ResizeAndCopyFrom(const typename ControlGrid<TDataType>::Pointer pOther) override
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Reverse the control grid in specific dimension
    virtual void Reverse(int idir)
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Transpose the control grid in
    virtual void Transpose(int idir, int jdir)
    {
        KRATOS_ERROR << "Error calling base class function";
    }

    /// Create a connectivity matrix for the structured control grid
    virtual void CreateConnectivity(std::size_t offset, std::vector<std::vector<std::size_t> >& connectivities) const
    {
        KRATOS_ERROR << "Error calling base class function";
    }

private:

    DataContainerType mData;

    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType);
        rSerializer.save("Data", mData);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType);
        rSerializer.load("Data", mData);
    }
    ///@}
};

/**
Class for control value container by a regular grid
*/
template<int TDim, typename TDataType> class StructuredControlGrid;

template<typename TDataType>
class StructuredControlGrid<0, TDataType> : public BaseStructuredControlGrid<TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(StructuredControlGrid);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const StructuredControlGrid> ConstPointer;
#endif

    // type definitions
    typedef BaseStructuredControlGrid<TDataType> BaseType;
    typedef typename BaseType::DataContainerType DataContainerType;
    typedef typename BaseType::DataType DataType;

    /// Constructor with size
    StructuredControlGrid(const std::vector<std::size_t>& sizes) : BaseType()
    {
        BaseType::Data().resize(1);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType());
    }

    /// Constructor with size
    StructuredControlGrid(std::size_t n) : BaseType()
    {
        BaseType::Data().resize(1);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType());
    }

    /// Destructor
    ~StructuredControlGrid() override {}

    /// Create a new control grid pointer
    static typename StructuredControlGrid<0, TDataType>::Pointer Create(const std::vector<std::size_t>& sizes)
    {
        return typename StructuredControlGrid<0, TDataType>::Pointer(new StructuredControlGrid<0, TDataType>(sizes));
    }

    /// resize the grid
    void Resize(std::size_t new_size)
    {
        // DO NOTHING
    }

    /// resize the grid
    void resize(std::size_t new_size)
    {
        // DO NOTHING
    }

    /// Get the size of underlying data
    std::size_t Size() const {return BaseType::Data().size();}

    /// Get the size of the grid is specific dimension
    std::size_t Size(std::size_t dim) const {return 1;}

    /// Get the global index of the grid value providing the index vector
    std::size_t GetIndex(const std::vector<std::size_t>& index) const override
    {
        KRATOS_ERROR << "GetIndex in 0D is meaningless";
    }

    /// Get the value at specific grid point
    TDataType& GetValue(std::size_t i)
    {
        if (i != 0)
            KRATOS_ERROR << "Value is only defined at position 0";
        return BaseType::Data()[i];
    }

    /// Get the value at specific grid point
    const TDataType& GetValue(std::size_t i) const
    {
        if (i != 0)
            KRATOS_ERROR << "Value is only defined at position 0";
        return BaseType::Data()[i];
    }

    /// Set the value at specific grid point
    void SetValue(std::size_t i, const TDataType& value)
    {
        if (i != 0)
            KRATOS_ERROR << "Value is only defined at position 0";
        BaseType::Data()[i] = value;
    }

    // overload operator ()
    TDataType& operator() (std::size_t i)
    {
        if (i != 0)
            KRATOS_ERROR << "Value is only defined at position 0";
        return BaseType::Data()[i];
    }

    // overload operator ()
    const TDataType& operator() (std::size_t i) const
    {
        if (i != 0)
            KRATOS_ERROR << "Value is only defined at position 0";
        return BaseType::Data()[i];
    }

    /// Copy the data the other grid. The size of two grids must be equal.
    void CopyFrom(const ControlGrid<TDataType>& rOther) override
    {
        BaseType::CopyFrom(rOther);
    }

    /// Copy the data the other grid. The size of two grids must be equal.
    void CopyFrom(const typename ControlGrid<TDataType>::Pointer pOther) override
    {
        BaseType::CopyFrom(pOther);
    }

    /// Copy the data the other grid. The size of two grids must be equal.
    virtual void CopyFrom(const StructuredControlGrid<0, TDataType>& rOther)
    {
        for (std::size_t i = 0; i < this->Size(); ++i)
        {
            this->SetValue(i, rOther.GetValue(i));
        }
    }

    /// Copy the data the other grid
    virtual void CopyFrom(const typename StructuredControlGrid<1, TDataType>::Pointer pOther)
    {
        this->CopyFrom(*pOther);
    }

    /// Copy the data the other grid. In the case that the source has different size, the grid is resized.
    virtual void ResizeAndCopyFrom(const StructuredControlGrid<0, TDataType>& rOther)
    {
        for (std::size_t i = 0; i < this->Size(); ++i)
        {
            this->SetValue(i, rOther.GetValue(i));
        }
    }

    /// Copy the data the other grid. In the case that the source has different size, the grid is resized.
    virtual void ResizeAndCopyFrom(const typename StructuredControlGrid<0, TDataType>::Pointer pOther)
    {
        this->ResizeAndCopyFrom(*pOther);
    }

    /// Reverse the control grid in specific dimension
    void Reverse(int idir) override
    {
        KRATOS_ERROR << "Reverse does not make sense for 0D grid";
    }

    /// Transpose the control grid
    void Transpose(int idir, int jdir) override
    {
        KRATOS_ERROR << "Tranpose does not make sense for 0D grid";
    }

    /// Create a connectivity matrix for the structured control grid
    void CreateConnectivity(std::size_t offset, std::vector<std::vector<std::size_t> >& connectivities) const override
    {
        connectivities.clear();
        connectivities.resize(1);
        connectivities[0].resize(1);
        connectivities[0][0] = 0;
    }

    /// Overload assignment operator
    StructuredControlGrid<0, TDataType>& operator=(const StructuredControlGrid<0, TDataType>& rOther)
    {
        BaseType::operator=(rOther);
        this->CopyFrom(rOther);
        return *this;
    }

    /// Clone this grid
    typename ControlGrid<TDataType>::Pointer Clone() const override
    {
        typename StructuredControlGrid<0, TDataType>::Pointer pNewControlGrid = typename StructuredControlGrid<0, TDataType>::Pointer(new StructuredControlGrid<0, TDataType>(1));
        *pNewControlGrid = *this;
        return pNewControlGrid;
    }

    std::string Type() const override
    {
        return StaticType();
    }

    /// Get the static type of the control grid
    static std::string StaticType()
    {
        return "StructuredControlGrid0D";
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "StructuredGrid<0> " << BaseType::Name() << "[1]";
    }

    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << " Data:\n (";
        for (std::size_t i = 0; i < BaseType::Data().size(); ++i)
        {
            rOStream << " " << BaseType::Data()[i];
        }
        rOStream << ")" << std::endl;
    }
};

template<typename TDataType>
class StructuredControlGrid<1, TDataType> : public BaseStructuredControlGrid<TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(StructuredControlGrid);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const StructuredControlGrid> ConstPointer;
#endif

    // type definitions
    typedef BaseStructuredControlGrid<TDataType> BaseType;
    typedef typename BaseType::DataContainerType DataContainerType;
    typedef typename BaseType::DataType DataType;

    /// Constructor with size
    StructuredControlGrid(const std::vector<std::size_t>& sizes) : BaseType(), mSize(sizes[0])
    {
        BaseType::Data().resize(sizes[0]);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType());
    }

    /// Constructor with size
    StructuredControlGrid(std::size_t n) : BaseType(), mSize(n)
    {
        BaseType::Data().resize(n);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType());
    }

    /// Copy constructor
    StructuredControlGrid(const StructuredControlGrid& rOther)
    : BaseType(rOther), mSize(rOther.mSize)
    {}

    /// Destructor
    ~StructuredControlGrid() override {}

    /// Create a new control grid pointer
    static typename StructuredControlGrid<1, TDataType>::Pointer Create(const std::vector<std::size_t>& sizes)
    {
        return typename StructuredControlGrid<1, TDataType>::Pointer(new StructuredControlGrid<1, TDataType>(sizes));
    }

    /// resize the grid
    void Resize(std::size_t new_size)
    {
        resize(new_size);
    }

    /// resize the grid
    void resize(std::size_t new_size)
    {
        mSize = new_size;
        BaseType::Data().resize(new_size);
    }

    /// Get the size of underlying data
    std::size_t Size() const {return BaseType::Data().size();}

    /// Get the size of the grid is specific dimension
    std::size_t Size(std::size_t dim) const {return mSize;}

    /// Get the global index of the grid value providing the index vector
    std::size_t GetIndex(const std::vector<std::size_t>& index) const override
    {
        return index[0];
    }

    /// Get the value at specific grid point
    TDataType& GetValue(std::size_t i)
    {
        return BaseType::Data()[i];
    }

    /// Get the value at specific grid point
    const TDataType& GetValue(std::size_t i) const
    {
        return BaseType::Data()[i];
    }

    /// Set the value at specific grid point
    void SetValue(std::size_t i, const TDataType& value)
    {
        BaseType::Data()[i] = value;
    }

    // overload operator ()
    TDataType& operator() (std::size_t i) {return BaseType::Data()[i];}

    // overload operator ()
    const TDataType& operator() (std::size_t i) const {return BaseType::Data()[i];}

    /// Copy the data the other grid. The size of two grids must be equal.
    void CopyFrom(const ControlGrid<TDataType>& rOther) override
    {
        BaseType::CopyFrom(rOther);
    }

    /// Copy the data the other grid. The size of two grids must be equal.
    void CopyFrom(const typename ControlGrid<TDataType>::Pointer pOther) override
    {
        BaseType::CopyFrom(pOther);
    }

    /// Copy the data the other grid. The size of two grids must be equal.
    virtual void CopyFrom(const StructuredControlGrid<1, TDataType>& rOther)
    {
        if (rOther.Size() != this->Size())
        {
            KRATOS_ERROR << "The size of the grid is incompatible"
                         << ", this size = " << this->Size()
                         << ", other size = " << rOther.Size();
        }

        for (std::size_t i = 0; i < this->Size(); ++i)
        {
            this->SetValue(i, rOther.GetValue(i));
        }
    }

    /// Copy the data the other grid
    virtual void CopyFrom(const typename StructuredControlGrid<1, TDataType>::Pointer pOther)
    {
        this->CopyFrom(*pOther);
    }

    /// Copy the data the other grid. In the case that the source has different size, the grid is resized.
    virtual void ResizeAndCopyFrom(const StructuredControlGrid<1, TDataType>& rOther)
    {
        if (rOther.Size() != this->Size())
        {
            BaseType::Data().resize(rOther.Size());
        }
        for (std::size_t i = 0; i < this->Size(); ++i)
        {
            this->SetValue(i, rOther.GetValue(i));
        }
    }

    /// Copy the data the other grid. In the case that the source has different size, the grid is resized.
    virtual void ResizeAndCopyFrom(const typename StructuredControlGrid<1, TDataType>::Pointer pOther)
    {
        this->ResizeAndCopyFrom(*pOther);
    }

    /// Reverse the control grid in specific dimension
    void Reverse(int idir) override
    {
        if (idir == 0)
        {
            std::reverse(BaseType::Data().begin(), BaseType::Data().end());
        }
        else
            KRATOS_ERROR << "Invalid direction " << idir;
    }

    /// Transpose the control grid
    void Transpose(int idir, int jdir) override
    {
        KRATOS_ERROR << "Tranpose does not make sense for 1D grid";
    }

    /// Create a connectivity matrix for the structured control grid
    void CreateConnectivity(std::size_t offset, std::vector<std::vector<std::size_t> >& connectivities) const override
    {
        connectivities.clear();
        connectivities.resize(this->Size() - 1);

        for (std::size_t i = 0; i < this->Size() - 1; ++i)
        {
            connectivities[i].resize(2);
            connectivities[i][0] = i + offset;
            connectivities[i][1] = i + 1 + offset;
        }
    }

    /// Overload assignment operator
    StructuredControlGrid<1, TDataType>& operator=(const StructuredControlGrid<1, TDataType>& rOther)
    {
        BaseType::operator=(rOther);
        this->mSize = rOther.mSize;
        this->CopyFrom(rOther);
        return *this;
    }

    /// Clone this grid
    typename ControlGrid<TDataType>::Pointer Clone() const override
    {
        typename StructuredControlGrid<1, TDataType>::Pointer pNewControlGrid = typename StructuredControlGrid<1, TDataType>::Pointer(new StructuredControlGrid<1, TDataType>(mSize));
        *pNewControlGrid = *this;
        return pNewControlGrid;
    }

    std::string Type() const override
    {
        return StaticType();
    }

    /// Get the static type of the control grid
    static std::string StaticType()
    {
        return "StructuredControlGrid1D";
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "StructuredGrid<1> " << BaseType::Name() << "[" << mSize << "]";
    }

    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << " Data:\n (";
        for (std::size_t i = 0; i < BaseType::Data().size(); ++i)
        {
            rOStream << " " << BaseType::Data()[i];
        }
        rOStream << ")" << std::endl;
    }

private:
    std::size_t mSize;

    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType);
        rSerializer.save("Size", mSize);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType);
        rSerializer.load("Size", mSize);
    }
    ///@}
};

template<typename TDataType>
class StructuredControlGrid<2, TDataType> : public BaseStructuredControlGrid<TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(StructuredControlGrid);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const StructuredControlGrid> ConstPointer;
#endif

    // type definitions
    typedef BaseStructuredControlGrid<TDataType> BaseType;
    typedef typename BaseType::DataContainerType DataContainerType;
    typedef typename BaseType::DataType DataType;

    /// Constructor with size
    StructuredControlGrid(const std::vector<std::size_t>& sizes) : BaseType(), mSize{sizes[0], sizes[1]}
    {
        BaseType::Data().resize(sizes[0]*sizes[1]);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType());
    }

    /// Constructor with size
    StructuredControlGrid(std::size_t m, std::size_t n) : mSize{m, n}
    {
        BaseType::Data().resize(m * n);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType());
    }

    /// Copy constructor
    StructuredControlGrid(const StructuredControlGrid& rOther)
    : BaseType(rOther)
    {
        mSize[0] = rOther.mSize[0];
        mSize[1] = rOther.mSize[1];
    }

    /// Destructor
    ~StructuredControlGrid() override {}

    /// Create a new control grid pointer
    static typename StructuredControlGrid<2, TDataType>::Pointer Create(const std::vector<std::size_t>& sizes)
    {
        return typename StructuredControlGrid<2, TDataType>::Pointer(new StructuredControlGrid<2, TDataType>(sizes));
    }

    /// resize the grid
    void Resize(std::size_t new_size1, std::size_t new_size2)
    {
        resize(new_size1, new_size2);
    }

    /// resize the grid
    void resize(std::size_t new_size1, std::size_t new_size2)
    {
        mSize[0] = new_size1;
        mSize[1] = new_size2;
        BaseType::Data().resize(new_size1 * new_size2);
    }

    /// Get the size of underlying data
    std::size_t Size() const {return BaseType::Data().size();}

    /// Get the size of the grid is specific dimension
    std::size_t Size(std::size_t dim) const {return mSize[dim];}

    /// Get the global index of the grid value providing the index vector
    std::size_t GetIndex(const std::vector<std::size_t>& index) const override
    {
        return index[1] * mSize[0] + index[0];
    }

    /// Get the value at specific grid point
    TDataType& GetValue(std::size_t i, std::size_t j)
    {
        return BaseType::Data()[j * mSize[0] + i];
    }

    /// Get the value at specific grid point
    const TDataType& GetValue(std::size_t i, std::size_t j) const
    {
        return BaseType::Data()[j * mSize[0] + i];
    }

    /// Set the value at specific grid point
    void SetValue(std::size_t i, std::size_t j, const TDataType& value)
    {
        BaseType::Data()[j * mSize[0] + i] = value;
    }

    // overload operator ()
    TDataType& operator() (std::size_t i, std::size_t j) {return BaseType::Data()[j * mSize[0] + i];}

    // overload operator ()
    const TDataType& operator() (std::size_t i, std::size_t j) const {return BaseType::Data()[j * mSize[0] + i];}

    /// Copy the data the other grid. The size of two grids must be equal.
    void CopyFrom(const ControlGrid<TDataType>& rOther) override
    {
        BaseType::CopyFrom(rOther);
    }

    /// Copy the data the other grid. The size of two grids must be equal.
    void CopyFrom(const typename ControlGrid<TDataType>::Pointer pOther) override
    {
        BaseType::CopyFrom(pOther);
    }

    /// Copy the data the other grid
    virtual void CopyFrom(const StructuredControlGrid<2, TDataType>& rOther)
    {
        if ( ( rOther.Size(0) != this->Size(0) ) || ( rOther.Size(1) != this->Size(1) ) )
        {
            KRATOS_ERROR << "The size of the grid is incompatible"
                         << ", this size = (" << this->Size(0) << ", " << this->Size(1) << ")"
                         << ", other size = (" << rOther.Size(0) << ", " << rOther.Size(1) << ")";
        }

        for (std::size_t i = 0; i < this->Size(0); ++i)
        {
            for (std::size_t j = 0; j < this->Size(1); ++j)
            {
                this->SetValue(i, j, rOther.GetValue(i, j));
            }
        }
    }

    /// Copy the data the other grid
    virtual void CopyFrom(const typename StructuredControlGrid<2, TDataType>::Pointer pOther)
    {
        this->CopyFrom(*pOther);
    }

    /// Copy the data from a column of structured control_grid in 1D
    /// If dir==0, the column of grid will be copied along u-direction
    /// If dir==1, the column of grid will be copied along v-direction
    void CopyFrom(int dir, const std::vector<typename StructuredControlGrid<1, TDataType>::Pointer>& pOthers)
    {
        if (pOthers.size() != this->Size(dir))
        {
            KRATOS_ERROR << "The size is incompatible";
        }

        if (dir == 0)
        {
            for (std::size_t i = 0; i < this->Size(0); ++i)
            {
                for (std::size_t j = 0; j < this->Size(1); ++j)
                {
                    this->SetValue(i, j, pOthers[i]->GetValue(j));
                }
            }
        }
        else if (dir == 1)
        {
            for (std::size_t i = 0; i < this->Size(0); ++i)
            {
                for (std::size_t j = 0; j < this->Size(1); ++j)
                {
                    this->SetValue(i, j, pOthers[j]->GetValue(i));
                }
            }
        }
    }

    /// Copy the data the other grid. In the case that the source has different size, the grid is resized.
    virtual void ResizeAndCopyFrom(const StructuredControlGrid<2, TDataType>& rOther)
    {
        if ( ( rOther.Size(0) != this->Size(1) ) || ( rOther.Size(1) != this->Size(1) ) )
        {
            BaseType::Data().resize(rOther.Size(0)*rOther.Size(1));
        }
        for (std::size_t i = 0; i < this->Size(0); ++i)
        {
            for (std::size_t j = 0; j < this->Size(1); ++j)
            {
                this->SetValue(i, j, rOther.GetValue(i, j));
            }
        }
    }

    /// Copy the data the other grid. In the case that the source has different size, the grid is resized.
    virtual void ResizeAndCopyFrom(const typename StructuredControlGrid<2, TDataType>::Pointer pOther)
    {
        this->ResizeAndCopyFrom(*pOther);
    }

    /// Reverse the control grid in specific dimension
    void Reverse(int idir) override
    {
        BSplinesIndexingUtility::Reverse<2, DataContainerType, std::size_t*>(BaseType::Data(), mSize, idir);
    }

    /// Transpose the control grid
    void Transpose(int idir, int jdir) override
    {
        BSplinesIndexingUtility::Transpose<2, DataContainerType, std::size_t*>(BaseType::Data(), mSize, idir, jdir);

        // swap size
        auto isize = mSize[idir];
        auto jsize = mSize[jdir];
        mSize[idir] = jsize;
        mSize[jdir] = isize;
    }

    /// Get the layer of control grid from the boundary, if the level = 0, the control grid on the boundary will be extracted.
    typename StructuredControlGrid<1, TDataType>::Pointer Get(const BoundarySide side, std::size_t level)
    {
        typename StructuredControlGrid<1, TDataType>::Pointer pControlGrid;

        if (side == _BLEFT_)
        {
            pControlGrid = typename StructuredControlGrid<1, TDataType>::Pointer( new StructuredControlGrid<1, TDataType>(this->Size(1)) );

            for (std::size_t i = 0; i < this->Size(1); ++i)
            {
                pControlGrid->SetValue(i, this->GetValue(0 + level, i));
            }
        }
        else if (side == _BRIGHT_)
        {
            pControlGrid = typename StructuredControlGrid<1, TDataType>::Pointer( new StructuredControlGrid<1, TDataType>(this->Size(1)) );

            for (std::size_t i = 0; i < this->Size(1); ++i)
            {
                pControlGrid->SetValue(i, this->GetValue(this->Size(0) - 1 - level, i));
            }
        }
        else if (side == _BTOP_)
        {
            pControlGrid = typename StructuredControlGrid<1, TDataType>::Pointer( new StructuredControlGrid<1, TDataType>(this->Size(0)) );

            for (std::size_t i = 0; i < this->Size(0); ++i)
            {
                pControlGrid->SetValue(i, this->GetValue(i, this->Size(1) - 1 - level));
            }
        }
        else if (side == _BBOTTOM_)
        {
            pControlGrid = typename StructuredControlGrid<1, TDataType>::Pointer( new StructuredControlGrid<1, TDataType>(this->Size(0)) );

            for (std::size_t i = 0; i < this->Size(0); ++i)
            {
                pControlGrid->SetValue(i, this->GetValue(i, 0 + level));
            }
        }
        else
        {
            KRATOS_ERROR << "Invalid side " << side;
        }

        return pControlGrid;
    }

    /// Create a connectivity for the structured control grid
    void CreateConnectivity(std::size_t offset, std::vector<std::vector<std::size_t> >& connectivities) const override
    {
        connectivities.clear();
        connectivities.resize((this->Size(0) - 1) * (this->Size(1) - 1));

        std::size_t cnt = 0;
        for (std::size_t i = 0; i < this->Size(0) - 1; ++i)
        {
            for (std::size_t j = 0; j < this->Size(1) - 1; ++j)
            {
                connectivities[cnt].resize(4);
                connectivities[cnt][0] = i     + j * this->Size(0)       + offset;
                connectivities[cnt][1] = i + 1 + j * this->Size(0)       + offset;
                connectivities[cnt][2] = i + 1 + (j + 1) * this->Size(0) + offset;
                connectivities[cnt][3] = i     + (j + 1) * this->Size(0) + offset;
                ++cnt;
            }
        }
    }

    /// Overload assignment operator
    StructuredControlGrid<2, TDataType>& operator=(const StructuredControlGrid<2, TDataType>& rOther)
    {
        BaseType::operator=(rOther);
        this->mSize[0] = rOther.mSize[0];
        this->mSize[1] = rOther.mSize[1];
        this->CopyFrom(rOther);
        return *this;
    }

    /// Clone this grid
    typename ControlGrid<TDataType>::Pointer Clone() const override
    {
        typename StructuredControlGrid<2, TDataType>::Pointer pNewControlGrid = typename StructuredControlGrid<2, TDataType>::Pointer(new StructuredControlGrid<2, TDataType>(mSize[0], mSize[1]));
        *pNewControlGrid = *this;
        return pNewControlGrid;
    }

    std::string Type() const override
    {
        return StaticType();
    }

    /// Get the static type of the control grid
    static std::string StaticType()
    {
        return "StructuredControlGrid2D";
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "StructuredGrid<2> " << BaseType::Name() << ": [" << mSize[0] << ", " << mSize[1] << "]";
    }

    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << " Data:\n (\n";
        for (std::size_t j = 0; j < mSize[1]; ++j)
        {
            rOStream << "  (";
            {
                for (std::size_t i = 0; i < mSize[0]; ++i)
                {
                    rOStream << " " << GetValue(i, j);
                }
            }
            rOStream << ")" << std::endl;
        }
        rOStream << " )" << std::endl;
    }

private:
    std::size_t mSize[2];

    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType);
        rSerializer.save("Size1", mSize[0]);
        rSerializer.save("Size2", mSize[1]);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType);
        rSerializer.load("Size1", mSize[0]);
        rSerializer.load("Size2", mSize[1]);
    }
    ///@}
};

template<typename TDataType>
class StructuredControlGrid<3, TDataType> : public BaseStructuredControlGrid<TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(StructuredControlGrid);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const StructuredControlGrid> ConstPointer;
#endif

    // type definitions
    typedef BaseStructuredControlGrid<TDataType> BaseType;
    typedef typename BaseType::DataContainerType DataContainerType;
    typedef typename BaseType::DataType DataType;

    /// Constructor with size
    StructuredControlGrid(const std::vector<std::size_t>& sizes) : BaseType(), mSize{sizes[0], sizes[1], sizes[2]}
    {
        BaseType::Data().resize(sizes[0]*sizes[1]*sizes[2]);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType());
    }

    /// Constructor with size
    StructuredControlGrid(std::size_t m, std::size_t n, std::size_t p) : mSize{m, n, p}
    {
        BaseType::Data().resize(m * n * p);
        std::fill(BaseType::Data().begin(), BaseType::Data().end(), TDataType());
    }

    /// Copy constructor
    StructuredControlGrid(const StructuredControlGrid& rOther)
    : BaseType(rOther)
    {
        mSize[0] = rOther.mSize[0];
        mSize[1] = rOther.mSize[1];
        mSize[2] = rOther.mSize[2];
    }

    /// Destructor
    ~StructuredControlGrid() override {}

    /// Create a new control grid pointer
    static typename StructuredControlGrid<3, TDataType>::Pointer Create(const std::vector<std::size_t>& sizes)
    {
        return typename StructuredControlGrid<3, TDataType>::Pointer(new StructuredControlGrid<3, TDataType>(sizes));
    }

    /// resize the grid
    void Resize(std::size_t new_size1, std::size_t new_size2, std::size_t new_size3)
    {
        resize(new_size1, new_size2, new_size3);
    }

    /// resize the grid
    void resize(std::size_t new_size1, std::size_t new_size2, std::size_t new_size3)
    {
        mSize[0] = new_size1;
        mSize[1] = new_size2;
        mSize[2] = new_size3;
        BaseType::Data().resize(new_size1 * new_size2 * new_size3);
    }

    /// Get the size of underlying data
    std::size_t Size() const {return BaseType::Data().size();}

    /// Get the size of the grid is specific dimension
    std::size_t Size(std::size_t dim) const {return mSize[dim];}

    /// Get the global index of the grid value providing the index vector
    std::size_t GetIndex(const std::vector<std::size_t>& index) const override
    {
        return (index[2] * mSize[1] + index[1]) * mSize[0] + index[0];
    }

    /// Get the value at specific grid point
    TDataType& GetValue(std::size_t i, std::size_t j, std::size_t k)
    {
        return BaseType::Data()[(k * mSize[1] + j) * mSize[0] + i];
    }

    /// Get the value at specific grid point
    const TDataType& GetValue(std::size_t i, std::size_t j, std::size_t k) const
    {
        return BaseType::Data()[(k * mSize[1] + j) * mSize[0] + i];
    }

    /// Set the value at specific grid point
    void SetValue(std::size_t i, std::size_t j, std::size_t k, const TDataType& value)
    {
        BaseType::Data()[(k * mSize[1] + j)*mSize[0] + i] = value;
    }

    // overload operator ()
    TDataType& operator() (std::size_t i, std::size_t j, std::size_t k)
    {
        return BaseType::Data()[(k * mSize[1] + j) * mSize[0] + i];
    }

    // overload operator ()
    const TDataType& operator() (std::size_t i, std::size_t j, std::size_t k) const
    {
        return BaseType::Data()[(k * mSize[1] + j) * mSize[0] + i];
    }

    /// Copy the data the other grid. The size of two grids must be equal.
    void CopyFrom(const ControlGrid<TDataType>& rOther) override
    {
        BaseType::CopyFrom(rOther);
    }

    /// Copy the data the other grid. The size of two grids must be equal.
    void CopyFrom(const typename ControlGrid<TDataType>::Pointer pOther) override
    {
        BaseType::CopyFrom(pOther);
    }

    /// Copy the data the other grid
    virtual void CopyFrom(const StructuredControlGrid<3, TDataType>& rOther)
    {
        if ( ( rOther.Size(0) != this->Size(0) )
                || ( rOther.Size(1) != this->Size(1) )
                || ( rOther.Size(2) != this->Size(2) ) )
        {
            KRATOS_ERROR << "The size of the grid is incompatible"
                         << ", this size = (" << this->Size(0) << ", " << this->Size(1) << ", " << this->Size(2) << ")"
                         << ", other size = (" << rOther.Size(0) << ", " << rOther.Size(1) << ", " << rOther.Size(2) << ")";
        }

        for (std::size_t i = 0; i < this->Size(0); ++i)
        {
            for (std::size_t j = 0; j < this->Size(1); ++j)
            {
                for (std::size_t k = 0; k < this->Size(2); ++k)
                {
                    this->SetValue(i, j, k, rOther.GetValue(i, j, k));
                }
            }
        }
    }

    /// Copy the data the other grid
    virtual void CopyFrom(const typename StructuredControlGrid<3, TDataType>::Pointer pOther)
    {
        this->CopyFrom(*pOther);
    }

    /// Copy the data from a column of structured control_grid in 2D
    /// If dir==0, the column of grid will be copied along u-direction
    /// If dir==1, the column of grid will be copied along v-direction
    /// If dir==2, the column of grid will be copied along w-direction
    void CopyFrom(int dir, const std::vector<typename StructuredControlGrid<2, TDataType>::Pointer>& pOthers)
    {
        if (pOthers.size() != this->Size(dir))
        {
            KRATOS_ERROR << "The size is incompatible";
        }

        if (dir == 0)
        {
            for (std::size_t i = 0; i < this->Size(0); ++i)
                for (std::size_t j = 0; j < this->Size(1); ++j)
                    for (std::size_t k = 0; k < this->Size(2); ++k)
                    {
                        this->SetValue(i, j, k, pOthers[i]->GetValue(j, k));
                    }
        }
        else if (dir == 1)
        {
            for (std::size_t i = 0; i < this->Size(0); ++i)
                for (std::size_t j = 0; j < this->Size(1); ++j)
                    for (std::size_t k = 0; k < this->Size(2); ++k)
                    {
                        this->SetValue(i, j, k, pOthers[j]->GetValue(i, k));
                    }
        }
        else if (dir == 2)
        {
            for (std::size_t i = 0; i < this->Size(0); ++i)
                for (std::size_t j = 0; j < this->Size(1); ++j)
                    for (std::size_t k = 0; k < this->Size(2); ++k)
                    {
                        this->SetValue(i, j, k, pOthers[k]->GetValue(i, j));
                    }
        }
    }

    /// Copy the data the other grid. In the case that the source has different size, the grid is resized.
    virtual void ResizeAndCopyFrom(const StructuredControlGrid<3, TDataType>& rOther)
    {
        if ( ( rOther.Size(0) != this->Size(1) )
                || ( rOther.Size(1) != this->Size(1) )
                || ( rOther.Size(2) != this->Size(2) ) )
        {
            BaseType::Data().resize(rOther.Size(0)*rOther.Size(1)*rOther.Size(2));
        }
        for (std::size_t i = 0; i < this->Size(0); ++i)
            for (std::size_t j = 0; j < this->Size(1); ++j)
                for (std::size_t k = 0; k < this->Size(2); ++k)
                {
                    this->SetValue(i, j, k, rOther.GetValue(i, j, k));
                }
    }

    /// Copy the data the other grid. In the case that the source has different size, the grid is resized.
    virtual void ResizeAndCopyFrom(const typename StructuredControlGrid<3, TDataType>::Pointer pOther)
    {
        this->ResizeAndCopyFrom(*pOther);
    }

    /// Reverse the control grid in specific dimension
    void Reverse(int idir) override
    {
        BSplinesIndexingUtility::Reverse<3, DataContainerType, std::size_t*>(BaseType::Data(), mSize, idir);
    }

    /// Transpose the control grid
    void Transpose(int idir, int jdir) override
    {
        BSplinesIndexingUtility::Transpose<3, DataContainerType, std::size_t*>(BaseType::Data(), mSize, idir, jdir);

        // swap size
        auto isize = mSize[idir];
        auto jsize = mSize[jdir];
        mSize[idir] = jsize;
        mSize[jdir] = isize;
    }

    /// Get the layer of control grid from the boundary, if the level = 0, the control grid on the boundary will be extracted.
    typename StructuredControlGrid<2, TDataType>::Pointer Get(const BoundarySide side, const std::size_t level) const
    {
        // TODO
        KRATOS_ERROR << "Not yet implemented";

        typename StructuredControlGrid<2, TDataType>::Pointer pControlGrid;

        if (side == _BLEFT_)
        {
        }
        else if (side == _BRIGHT_)
        {
        }
        else if (side == _BTOP_)
        {
        }
        else if (side == _BBOTTOM_)
        {
        }
        else if (side == _BFRONT_)
        {
        }
        else if (side == _BBACK_)
        {
        }

        return pControlGrid;
    }

    /// Create a connectivity for the structured control grid
    void CreateConnectivity(std::size_t offset, std::vector<std::vector<std::size_t> >& connectivities) const override
    {
        connectivities.clear();
        connectivities.resize((this->Size(0) - 1) * (this->Size(1) - 1) * (this->Size(2) - 1));

        std::size_t cnt = 0;
        for (std::size_t i = 0; i < this->Size(0) - 1; ++i)
        {
            for (std::size_t j = 0; j < this->Size(1) - 1; ++j)
            {
                for (std::size_t k = 0; k < this->Size(2) - 1; ++k)
                {
                    connectivities[cnt].resize(8);
                    connectivities[cnt][0] = i + j * this->Size(0) + k * this->Size(0) * this->Size(1) + offset;
                    connectivities[cnt][1] = i + 1 + j * this->Size(0) + k * this->Size(0) * this->Size(1) + offset;
                    connectivities[cnt][2] = i + 1 + (j + 1) * this->Size(0) + k * this->Size(0) * this->Size(1) + offset;
                    connectivities[cnt][3] = i + (j + 1) * this->Size(0) + k * this->Size(0) * this->Size(1) + offset;
                    connectivities[cnt][4] = connectivities[cnt][0] + this->Size(0) * this->Size(1);
                    connectivities[cnt][5] = connectivities[cnt][1] + this->Size(0) * this->Size(1);
                    connectivities[cnt][6] = connectivities[cnt][2] + this->Size(0) * this->Size(1);
                    connectivities[cnt][7] = connectivities[cnt][3] + this->Size(0) * this->Size(1);
                    ++cnt;
                }
            }
        }
    }

    /// Overload assignment operator
    StructuredControlGrid<3, TDataType>& operator=(const StructuredControlGrid<3, TDataType>& rOther)
    {
        BaseType::operator=(rOther);
        this->mSize[0] = rOther.mSize[0];
        this->mSize[1] = rOther.mSize[1];
        this->mSize[2] = rOther.mSize[2];
        this->CopyFrom(rOther);
        return *this;
    }

    /// Clone this grid
    typename ControlGrid<TDataType>::Pointer Clone() const override
    {
        typename StructuredControlGrid<3, TDataType>::Pointer pNewControlGrid = typename StructuredControlGrid<3, TDataType>::Pointer(new StructuredControlGrid<3, TDataType>(mSize[0], mSize[1], mSize[2]));
        *pNewControlGrid = *this;
        return pNewControlGrid;
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "StructuredGrid<3> " << BaseType::Name() << ": [" << mSize[0] << ", " << mSize[1] << ", " << mSize[2] << "]";
    }

    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << " Data:\n (";
        for (std::size_t k = 0; k < mSize[2]; ++k)
        {
            rOStream << " (";
            for (std::size_t j = 0; j < mSize[1]; ++j)
            {
                rOStream << " (";
                for (std::size_t i = 0; i < mSize[0]; ++i)
                {
                    rOStream << " " << GetValue(i, j, k);
                }
            }
            rOStream << ")" << std::endl;
        }
        rOStream << " )" << std::endl;
    }

    std::string Type() const override
    {
        return StaticType();
    }

    /// Get the static type of the control grid
    static std::string StaticType()
    {
        return "StructuredControlGrid3D";
    }

private:
    std::size_t mSize[3];

    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType);
        rSerializer.save("Size1", mSize[0]);
        rSerializer.save("Size2", mSize[1]);
        rSerializer.save("Size3", mSize[2]);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType);
        rSerializer.load("Size1", mSize[0]);
        rSerializer.load("Size2", mSize[1]);
        rSerializer.load("Size3", mSize[2]);
    }
    ///@}
};

/// output stream function
template<int TDim, typename TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const StructuredControlGrid<TDim, TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_REGULAR_CONTROL_GRID_H_INCLUDED defined
