//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Apr 2015 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_KNOT_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_KNOT_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"

namespace Kratos
{

/**
Represent a knot in isogeometric mesh topology.
A knot is determined by its value and its index in the knot vectors. The knot value is fixed when knot is constructed, but the index can be changed if new knot are added to the knot vector
The reason to keep the two stage definition is to control the repeated knot values. In case when the two knots have the same value, they will be differentiated by indexing
 */
template<typename TDataType>
class Knot
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Knot);
#ifdef SD_APP_FORWARD_COMPATIBILITY
    typedef Kratos::shared_ptr<const Knot> ConstPointer;
#endif

    /// Empty constructor
    Knot() : mIndex(-1), mValue(0), mIsActive(true)
    {}

    /// Default constructor
    Knot(const TDataType& Value) : mIndex(-1), mValue(Value), mIsActive(true)
    {}

    /// Copy constructor
    template<typename TOtherDataType>
    Knot(const Knot<TOtherDataType>& rOther)
    : mIndex(rOther.mIndex), mValue(rOther.mValue), mIsActive(rOther.mIsActive)
    {}

    /// Clone this knot
    Knot::Pointer Clone() const
    {
        Knot::Pointer pKnot = Knot::Pointer(new Knot(*this));
        return pKnot;
    }

    /// Get and Set for knot index
    std::size_t Index() const {return mIndex;}
    void UpdateIndex(std::size_t Index) {mIndex = Index;}

    /// Get the knot value
    TDataType& Value() {return mValue;}
    const TDataType& Value() const {return mValue;}

    /// Get and Set for mIsActive
    bool IsActive() const {return mIsActive;}
    void SetActive(const bool IsActive) {mIsActive = IsActive;}

    /// Overload operator assignment
    template<typename TOtherDataType>
    Knot& operator=(const Knot<TOtherDataType>& rOther)
    {
        this->mValue = rOther.mValue;
        this->mIndex = rOther.mIndex;
        this->mIsActive = rOther.mIsActive;
        return *this;
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "(" << Index() << ", " << Value() << ")";
    }

private:

    std::size_t mIndex;
    TDataType mValue;
    bool mIsActive;

    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        rSerializer.save( "Index", mIndex );
        rSerializer.save( "Value", mValue );
        rSerializer.save( "IsActive", mIsActive );
    }

    void load(Serializer& rSerializer)
    {
        rSerializer.load( "Index", mIndex );
        rSerializer.load( "Value", mValue );
        rSerializer.load( "IsActive", mIsActive );
    }
    ///@}
};

/// output stream function
template<class TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const Knot<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_KNOT_H_INCLUDED
