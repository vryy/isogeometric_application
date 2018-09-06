//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 11 May 2018 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_VALUE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_VALUE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"

namespace Kratos
{

/**
    Represent a control value in isogeometric mesh topology.
 */
template<typename TDataType, typename TWeightType>
class ControlValue
{
public:
    /// Type definitions
    typedef TDataType DataType;
    typedef TWeightType WeightType;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ControlValue);

    /// Default constructor
    ControlValue() : mWV(TDataType(0.0)), mW(0.0) {}

    /// Constructor with full coordinates
    ControlValue(const TDataType& wv, const TWeightType& w) : mWV(wv), mW(w) {}

    /// Destructor
    virtual ~ControlValue() {}

    /// Homogeneous value
    TDataType& WV() {return mWV;}
    const TDataType& WV() const {return mWV;}

    /// Inhomogeneous value
    TDataType V() const {return mWV/mW;}

    /// Weight
    TWeightType& W() {return mW;}
    const TWeightType& W() const {return mW;}

    /// Set the coordinate. The input is the physical coordinates in 3D space.
    void Set(const TDataType& _V, const TDataType& _W)
    {
        mWV = _W*_V;
        mW  = _W;
    }

    /// Add to the coordinate. The input is the increment of physical coordinates in 3D space.
    void Add(const TDataType& _V, const TDataType& _W)
    {
        mWV += _W*_V;
        mW  += _W;
    }

    /// Assignment operator
    ControlValue& operator=(const ControlValue& rOther)
    {
        this->mWV = rOther.mWV;
        this->mW = rOther.mW;
        return *this;
    }

    /// Addition operator
    ControlValue& operator+=(const ControlValue& rOther)
    {
        this->mWV += rOther.mWV;
        this->mW += rOther.mW;
        return *this;
    }

    /// Addition operator
    friend ControlValue operator+(ControlValue lhs, const ControlValue& rhs)
    {
        lhs += rhs;
        return lhs;
    }

    /// Multiplication operator
    template<typename TValueType>
    ControlValue& operator*=(const TValueType& alpha)
    {
        this->mWV *= alpha;
        this->mW *= alpha;
        return *this;
    }

    /// Multiplication operator
    template<typename TValueType>
    friend ControlValue operator*(const TValueType& alpha, ControlValue c)
    {
        c *= alpha;
        return c;
    }

    /// Multiplication operator
    template<typename TValueType>
    friend ControlValue operator*(ControlValue c, const TValueType& alpha)
    {
        c *= alpha;
        return c;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Control Value";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        // print the control point in homogeneous coordinates
        rOStream << "(V: " << V() << ", W: " << W() << ")";
    }

private:

    TDataType mWV;
    TWeightType mW;

    /// Serializer
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        rSerializer.save( "WV", mWV );
        rSerializer.save( "W", mW );
    }

    virtual void load(Serializer& rSerializer)
    {
        rSerializer.load( "WV", mWV );
        rSerializer.load( "W", mW );
    }
};

/// output stream function
template<typename TDataType, typename TWeightType>
inline std::ostream& operator <<(std::ostream& rOStream, const ControlValue<TDataType, TWeightType>& rThis)
{
    // rThis.PrintInfo(rOStream);
    // rOStream << ": ";
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_VALUE_H_INCLUDED

