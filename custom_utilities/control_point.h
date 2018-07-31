//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 4 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_POINT_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_POINT_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "custom_utilities/trans/transformation.h"
#include "custom_utilities/control_value.h"

namespace Kratos
{

/**
    Represent a control point in isogeometric mesh topology.
 */
template<typename TDataType>
class ControlPoint : public ControlValue<array_1d<TDataType, 3>, TDataType>
{
public:
    // Type definitions
    typedef ControlValue<array_1d<TDataType, 3>, TDataType> BaseType;
    typedef typename BaseType::DataType CoordinatesType;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ControlPoint);

    /// Default constructor
    ControlPoint() : BaseType() {}

    /// Constant constructor
    ControlPoint(const double& v)
    {
        BaseType::WV()[0] = v;
        BaseType::WV()[1] = v;
        BaseType::WV()[2] = v;
        BaseType::W() = v;
    }

    /// Constructor with full coordinates
    ControlPoint(const double& wx, const double& wy, const double& wz, const double& w)
    {
        BaseType::WV()[0] = wx;
        BaseType::WV()[1] = wy;
        BaseType::WV()[2] = wz;
        BaseType::W() = w;
    }

    /// Destructor
    virtual ~ControlPoint() {}

    /// homogeneous X-coordinate
    TDataType& WX() {return BaseType::WV()[0];}
    const TDataType& WX() const {return BaseType::WV()[0];}

    /// X-coordinate
    TDataType X() const {return BaseType::WV()[0] / BaseType::W();}

    /// homogeneous Y-coordinate
    TDataType& WY() {return BaseType::WV()[1];}
    const TDataType& WY() const {return BaseType::WV()[1];}

    /// Y-coordinate
    TDataType Y() const {return BaseType::WV()[1] / BaseType::W();}

    /// homogeneous Z-coordinate
    TDataType& WZ() {return BaseType::WV()[2];}
    const TDataType& WZ() const {return BaseType::WV()[2];}

    /// Z-coordinate
    TDataType Z() const {return BaseType::WV()[2] / BaseType::W();}

    /// Set the coordinate. The input is the physical coordinates in 3D space.
    void SetCoordinates(const TDataType& _X, const TDataType& _Y, const TDataType& _Z, const TDataType& _W)
    {
        BaseType::WV()[0] = _W*_X;
        BaseType::WV()[1] = _W*_Y;
        BaseType::WV()[2] = _W*_Z;
        BaseType::W()  = _W;
    }

    /// Add to the coordinate. The input is the increment of physical coordinates in 3D space.
    void AddCoordinates(const TDataType& _X, const TDataType& _Y, const TDataType& _Z, const TDataType& _W)
    {
        BaseType::WV()[0] += _W*_X;
        BaseType::WV()[1] += _W*_Y;
        BaseType::WV()[2] += _W*_Z;
        BaseType::W()  += _W;
    }

    // overload operator []
    TDataType& operator[] (const int& i)
    {
        if (i == 0) return WX();
        else if (i == 1) return WY();
        else if (i == 2) return WZ();
        else if (i == 3) return BaseType::W();
    }

    const TDataType& operator[] (const int& i) const
    {
        if (i == 0) return WX();
        else if (i == 1) return WY();
        else if (i == 2) return WZ();
        else if (i == 3) return BaseType::W();
    }

    // overload operator ()
    TDataType operator() (const int& i) const
    {
        if (i == 0) return X();
        else if (i == 1) return Y();
        else if (i == 2) return Z();
        else if (i == 3) return BaseType::W();
    }

    /// Assignment operator
    ControlPoint& operator=(const BaseType& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    /// Addition operator
    ControlPoint& operator+=(const BaseType& rOther)
    {
        BaseType::operator+=(rOther);
        return *this;
    }

    /// Addition operator
    friend ControlPoint operator+(ControlPoint lhs, const ControlPoint& rhs)
    {
        lhs += rhs;
        return lhs;
    }

    /// Multiplication operator
    ControlPoint& operator*=(const TDataType& alpha)
    {
        BaseType::operator*=(alpha);
        return *this;
    }

    /// Multiplication operator
    friend ControlPoint operator*(const TDataType& alpha, BaseType c)
    {
        ControlPoint p;
        p = alpha*c;
        return p;
    }

    /// Multiplication operator
    friend ControlPoint operator*(BaseType c, const TDataType& alpha)
    {
        ControlPoint p;
        p = alpha*c;
        return p;
    }

    /// Apply the homogeneous transformation to the control point
    void ApplyTransformation(const Transformation<TDataType>& trans)
    {
        TDataType res[4];
        for (std::size_t i = 0; i < 4; ++i)
            res[i] = trans(i, 0)*BaseType::WV()[0]
                   + trans(i, 1)*BaseType::WV()[1]
                   + trans(i, 2)*BaseType::WV()[2]
                   + trans(i, 3)*BaseType::W();
        BaseType::WV()[0] = res[0];
        BaseType::WV()[1] = res[1];
        BaseType::WV()[2] = res[2];
        BaseType::W() = res[3];
    }

    /// Multiplication operator
    friend ControlPoint operator*(const Transformation<TDataType>& trans, ControlPoint c)
    {
        c.ApplyTransformation(trans);
        return c;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Control Point";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
        // print the control point in homogeneous coordinates
        // rOStream << "(X: " << X() << ", Y: " << Y() << ", Z: " << Z() << ", W: " << BaseType::W() << ")";
        // rOStream << "(" << X() << ", " << Y() << ", " << Z() << ", " << BaseType::W() << ")";
        rOStream << "(" << WX() << ", " << WY() << ", " << WZ() << ", " << BaseType::W() << ")";
    }

private:

    /// Serializer
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
    }
};

/// output stream function
template<typename TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const ControlPoint<TDataType>& rThis)
{
    // rThis.PrintInfo(rOStream);
    // rOStream << ": ";
    rThis.PrintData(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_POINT_H_INCLUDED

