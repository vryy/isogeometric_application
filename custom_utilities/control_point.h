//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 4 Nov 2017 $
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
#include "brep_application/custom_algebra/trans/transformation.h"
#include "custom_utilities/control_value.h"

namespace Kratos
{

/**
 * Represent a control point in isogeometric mesh topology.
 */
template<typename TDataType, typename TWeightType = TDataType>
class ControlPoint : public ControlValue<array_1d<TDataType, 3>, TWeightType>
{
public:
    // Type definitions
    typedef TWeightType WeightType;
    typedef ControlValue<array_1d<TDataType, 3>, TWeightType> BaseType;
    typedef typename BaseType::DataType CoordinatesType;
    typedef TDataType CoordinateType;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ControlPoint);

    /// Default constructor
    ControlPoint() : BaseType() {}

    /// Constant constructor
    ControlPoint(TDataType v)
    {
        BaseType::WV()[0] = v;
        BaseType::WV()[1] = v;
        BaseType::WV()[2] = v;
        BaseType::W() = static_cast<TWeightType>(v);
    }

    /// Constructor with full coordinates
    ControlPoint(TDataType wx, TDataType wy, TDataType wz, TWeightType w)
    {
        BaseType::WV()[0] = wx;
        BaseType::WV()[1] = wy;
        BaseType::WV()[2] = wz;
        BaseType::W() = w;
    }

    /// Copy constructor
    ControlPoint(const ControlPoint& rOther) : BaseType(rOther) {}

    /// Destructor
    ~ControlPoint() override {}

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
    void SetCoordinates(const TDataType& _X, const TDataType& _Y, const TDataType& _Z, const TWeightType& _W)
    {
        BaseType::WV()[0] = _W * _X;
        BaseType::WV()[1] = _W * _Y;
        BaseType::WV()[2] = _W * _Z;
        BaseType::W()  = _W;
    }

    /// Add to the coordinate. The input is the increment of physical coordinates in 3D space.
    void AddCoordinates(const TDataType& _X, const TDataType& _Y, const TDataType& _Z, const TWeightType& _W)
    {
        BaseType::WV()[0] += _W * _X;
        BaseType::WV()[1] += _W * _Y;
        BaseType::WV()[2] += _W * _Z;
        BaseType::W()  += _W;
    }

    /// Apply the homogeneous transformation to the control point
    void ApplyTransformation(const Transformation<TDataType>& trans)
    {
        TDataType res[4];
        for (std::size_t i = 0; i < 4; ++i)
            res[i] = trans(i, 0) * BaseType::WV()[0]
                     + trans(i, 1) * BaseType::WV()[1]
                     + trans(i, 2) * BaseType::WV()[2]
                     + trans(i, 3) * BaseType::W();
        BaseType::WV()[0] = res[0];
        BaseType::WV()[1] = res[1];
        BaseType::WV()[2] = res[2];
        if constexpr (std::is_same<TWeightType, TDataType>::value)
        {
            BaseType::W() = res[3];
        }
        else if constexpr (std::is_same<TDataType, KRATOS_COMPLEX_TYPE>::value
                        && std::is_same<TWeightType, KRATOS_DOUBLE_TYPE>::value)
        {
            if (res[3].imag() == 0)
                BaseType::W() = res[3].real();
            else
                // here we shall throw and error because the transformation leads to complex weight.
                // Which is not meaningful because the weight is assumed to be double.
                // The typical transformation like translation, rotation shall preserve the data type of weight.
                KRATOS_ERROR << "Encounted complex weight after transformation. This is not allowed.";
        }
        else
            KRATOS_ERROR << "Unimplemented case";
    }

    /// Copy the coordinates of the control point to another point
    template<typename TPointType>
    void Copy(TPointType& rPoint) const
    {
        rPoint[0] = X();
        rPoint[1] = Y();
        rPoint[2] = Z();
    }

    // overload operator []
    TDataType& operator[] (int i)
    {
        if (i == 0) { return WX(); }
        else if (i == 1) { return WY(); }
        else if (i == 2) { return WZ(); }
        else if (i == 3) { return BaseType::W(); }
        else
            KRATOS_ERROR << "Out of bound access at i = " << i;
    }

    const TDataType& operator[] (int i) const
    {
        if (i == 0) { return WX(); }
        else if (i == 1) { return WY(); }
        else if (i == 2) { return WZ(); }
        else if (i == 3) { return BaseType::W(); }
        else
            KRATOS_ERROR << "Out of bound access at i = " << i;
    }

    // overload operator ()
    TDataType operator() (int i) const
    {
        if (i == 0) { return X(); }
        else if (i == 1) { return Y(); }
        else if (i == 2) { return Z(); }
        else if (i == 3) { return BaseType::W(); }
        else
            KRATOS_ERROR << "Out of bound access at i = " << i;
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
        p = alpha * c;
        return p;
    }

    /// Multiplication operator
    friend ControlPoint operator*(BaseType c, const TDataType& alpha)
    {
        ControlPoint p;
        p = alpha * c;
        return p;
    }

    /// Multiplication operator
    friend ControlPoint operator*(const Transformation<TDataType>& trans, ControlPoint c)
    {
        c.ApplyTransformation(trans);
        return c;
    }

    /// Compute distance to other control point
    TDataType Distance(const ControlPoint& rOther) const
    {
        return sqrt(pow(X() - rOther.X(), 2) + pow(Y() - rOther.Y(), 2) + pow(Z() - rOther.Z(), 2));
    }

    /// Information
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Control Point";
    }

    void PrintData(std::ostream& rOStream) const override
    {
        // print the control point in homogeneous coordinates
        // rOStream << "(X: " << X() << ", Y: " << Y() << ", Z: " << Z() << ", W: " << BaseType::W() << ")";
        // rOStream << "(" << X() << ", " << Y() << ", " << Z() << ", " << BaseType::W() << ")";
        rOStream << "(" << WX() << ", " << WY() << ", " << WZ() << ", " << BaseType::W() << ")";
    }
};

/// output stream function
template<typename TDataType, typename TWeightType>
inline std::ostream& operator <<(std::ostream& rOStream, const ControlPoint<TDataType, TWeightType>& rThis)
{
    // rThis.PrintInfo(rOStream);
    // rOStream << ": ";
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_CONTROL_POINT_H_INCLUDED
