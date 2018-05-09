//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_IGA_DEFINE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_IGA_DEFINE_H_INCLUDED

namespace Kratos
{

enum BoundarySide
{
    _LEFT_   = 0,
    _RIGHT_  = 1,
    _TOP_    = 2,
    _BOTTOM_ = 3,
    _FRONT_  = 4,
    _BACK_   = 5,
    _NUMBER_OF_BOUNDARY_SIDE = 6
};

#define BOUNDARY_FLAG(x) (1 << (x+1))

enum IsogeometricEchoFlags
{
    ECHO_REFINEMENT   = 0b0000000000000001,
};

enum PreElementType
{
    _NURBS_ = 0,
    _BEZIER_ = 1
};

enum PostElementType
{
    _TRIANGLE_ = 0,
    _QUADRILATERAL_ = 1,
    _TETRAHEDRA_ = 2,
    _HEXAHEDRA_ = 3
};

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_IGA_DEFINE_H_INCLUDED defined

