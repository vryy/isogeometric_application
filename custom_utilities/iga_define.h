//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_IGA_DEFINE_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_IGA_DEFINE_H_INCLUDED

#include <cstring>

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

enum BoundaryFlag
{
    _FLEFT_   = BOUNDARY_FLAG(_LEFT_),
    _FRIGHT_  = BOUNDARY_FLAG(_RIGHT_),
    _FTOP_    = BOUNDARY_FLAG(_TOP_),
    _FBOTTOM_ = BOUNDARY_FLAG(_BOTTOM_),
    _FFRONT_  = BOUNDARY_FLAG(_FRONT_),
    _FBACK_   = BOUNDARY_FLAG(_BACK_)
};

inline std::string BoundarySideName(const BoundarySide& side)
{
    switch(side)
    {
        case _LEFT_:    return "left";
        case _RIGHT_:   return "right";
        case _TOP_:     return "top";
        case _BOTTOM_: return "bottom";
        case _FRONT_:   return "front";
        case _BACK_:    return "back";
        default:        return "inner";
    }
}

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

