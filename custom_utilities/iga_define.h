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
#include <vector>

namespace Kratos
{

#if defined(__GNUC__) || defined(__clang__)
#define IGA_DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define IGA_DEPRECATED __declspec(deprecated)
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define IGA_DEPRECATED
#endif

enum BoundarySide
{
    _BLEFT_   = 0,
    _BRIGHT_  = 1,
    _BTOP_    = 2,
    _BBOTTOM_ = 3,
    _BFRONT_  = 4,
    _BBACK_   = 5,
    _NUMBER_OF_BOUNDARY_SIDE = 6
};

enum BoundarySide2D
{
    _B2LEFT_   = 0,
    _B2RIGHT_  = 1,
    _B2TOP_    = 2,
    _B2BOTTOM_ = 3,
};

enum BoundarySide3D
{
    _B3LEFT_   = 0,
    _B3RIGHT_  = 1,
    _B3TOP_    = 2,
    _B3BOTTOM_ = 3,
    _B3FRONT_  = 4,
    _B3BACK_   = 5,
};

#define BOUNDARY_FLAG(x) (1 << (x+1))

enum BoundaryDirection
{
    _FORWARD_ = 0,
    _REVERSED_ = 1
};

enum BoundaryFlag
{
    _FLEFT_   = BOUNDARY_FLAG(_BLEFT_),
    _FRIGHT_  = BOUNDARY_FLAG(_BRIGHT_),
    _FTOP_    = BOUNDARY_FLAG(_BTOP_),
    _FBOTTOM_ = BOUNDARY_FLAG(_BBOTTOM_),
    _FFRONT_  = BOUNDARY_FLAG(_BFRONT_),
    _FBACK_   = BOUNDARY_FLAG(_BBACK_)
};

inline std::string BoundarySideName(const BoundarySide& side)
{
    switch(side)
    {
        case _BLEFT_:    return "left";
        case _BRIGHT_:   return "right";
        case _BTOP_:     return "top";
        case _BBOTTOM_:  return "bottom";
        case _BFRONT_:   return "front";
        case _BBACK_:    return "back";
        default:        return "inner";
    }
}

template<int TDim>
struct ParameterDirection
{
    static std::vector<int> Get(const BoundarySide& side)
    {
        std::stringstream ss;
        ss << __FUNCTION__ << " is not implemented for " << TDim << "D";
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }
};

template<>
struct ParameterDirection<2>
{
    static int Get_(const BoundarySide& side)
    {
        switch(side)
        {
            case _BLEFT_:    return 1;
            case _BRIGHT_:   return 1;
            case _BTOP_:     return 0;
            case _BBOTTOM_:  return 0;
            default:        return -1;
        }
        return -1;
    }

    static std::vector<int> Get(const BoundarySide& side)
    {
        return std::vector<int>{Get_(side)};
    }
};

template<>
struct ParameterDirection<3>
{
    static std::vector<int> Get(const BoundarySide& side)
    {
        switch(side)
        {
            case _BLEFT_:    return std::vector<int>{1, 2};
            case _BRIGHT_:   return std::vector<int>{1, 2};
            case _BFRONT_:   return std::vector<int>{0, 2};
            case _BBACK_:    return std::vector<int>{0, 2};
            case _BTOP_:     return std::vector<int>{0, 1};
            case _BBOTTOM_:  return std::vector<int>{0, 1};
            default:        return std::vector<int>{-1, -1};
        }
        return std::vector<int>{-1, -1};
    }
};

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

