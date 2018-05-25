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

enum BoundaryDirection
{
    _FORWARD_ = 0,
    _REVERSED_ = 1
};

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
        case _BOTTOM_:  return "bottom";
        case _FRONT_:   return "front";
        case _BACK_:    return "back";
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
            case _LEFT_:    return 1;
            case _RIGHT_:   return 1;
            case _TOP_:     return 0;
            case _BOTTOM_:  return 0;
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
            case _LEFT_:    return std::vector<int>{1, 2};
            case _RIGHT_:   return std::vector<int>{1, 2};
            case _FRONT_:   return std::vector<int>{0, 2};
            case _BACK_:    return std::vector<int>{0, 2};
            case _TOP_:     return std::vector<int>{0, 1};
            case _BOTTOM_:  return std::vector<int>{0, 1};
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

