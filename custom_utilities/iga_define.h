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
#if defined(ISOGEOMETRIC_USE_STD_ANY)
#include <any>
#include <utility>
#elif defined(ISOGEOMETRIC_USE_BOOST_ANY)
#include <boost/any.hpp>
#include <boost/move/utility.hpp>
#else
#include <boost/any.hpp>
#include <boost/move/utility.hpp>
#endif


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

enum ParametricAxis
{
    _PA_U_   = 0,
    _PA_V_   = 1,
    _PA_W_   = 2,
    _NUMBER_OF_PARAMETRIC_AXIS = 3
};

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
    _REVERSED_ = 1,
    _UNDEFINED_DIR_ = -1
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
            case _BLEFT_:    return _PA_V_;
            case _BRIGHT_:   return _PA_V_;
            case _BTOP_:     return _PA_U_;
            case _BBOTTOM_:  return _PA_U_;
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
            case _BLEFT_:    return std::vector<int>{_PA_V_, _PA_W_}; // v, w
            case _BRIGHT_:   return std::vector<int>{_PA_V_, _PA_W_};
            case _BFRONT_:   return std::vector<int>{_PA_U_, _PA_W_}; // u, w
            case _BBACK_:    return std::vector<int>{_PA_U_, _PA_W_};
            case _BTOP_:     return std::vector<int>{_PA_U_, _PA_V_}; // u, v
            case _BBOTTOM_:  return std::vector<int>{_PA_U_, _PA_V_};
            default:        return std::vector<int>{-1, -1};
        }
        return std::vector<int>{-1, -1};
    }
};

struct ReversedBoundarySide
{
    static BoundarySide Get(const BoundarySide& side)
    {
        switch(side)
        {
            case _BLEFT_:    return _BRIGHT_;
            case _BRIGHT_:   return _BLEFT_;
            case _BFRONT_:   return _BBACK_;
            case _BBACK_:    return _BFRONT_;
            case _BTOP_:     return _BBOTTOM_;
            case _BBOTTOM_:  return _BTOP_;
        }
    }
};

enum IsogeometricEchoFlags
{
    ECHO_REFINEMENT   = 0b0000000000000001,
    ECHO_REFINEMENT_DETAIL   = 0b0000000000000010,
};

struct IsogeometricEchoCheck
{
    static bool Has(const int& echo_level, const int& echo_flag)
    {
        return ((echo_level & echo_flag) == echo_flag);
    }
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

/**
 * Macro to define associated pointers
 */
#define ISOGEOMETRIC_CLASS_POINTER_DEFINITION(a) typedef Kratos::shared_ptr<a > Pointer; \
typedef Kratos::shared_ptr<a > SharedPointer; \
typedef Kratos::shared_ptr<const a > ConstPointer; \
typedef Kratos::weak_ptr<a > WeakPointer; \
typedef Kratos::unique_ptr<a > UniquePointer

/**
 * Helper struct to extract the pointer type
 * One case use typename Isogeometric_Pointer_Helper<TType>::Pointer as replacement for typename TType::Pointer
 * This is useful when TType is used in nested template argument, e.g. template<..., typename type_t<TType> >
 */
template<class TType>
struct Isogeometric_Pointer_Helper
{
    ISOGEOMETRIC_CLASS_POINTER_DEFINITION(TType);
};

#if defined(ISOGEOMETRIC_USE_STD_ANY)
using any = std::any;
using bad_any_cast = std::bad_any_cast;
template<typename T, typename...Args>
T any_cast(Args &&...args) {
    return std::any_cast<T>(std::forward<Args>(args)...);
}
#elif defined(ISOGEOMETRIC_USE_BOOST_ANY)
using any = boost::any;
using bad_any_cast = boost::bad_any_cast;
template<typename T, typename...Args>
T any_cast(Args &&...args) {
    return boost::any_cast<T>(boost::forward<Args>(args)...);
}
#else
using any = boost::any;
using bad_any_cast = boost::bad_any_cast;
template<typename T, typename...Args>
T any_cast(Args &&...args) {
    return boost::any_cast<T>(boost::forward<Args>(args)...);
}
#endif

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_IGA_DEFINE_H_INCLUDED defined

