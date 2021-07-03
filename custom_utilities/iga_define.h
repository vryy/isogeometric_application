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
#include "includes/define.h"
#ifdef SD_APP_FORWARD_COMPATIBILITY
#include <memory>
#else
#include <boost/pointer_cast.hpp>
#endif

namespace iga
{
    // alias for dynamic_pointer_cast
    #ifdef SD_APP_FORWARD_COMPATIBILITY
    template <class T, class U>
    inline std::shared_ptr<T> dynamic_pointer_cast(const std::shared_ptr<U> &r) noexcept
    {
        return std::dynamic_pointer_cast<T>(r);
    }
    #else
    template <class T, class U>
    inline boost::shared_ptr<T> dynamic_pointer_cast(const boost::shared_ptr<U> &r) noexcept
    {
        return boost::dynamic_pointer_cast<T>(r);
    }
    #endif

    // alias for make_shared
    #ifdef SD_APP_FORWARD_COMPATIBILITY
    template<typename C, typename...Args>
    inline std::shared_ptr<C> make_shared(Args &&...args)
    {
        return std::make_shared<C>(std::forward<Args>(args)...);
    }
    #else
    template<typename C, typename...Args>
    inline boost::shared_ptr<C> make_shared(Args &&...args)
    {
        return boost::make_shared<C>(std::forward<Args>(args)...);
    }
    #endif

    // wrapper for object and pointer to that object
    #ifdef SD_APP_FORWARD_COMPATIBILITY
    template<typename T, typename PT = Kratos::shared_ptr<T> >
    #else
    template<typename T, typename PT = boost::shared_ptr<T> >
    #endif
    struct Wrapper
    {
        Wrapper(PT p) : mp(p) {}
        PT& GetPointer() {return mp;}
        T& GetReference() {return *mp;}
        PT mp;
    };
} // namespace iga

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

struct IsogeometricEcho
{
    KRATOS_CLASS_POINTER_DEFINITION(IsogeometricEcho);

    IsogeometricEcho() : mEchoLevel(0) {}
    ~IsogeometricEcho() {}

    IsogeometricEcho(IsogeometricEcho const& rOther)
    : mEchoLevel(rOther.mEchoLevel)
    {}

    void SetEchoLevel(const int& level) {mEchoLevel = level;}

    const int& GetEchoLevel() const {return mEchoLevel;}

    static bool Has(const int& echo_level, const int& echo_flag)
    {
        return ((echo_level & echo_flag) == echo_flag);
    }

    int mEchoLevel;
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
 * Helper struct to extract the pointer type
 * One case use typename Isogeometric_Pointer_Helper<TType>::Pointer as replacement for typename TType::Pointer
 * This is useful when TType is used in nested template argument, e.g. template<..., typename type_t<TType> >
 */
template<class TType>
struct Isogeometric_Pointer_Helper
{
    KRATOS_CLASS_POINTER_DEFINITION(TType);
};

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_IGA_DEFINE_H_INCLUDED defined
