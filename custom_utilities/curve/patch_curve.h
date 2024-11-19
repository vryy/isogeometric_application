//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         isogeometric_application/LICENSE.txt
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//  Date:            14 Nov 2024
//

#if !defined(KRATOS_PATCH_CURVE_H_INCLUDED )
#define  KRATOS_PATCH_CURVE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_utilities/patch.h"
#include "brep_application/custom_algebra/curve/curve.h"

namespace Kratos
{
///@addtogroup BRepApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/** Curve representation on the patch
 */
template<int TDim>
class PatchCurve : public Curve
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PatchCurve
    KRATOS_CLASS_POINTER_DEFINITION(PatchCurve);

    typedef FunctionR1R3 SuperType;

    typedef Curve BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PatchCurve(typename Patch<TDim>::ConstPointer pPatch)
        : BaseType()
        , mpPatch(pPatch)
    {
    }

    /// Copy constructor.
    PatchCurve(PatchCurve const& rOther)
        : BaseType(rOther)
        , mpPatch(rOther.mpPatch)
        , mLocalCoords(rOther.mLocalCoords)
        , mLocalDirs(rOther.mLocalDirs)
    {}

    /// Destructor.
    ~PatchCurve() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Set the evaluation information
    void SetLocalCoord(const int ldir, const int dir, const double v)
    {
        if (ldir > TDim - 2)
        {
            KRATOS_ERROR << "Invalid parameter direction " << ldir;
        }

        mLocalDirs[ldir] = dir;
        mLocalCoords[ldir] = v;
    }

    /// inherit from Function
    SuperType::Pointer CloneFunction() const final
    {
        return SuperType::Pointer(new PatchCurve(*this));
    }

    /// inherit from Curve
    Curve::Pointer Clone() const final
    {
        return BaseType::Pointer(new PatchCurve(*this));
    }

    /// inherit from Function
    OutputType GetValue(const InputType& t) const final
    {
        OutputType P;

        std::vector<double> xi = CreateCoordinates(t);

        const auto pGridFunc = mpPatch->pGetGridFunction(CONTROL_POINT_COORDINATES);

        noalias(P) = pGridFunc->GetValue(xi);

        return std::move(P);
    }

    /// inherit from Function
    OutputType GetDerivative(const int& component, const InputType& t) const final
    {
        OutputType P;

        std::vector<double> xi = CreateCoordinates(t);

        const auto pGridFunc = mpPatch->pGetGridFunction(CONTROL_POINT_COORDINATES);
        const auto Cvec = pGridFunc->GetDerivative(xi);

        noalias(P) = Cvec[0];

        return std::move(P);
    }

    /// inherit from Function
    OutputType GetSecondDerivative(const int& component_1, const int& component_2, const InputType& t) const final
    {
        OutputType P;

        std::vector<double> xi = CreateCoordinates(t);

        const auto pGridFunc = mpPatch->pGetGridFunction(CONTROL_POINT_COORDINATES);
        const auto Cvec = pGridFunc->GetSecondDerivative(xi);

        noalias(P) = Cvec[0];

        return std::move(P);
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const final
    {
        return "Patch Curve";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    typename Patch<TDim>::ConstPointer mpPatch;
    std::array < double, TDim - 1 > mLocalCoords;
    std::array < int, TDim - 1 > mLocalDirs;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    std::vector<double> CreateCoordinates(const InputType t) const
    {
        std::vector<double> xi(TDim);

        if constexpr (TDim == 1)
        {
            xi[0] = t;
        }
        else if constexpr (TDim == 2)
        {
            xi[mLocalDirs[0]] = mLocalCoords[0];
            xi[1 - mLocalDirs[0]] = t;
        }
        else if constexpr (TDim == 3)
        {
            xi[mLocalDirs[0]] = mLocalCoords[0];
            xi[mLocalDirs[1]] = mLocalCoords[1];
            xi[3 - mLocalDirs[0] - mLocalDirs[1]] = t;
        }

        return std::move(xi);
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    PatchCurve& operator=(PatchCurve const& rOther);

    ///@}

}; // Class PatchCurve

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_PATCH_CURVE_H_INCLUDED  defined
