//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         brep_application/LICENSE.txt
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//  Date:            11 Oct 2024
//

#if !defined(KRATOS_MULTIPATCH_NORMAL_LEVEL_SET_H_INCLUDED )
#define  KRATOS_MULTIPATCH_NORMAL_LEVEL_SET_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/iga_define.h"
#include "custom_utilities/multipatch.h"
#include "custom_utilities/isogeometric_projection_utility.h"
#include "custom_algebra/level_set/level_set.h"
#include "custom_algebra/curve/curve.h"
#include "custom_algebra/curve/linear_curve.h"

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

/// MultiPatch level set
/** A level set using multi-patch as level boundary
 * It uses normal projection to evaluate the level set value
 */
template<int TDim>
class MultiPatchNormalLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MultiPatchNormalLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(MultiPatchNormalLevelSet);

    typedef LevelSet BaseType;

    typedef BaseType::PointType PointType;

    typedef MultiPatch<TDim> MultiPatchType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MultiPatchNormalLevelSet(typename MultiPatchType::Pointer pMultiPatch)
        : BaseType(), mpMultiPatch(pMultiPatch)
        , mProjectionTolerance(1e-10), mMaxIterations(30)
    {
        this->SetEchoLevel(1);
        for (int i = 0; i < TDim; ++i)
            mnsampling[i] = 5;
        mInnerPoint.clear();
        mpCurve = Curve::Pointer(new LinearCurve(0.0, 0.0, 0.0, 0.0, 0.0, 1.0));
    }

    /// Copy constructor.
    MultiPatchNormalLevelSet(MultiPatchNormalLevelSet const& rOther)
        : BaseType(rOther)
        , mpMultiPatch(rOther.mpMultiPatch)
        , mnsampling(rOther.mnsampling)
        , mProjectionTolerance(rOther.mProjectionTolerance)
        , mMaxIterations(rOther.mMaxIterations)
        , mInnerPoint(rOther.mInnerPoint)
        , mpCurve(rOther.mpCurve)
    {}

    /// Destructor.
    virtual ~MultiPatchNormalLevelSet() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void SetPredictionSampling(int i, int nsampling)
    {
        if (i < TDim && i > -1)
            mnsampling[i] = nsampling;
    }

    void SetProjectionTolerance(double Tol)
    {
        mProjectionTolerance = Tol;
    }

    void SetMaxIterations(int max_iters)
    {
        mMaxIterations = max_iters;
    }

    template<typename TPointType>
    void SetInnerPoint(const TPointType& rPoint)
    {
        for (int i = 0; i < TDim; ++i)
            mInnerPoint[i] = rPoint[i];
    }

    void SetCurve(Curve::Pointer pCurve)
    {
        mpCurve = pCurve;
    }

    LevelSet::Pointer CloneLevelSet() const final
    {
        return LevelSet::Pointer(new MultiPatchNormalLevelSet(*this));
    }

    std::size_t WorkingSpaceDimension() const final
    {
        return 3;
    }

    /// inherit from LevelSet
    double GetValue(const PointType& P) const final
    {
        std::vector<double> local_point(TDim);
        PointType global_point;
        int target_patch_id;
        int error_code = IsogeometricProjectionUtility<PointType, TDim>::ComputeNormalProjection(P,
                                local_point, global_point, target_patch_id,
                                mpMultiPatch, mProjectionTolerance, mMaxIterations,
                                mnsampling,
                                this->GetEchoLevel() - 2);

        if (this->GetEchoLevel() > 1)
        {
            std::cout << "local_point: " << local_point[0] << ", " << local_point[1] << std::endl;
            KRATOS_WATCH(global_point)
            KRATOS_WATCH(target_patch_id)
        }

        if (error_code != 0)
        {
            if (this->GetEchoLevel() > 0)
            {
                std::cout << "WARNING!!!Error computing the normal projection point on multipatch" << std::endl;
                std::cout << " "; KRATOS_WATCH(P)
                std::cout << " "; KRATOS_WATCH(error_code)
                std::cout << " local_point: " << local_point[0] << ", " << local_point[1] << std::endl;
                std::cout << " "; KRATOS_WATCH(global_point)
            }
        }

        // check with inner point
        double test;

        if constexpr (TDim == 1)
        {
            test = inner_prod(P - global_point, mInnerPoint - global_point);
        }
        else if constexpr (TDim == 2)
        {
            PointType projection;
            mpCurve->ProjectOnCurve(P, projection);
            test = inner_prod(P - global_point, projection - global_point);
        }

        if (test < 0)
            return norm_2(P - global_point);
        else
            return -norm_2(P - global_point);
    }

    /// inherit from BRep
    /// projects a point on the surface of level_set
    int ProjectOnSurface(const PointType& P, PointType& Proj) const final
    {
        std::vector<double> local_point(TDim);
        int target_patch_id;
        int error_code = IsogeometricProjectionUtility<PointType, TDim>::ComputeNormalProjection(P,
                                local_point, Proj, target_patch_id,
                                mpMultiPatch, mProjectionTolerance, mMaxIterations,
                                mnsampling,
                                this->GetEchoLevel() - 2);
        return error_code;
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
        return "MultiPatch Level Set";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        rOStream << "MultiPatch: " << *mpMultiPatch;
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

    typename MultiPatchType::Pointer mpMultiPatch;
    std::array<unsigned int, TDim> mnsampling;
    double mProjectionTolerance;
    int mMaxIterations;
    PointType mInnerPoint;
    Curve::Pointer mpCurve;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

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
    MultiPatchNormalLevelSet& operator=(MultiPatchNormalLevelSet const& rOther);

    ///@}

}; // Class MultiPatchNormalLevelSet

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MULTIPATCH_NORMAL_LEVEL_SET_H_INCLUDED  defined
