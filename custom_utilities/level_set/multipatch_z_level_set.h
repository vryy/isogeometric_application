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
//  Date:            22 Jan 2021
//

#if !defined(KRATOS_MULTIPATCH_Z_LEVEL_SET_H_INCLUDED )
#define  KRATOS_MULTIPATCH_Z_LEVEL_SET_H_INCLUDED

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

/// MultiPatch Z level set
/** A level set using multi-patch as level boundary
level_set > 0 if the vertical projection on the multipatch is below the point, and vice versa
*/
class MultiPatchZLevelSet : public LevelSet, public IsogeometricEcho
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MultiPatchZLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(MultiPatchZLevelSet);

    typedef LevelSet BaseType;

    typedef BaseType::PointType PointType;

    typedef MultiPatch<2> MultiPatchType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MultiPatchZLevelSet(typename MultiPatchType::Pointer pMultiPatch)
        : BaseType(), mpMultiPatch(pMultiPatch)
        , mnsampling1(5), mnsampling2(5), mProjectionTolerance(1.0e-10), mMaxIterations(30)
    {
        this->SetEchoLevel(1);
    }

    /// Copy constructor.
    MultiPatchZLevelSet(MultiPatchZLevelSet const& rOther)
        : BaseType(rOther), IsogeometricEcho(rOther)
        , mpMultiPatch(rOther.mpMultiPatch)
        , mnsampling1(rOther.mnsampling1)
        , mnsampling2(rOther.mnsampling2)
        , mProjectionTolerance(rOther.mProjectionTolerance)
        , mMaxIterations(rOther.mMaxIterations)
    {}

    /// Destructor.
    virtual ~MultiPatchZLevelSet() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void SetPredictionSampling(std::size_t nsampling1, std::size_t nsampling2)
    {
        mnsampling1 = nsampling1;
        mnsampling2 = nsampling2;
    }

    void SetProjectionTolerance(double Tol)
    {
        mProjectionTolerance = Tol;
    }

    void SetMaxIterations(int max_iters)
    {
        mMaxIterations = max_iters;
    }

    LevelSet::Pointer CloneLevelSet() const final
    {
        return LevelSet::Pointer(new MultiPatchZLevelSet(*this));
    }

    std::size_t WorkingSpaceDimension() const final
    {
        return 3;
    }

    /// inherit from LevelSet
    double GetValue(const PointType& P) const final
    {
        std::vector<double> local_point(2);
        PointType global_point;
        int target_patch_id;
        int error_code = IsogeometricProjectionUtility::ComputeVerticalProjection(P,
                                    local_point, global_point, target_patch_id,
                                    mpMultiPatch, mProjectionTolerance, mMaxIterations,
                                    mnsampling1, mnsampling2,
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
                std::cout << "WARNING!!!Error computing the vertical projection point on multipatch" << std::endl;
                std::cout << " "; KRATOS_WATCH(P)
                std::cout << " "; KRATOS_WATCH(error_code)
                std::cout << " local_point: " << local_point[0] << ", " << local_point[1] << std::endl;
                std::cout << " "; KRATOS_WATCH(global_point)
            }
        }

        return (P[2] - global_point[2]);
    }

    /// inherit from BRep
    /// Check if a set of points is cut by the level set
    int CutStatus(const std::vector<PointType>& r_points) const final
    {
        typedef typename MultiPatchType::patch_ptr_iterator patch_ptr_iterator;

        std::vector<std::size_t> in_list, out_list, on_list;
        std::vector<double> bounding_box;
        bool is_above, is_below;
        for (std::size_t v = 0; v < r_points.size(); ++v)
        {
            // first simply check with the bounding box of each patch
            is_above = true; is_below = true;
            for (patch_ptr_iterator it = mpMultiPatch->Patches().ptr_begin(); it != mpMultiPatch->Patches().ptr_end(); ++it)
            {
                (*it)->GetBoundingBox(bounding_box);

                if (r_points[v][2] < bounding_box[5])
                {
                    is_above = false;
                }

                if (r_points[v][2] > bounding_box[4])
                {
                    is_below = false;
                }

                if (!is_above && !is_below)
                {
                    break;
                }
            }

            if (is_above)
            {
                out_list.push_back(v);
            }

            if (is_below)
            {
                in_list.push_back(v);
            }

            // if the above or below state cannot be clearly determined, then a projection is necessary
            if (!is_above && !is_below)
            {
                double phi = this->GetValue(r_points[v]);
                if (phi < -this->GetTolerance())
                {
                    in_list.push_back(v);
                }
                else if (phi > this->GetTolerance())
                {
                    out_list.push_back(v);
                }
                else
                {
                    on_list.push_back(v);
                }
            }
        }

        int stat;
        if (in_list.size() == 0 && out_list.size() == 0)
        {
            for (std::size_t v = 0; v < r_points.size(); ++v)
                KRATOS_WATCH(r_points[v])
                KRATOS_WATCH(in_list.size())
                KRATOS_WATCH(out_list.size())
                KRATOS_WATCH(on_list.size())
                KRATOS_WATCH(this->GetTolerance())
                KRATOS_THROW_ERROR(std::logic_error, "!!!FATAL ERROR!!!The geometry is degenerated. We won't handle it.", "")
            }
        else
        {
            if (in_list.size() == 0)
            {
                stat = BRep::_OUT;
                return stat;
            }

            if (out_list.size() == 0)
            {
                stat = BRep::_IN;
                return stat;
            }

            stat = BRep::_CUT;
            return stat;
        }

        return -99; // can't come here. Just to silence the compiler.
    }

    /// inherit from BRep
    /// projects a point on the surface of level_set
    int ProjectOnSurface(const PointType& P, PointType& Proj) const final
    {
        std::vector<double> local_point(2);
        int target_patch_id;
        int error_code = IsogeometricProjectionUtility::ComputeVerticalProjection(P,
                                        local_point, Proj, target_patch_id,
                                        mpMultiPatch, mProjectionTolerance, mMaxIterations,
                                        mnsampling1, mnsampling2,
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
        return "MultiPatch Z Level Set";
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
    std::size_t mnsampling1, mnsampling2;
    double mProjectionTolerance;
    int mMaxIterations;

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
    MultiPatchZLevelSet& operator=(MultiPatchZLevelSet const& rOther);

    ///@}

}; // Class MultiPatchZLevelSet

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MULTIPATCH_Z_LEVEL_SET_H_INCLUDED  defined
