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
#include "custom_utilities/level_set/isogeometric_projection_utility.h"
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
    , mnsampling1(5), mnsampling2(5), mTolerance(1.0e-10), mMaxIterations(30)
    {
        this->SetEchoLevel(1);
    }

    /// Copy constructor.
    MultiPatchZLevelSet(MultiPatchZLevelSet const& rOther)
    : BaseType(rOther), IsogeometricEcho(rOther)
    , mpMultiPatch(rOther.mpMultiPatch)
    , mnsampling1(rOther.mnsampling1)
    , mnsampling2(rOther.mnsampling2)
    {}

    /// Destructor.
    virtual ~MultiPatchZLevelSet() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    void SetPredictionSampling(const std::size_t& nsampling1, const std::size_t& nsampling2)
    {
        mnsampling1 = nsampling1;
        mnsampling2 = nsampling2;
    }


    void SetTolerance(const double& Tol)
    {
        mTolerance = Tol;
    }


    void SetMaxIterations(const int& max_iters)
    {
        mMaxIterations = max_iters;
    }


    virtual LevelSet::Pointer CloneLevelSet() const
    {
        return LevelSet::Pointer(new MultiPatchZLevelSet(*this));
    }


    virtual std::size_t WorkingSpaceDimension() const
    {
        return 3;
    }


    virtual double GetValue(const PointType& P) const
    {
        std::vector<double> local_point(2);
        PointType global_point;
        int target_patch_id;
        int error_code = IsogeometricProjectionUtility::ComputeVerticalProjection(P,
            local_point, global_point, target_patch_id,
            mpMultiPatch, mTolerance, mMaxIterations,
            mnsampling1, mnsampling2,
            this->GetEchoLevel()-2);
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
                std::cout << "WARNING!!!Error computing the vertical point projection point on multipatch" << std::endl;
                std::cout << " "; KRATOS_WATCH(P)
                std::cout << " "; KRATOS_WATCH(error_code)
                std::cout << " local_point: " << local_point[0] << ", " << local_point[1] << std::endl;
                std::cout << " "; KRATOS_WATCH(global_point)
            }
        }
        return (P[2] - global_point[2]);
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
    virtual std::string Info() const
    {
        return "MultiPatch Z Level Set";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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
    double mTolerance;
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


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                MultiPatchZLevelSet& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const MultiPatchZLevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MULTIPATCH_Z_LEVEL_SET_H_INCLUDED  defined
