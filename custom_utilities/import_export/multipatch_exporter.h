//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_EXPORTER_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_EXPORTER_H_INCLUDED

// System includes
#include <vector>
#include <fstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/patch.h"

namespace Kratos
{

/**
Abstract class to export the MultiPatch Geometry to various visualization framework
 */
template<int TDim>
class MultiPatchExporter
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiPatchExporter);

    /// Type definition
    typedef KnotArray1D<double> knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /// Default constructor
    MultiPatchExporter() : mAccuracy(15) {}

    /// Destructor
    virtual ~MultiPatchExporter() {}

    /// Set the accuracy
    void SetAccuracy(const std::size_t& Accuracy) {mAccuracy = Accuracy;}

    /// Get the accuracy
    const std::size_t& Accuracy() const {return mAccuracy;}

    /// Export a single patch
    virtual void Export(typename Patch<TDim>::Pointer pPatch, std::ostream& rOStream) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Export a multipatch
    virtual void Export(typename MultiPatch<TDim>::Pointer pMultiPatch, std::ostream& rOStream) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiPatchExporter";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:

    std::size_t mAccuracy; // the number of digits after comma

}; // end class MultiPatchExporter

/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const MultiPatchExporter<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_EXPORTER_H_INCLUDED defined

