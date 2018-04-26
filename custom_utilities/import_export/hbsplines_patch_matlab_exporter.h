//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_NURBS_PATCH_MATLAB_EXPORTER_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_NURBS_PATCH_MATLAB_EXPORTER_H_INCLUDED

// System includes
#include <vector>
#include <fstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"
#include "custom_utilities/import_export/multipatch_exporter.h"

namespace Kratos
{

class HBSplinesPatchMatlabExporter_Helper
{
public:
    template<int TDim>
    static void WriteMatlabControlPoints(std::ostream& rOStream, typename Patch<TDim>::Pointer pPatch, const std::string& var_name)
    {
        KRATOS_THROW_ERROR(std::logic_error, "WriteMatlabControlPoints is not implemented for dimension", TDim)
    }
};

class HBSplinesPatchMatlabExporter
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(HBSplinesPatchMatlabExporter);

    template<int TDim>
    static void Export(typename Patch<TDim>::Pointer pPatch, const std::string& filename)
    {
        HBSplinesPatchMatlabExporterWriter<TDim> dummy;
        dummy.Export(pPatch, filename);
    }

    template<int TDim>
    static void Export(typename MultiPatch<TDim>::Pointer pMultiPatch, const std::string& filename)
    {
        HBSplinesPatchMatlabExporterWriter<TDim> dummy;
        dummy.Export(pMultiPatch, filename);
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "HBSplinesPatchMatlabExporter";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};


/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const HBSplinesPatchMatlabExporter& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}



} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_NURBS_PATCH_MATLAB_EXPORTER_H_INCLUDED defined

