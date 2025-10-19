//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 8 May 2018 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_HBSPLINES_PATCH_MATLAB_EXPORTER_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_HBSPLINES_PATCH_MATLAB_EXPORTER_H_INCLUDED

// System includes
#include <ctime>
#include <vector>
#include <fstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/hbsplines/hbsplines_fespace.h"
#include "custom_utilities/import_export/multi_pbsplines_patch_matlab_exporter.h"

namespace Kratos
{

class MultiHBSplinesPatchMatlabExporter
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MultiHBSplinesPatchMatlabExporter);

    template<int TDim>
    static void Export(typename Patch<TDim>::Pointer pPatch, const std::string& filename)
    {
        std::ofstream outfile;
        outfile.open(filename, std::ios::out);

        std::time_t t = std::time(0);   // get time now
        std::tm* now = std::localtime(&t);

        outfile << "%% hierarchical B-Splines mesh information, (c) Hoang Giang Bui, " << (now->tm_year + 1900) << "\n";
        outfile << "clc\n";
        outfile << "clear\n";
        outfile << "%close all\n";
        outfile << "%hold on\n";
        outfile << "%axis equal\n\n";

        MultiPBSplinesPatchMatlabExporter<HBSplinesFESpace<TDim> > dummy;
        dummy.Export(pPatch, outfile);

        outfile.close();
        std::cout << "Export patch information to " << filename << " completed" << std::endl;
    }

    template<int TDim>
    static void Export(typename MultiPatch<TDim>::Pointer pMultiPatch, const std::string& filename)
    {
        std::ofstream outfile;
        outfile.open(filename, std::ios::out);

        std::time_t t = std::time(0);   // get time now
        std::tm* now = std::localtime(&t);

        outfile << "%% hierarchical B-Splines mesh information, (c) Hoang Giang Bui, " << (now->tm_year + 1900) << "\n";
        outfile << "clc\n";
        outfile << "clear\n";
        outfile << "%close all\n";
        outfile << "%hold on\n";
        outfile << "%axis equal\n\n";

        MultiPBSplinesPatchMatlabExporter<HBSplinesFESpace<TDim> > dummy;
        dummy.Export(pMultiPatch, outfile);

        outfile.close();
        std::cout << "Export multipatch information to " << filename << " completed" << std::endl;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiHBSplinesPatchMatlabExporter";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const MultiHBSplinesPatchMatlabExporter& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_HBSPLINES_PATCH_MATLAB_EXPORTER_H_INCLUDED defined
