//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 8 May 2018 $
//   Revision:            $Revision: 1.0 $
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
#include "custom_utilities/import_export/multipatch_exporter.h"

namespace Kratos
{

template<int TDim>
class MultiHBSplinesPatchMatlabExporterWriter : public MultiPatchExporter<TDim>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MultiHBSplinesPatchMatlabExporterWriter);

    virtual void Export(typename Patch<TDim>::Pointer pPatch, std::ostream& rOStream)
    {
        if (pPatch->pFESpace()->Type() != HBSplinesFESpace<TDim>::StaticType())
            KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "only support the hierarchical B-Splines patch")

        // Type definitions
        typedef typename HBSplinesFESpace<TDim>::bf_t bf_t;
        typedef typename HBSplinesFESpace<TDim>::bf_container_t bf_container_t;
        typedef typename HBSplinesFESpace<TDim>::cell_container_t cell_container_t;
        typedef ControlPoint<double> ControlPointType;

        std::size_t patch_id = pPatch->Id();
        rOStream << "%%Information on hierarchical B-Splines patch " << patch_id << "\n\n";

        /* export the basis function information */
        std::size_t cnt = 0;
        std::vector<std::vector<double> > local_knots(TDim);
        double min_xi = static_cast<double>(INT_MAX);
        double max_xi = -min_xi;
        double min_eta = min_xi;
        double max_eta = -min_eta;
        double min_zeta = min_xi;
        double max_zeta = -min_zeta;

        // extract the hierarchical B-Splines space
        typename HBSplinesFESpace<TDim>::Pointer pFESpace = boost::dynamic_pointer_cast<HBSplinesFESpace<TDim> >(pPatch->pFESpace());
        if (pFESpace == NULL)
            KRATOS_THROW_ERROR(std::runtime_error, "The cast to HBSplinesFESpace is failed.", "")

        rOStream << "P" << patch_id << "_params.p1 = " << pFESpace->Order(0) << ";\n";
        rOStream << "P" << patch_id << "_params.p2 = " << pFESpace->Order(1) << ";\n";
        if (TDim == 3)
            rOStream << "P" << patch_id << "_params.p3 = " << pFESpace->Order(2) << ";\n";
        rOStream << "\n";

        for (typename bf_container_t::iterator it_bf = pFESpace->bf_begin(); it_bf != pFESpace->bf_end(); ++it_bf)
        {
            (*it_bf)->LocalKnots(0, local_knots[0]);
            (*it_bf)->LocalKnots(1, local_knots[1]);
            if (TDim == 3) (*it_bf)->LocalKnots(2, local_knots[2]);

            double min_xi_bf = *std::min_element(local_knots[0].begin(), local_knots[0].end());
            double max_xi_bf = *std::max_element(local_knots[0].begin(), local_knots[0].end());
            if (min_xi_bf < min_xi) min_xi = min_xi_bf;
            if (max_xi_bf > max_xi) max_xi = max_xi_bf;

            double min_eta_bf = *std::min_element(local_knots[1].begin(), local_knots[1].end());
            double max_eta_bf = *std::max_element(local_knots[1].begin(), local_knots[1].end());
            if (min_eta_bf < min_eta) min_eta = min_eta_bf;
            if (max_eta_bf > max_eta) max_eta = max_eta_bf;

            if(TDim == 3)
            {
                double min_zeta_bf = *std::min_element(local_knots[2].begin(), local_knots[2].end());
                double max_zeta_bf = *std::max_element(local_knots[2].begin(), local_knots[2].end());
                if (min_zeta_bf < min_zeta) min_zeta = min_zeta_bf;
                if (max_zeta_bf > max_zeta) max_zeta = max_zeta_bf;
            }

            ++cnt;

            rOStream << "P" << patch_id << "_Xi{" << cnt << "} = [";
            for(std::size_t i = 0; i < local_knots[0].size(); ++i)
                rOStream << " " << local_knots[0][i];
            rOStream << "];\n";

            rOStream << "P" << patch_id << "_Eta{" << cnt << "} = [";
            for(std::size_t i = 0; i < local_knots[1].size(); ++i)
                rOStream << " " << local_knots[1][i];
            rOStream << "];\n";

            if(TDim == 3)
            {
                rOStream << "P" << patch_id << "_Zeta{" << cnt << "} = [";
                for(std::size_t i = 0; i < local_knots[2].size(); ++i)
                    rOStream << " " << local_knots[2][i];
                rOStream << "];\n";
            }

            ControlPointType C = (*it_bf)->GetValue(CONTROL_POINT);

            rOStream << "P" << patch_id << "_P(" << cnt << ",:) = [" << C.X() << " " << C.Y() << " " << C.Z() << "];\n";
            rOStream << "P" << patch_id << "_W(" << cnt << ") = " << C.W() << ";\n";
            rOStream << "P" << patch_id << "_Id(" << cnt << ") = " << (*it_bf)->Id() << ";\n";
            rOStream << "P" << patch_id << "_EqId(" << cnt << ") = " << (*it_bf)->EquationId() << ";\n";
            rOStream << std::endl;
        }

        /* export the cell information */
        pFESpace->UpdateCells();

        cnt = 0;
        for(typename cell_container_t::iterator it_cell = pFESpace->pCellManager()->begin(); it_cell != pFESpace->pCellManager()->end(); ++it_cell)
        {
            ++cnt;

            // write the boundary of the cell
            rOStream << "P" << patch_id << "_S{" << cnt << "} = [" << (*it_cell)->LeftValue() << " " << (*it_cell)->RightValue() << ";";
            rOStream << (*it_cell)->DownValue() << " " << (*it_cell)->UpValue() << "];\n";

            // write the extraction operator
            Matrix C = (*it_cell)->GetExtractionOperator();

            rOStream << "P" << patch_id << "_C{" << cnt << "} = [";
            for(std::size_t i = 0; i < C.size1(); ++i)
            {
                for(std::size_t j = 0;  j < C.size2(); ++ j)
                    rOStream << " " << C(i, j);
                rOStream << ";";
            }
            rOStream << "];\n";

            // write the supported basis functions
            const std::vector<std::size_t>& bfs = (*it_cell)->GetSupportedAnchors();
            rOStream << "P" << patch_id << "_N{" << cnt << "} = [";
            for(std::size_t i = 0; i < bfs.size(); ++i)
                rOStream << " " << bfs[i];
            rOStream << "];\n" << std::endl;
        }

        rOStream << std::endl;
    }

    virtual void Export(typename MultiPatch<TDim>::Pointer pMultiPatch, std::ostream& rOStream)
    {
        typedef typename MultiPatch<TDim>::PatchContainerType PatchContainerType;

        for (typename PatchContainerType::ptr_iterator it = pMultiPatch->Patches().ptr_begin(); it != pMultiPatch->Patches().ptr_end(); ++it)
        {
            this->Export(*it, rOStream);
        }
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiHBSplinesPatchMatlabExporterWriter";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};

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
/*        outfile << "close all\n";*/
/*        outfile << "hold on\n";*/
/*        outfile << "axis equal\n\n";*/

        MultiHBSplinesPatchMatlabExporterWriter<TDim> dummy;
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
/*        outfile << "close all\n";*/
/*        outfile << "hold on\n";*/
/*        outfile << "axis equal\n\n";*/

        MultiHBSplinesPatchMatlabExporterWriter<TDim> dummy;
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

