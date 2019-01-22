//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 8 May 2018 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_PBSPLINES_PATCH_MATLAB_EXPORTER_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_PBSPLINES_PATCH_MATLAB_EXPORTER_H_INCLUDED

// System includes
#include <ctime>
#include <vector>
#include <fstream>
#include <iomanip>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/control_point.h"
#include "custom_utilities/import_export/multipatch_exporter.h"

namespace Kratos
{

template<class TFESpaceType>
class MultiPBSplinesPatchMatlabExporter : public MultiPatchExporter<TFESpaceType::Dim()>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MultiPBSplinesPatchMatlabExporter);

    virtual void Export(typename Patch<TFESpaceType::Dim()>::Pointer pPatch, std::ostream& rOStream)
    {
        // extract the point-based B-Splines space
        typename TFESpaceType::Pointer pFESpace = boost::dynamic_pointer_cast<TFESpaceType>(pPatch->pFESpace());
        if (pFESpace == NULL)
        {
            std::stringstream ss;
            ss << "The cast to " << TFESpaceType::StaticType() << " is failed.";
            KRATOS_THROW_ERROR(std::runtime_error, ss.str(), "")
        }

        // Type definitions
        typedef typename TFESpaceType::bf_t bf_t;
        typedef typename TFESpaceType::bf_container_t bf_container_t;
        typedef typename TFESpaceType::cell_container_t cell_container_t;
        typedef ControlPoint<double> ControlPointType;

        // set the export accuracy
        rOStream << std::setprecision(15);

        std::size_t patch_id = pPatch->Id();
        rOStream << "%%Information on " << pPatch->pFESpace()->Type() << " patch " << patch_id << "\n\n";

        /* export the basis function information */
        std::size_t cnt = 0;
        std::vector<std::vector<double> > local_knots(TFESpaceType::Dim());
        double min_xi = static_cast<double>(INT_MAX);
        double max_xi = -min_xi;
        double min_eta = min_xi;
        double max_eta = -min_eta;
        double min_zeta = min_xi;
        double max_zeta = -min_zeta;

        rOStream << "% Degree" << std::endl;
        rOStream << "P" << patch_id << "_params.p1 = " << pFESpace->Order(0) << ";\n";
        if (TFESpaceType::Dim() > 1)
            rOStream << "P" << patch_id << "_params.p2 = " << pFESpace->Order(1) << ";\n";
        if (TFESpaceType::Dim() > 2)
            rOStream << "P" << patch_id << "_params.p3 = " << pFESpace->Order(2) << ";\n";
        rOStream << "\n";

        for (typename bf_container_t::iterator it_bf = pFESpace->bf_begin(); it_bf != pFESpace->bf_end(); ++it_bf)
        {
            (*it_bf)->LocalKnots(0, local_knots[0]);
            if (TFESpaceType::Dim() > 1) (*it_bf)->LocalKnots(1, local_knots[1]);
            if (TFESpaceType::Dim() > 2) (*it_bf)->LocalKnots(2, local_knots[2]);

            double min_xi_bf = *std::min_element(local_knots[0].begin(), local_knots[0].end());
            double max_xi_bf = *std::max_element(local_knots[0].begin(), local_knots[0].end());
            if (min_xi_bf < min_xi) min_xi = min_xi_bf;
            if (max_xi_bf > max_xi) max_xi = max_xi_bf;

            if (TFESpaceType::Dim() > 1)
            {
                double min_eta_bf = *std::min_element(local_knots[1].begin(), local_knots[1].end());
                double max_eta_bf = *std::max_element(local_knots[1].begin(), local_knots[1].end());
                if (min_eta_bf < min_eta) min_eta = min_eta_bf;
                if (max_eta_bf > max_eta) max_eta = max_eta_bf;
            }

            if(TFESpaceType::Dim() > 2)
            {
                double min_zeta_bf = *std::min_element(local_knots[2].begin(), local_knots[2].end());
                double max_zeta_bf = *std::max_element(local_knots[2].begin(), local_knots[2].end());
                if (min_zeta_bf < min_zeta) min_zeta = min_zeta_bf;
                if (max_zeta_bf > max_zeta) max_zeta = max_zeta_bf;
            }

            ++cnt;

            rOStream << "% basis function " << (*it_bf)->Id() << std::endl;
            rOStream << "P" << patch_id << "_Xi{" << cnt << "} = [";
            for(std::size_t i = 0; i < local_knots[0].size(); ++i)
                rOStream << " " << local_knots[0][i];
            rOStream << "];\n";

            if (TFESpaceType::Dim() > 1)
            {
                rOStream << "P" << patch_id << "_Eta{" << cnt << "} = [";
                for(std::size_t i = 0; i < local_knots[1].size(); ++i)
                    rOStream << " " << local_knots[1][i];
                rOStream << "];\n";
            }

            if(TFESpaceType::Dim() > 2)
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
            rOStream << "% cell " << cnt << " information" << std::endl;
            rOStream << "P" << patch_id << "_CId{" << cnt << "} = " << (*it_cell)->Id() << ";\n";
            rOStream << "P" << patch_id << "_S{" << cnt << "} = [" << (*it_cell)->XiMinValue() << " " << (*it_cell)->XiMaxValue();
            if (TFESpaceType::Dim() > 1) rOStream << "; " << (*it_cell)->EtaMinValue() << " " << (*it_cell)->EtaMaxValue();
            if (TFESpaceType::Dim() > 2) rOStream << "; " << (*it_cell)->ZetaMinValue() << " " << (*it_cell)->ZetaMaxValue();
            rOStream << "];\n";

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

        /* visualization */
        rOStream << "% visualize the geometry\n";
        rOStream << "P" << patch_id << "_params.min_xi = " << min_xi << ";\n";
        rOStream << "P" << patch_id << "_params.max_xi = " << 0.999999*max_xi << ";\n";
        rOStream << "P" << patch_id << "_params.num_points1 = 100;\n";
        if (TFESpaceType::Dim() > 1)
        {
            rOStream << "P" << patch_id << "_params.min_eta = " << min_eta << ";\n";
            rOStream << "P" << patch_id << "_params.max_eta = " << 0.999999*max_eta << ";\n";
            rOStream << "P" << patch_id << "_params.num_points2 = 100;\n";
        }
        if (TFESpaceType::Dim() > 2)
        {
            rOStream << "P" << patch_id << "_params.min_zeta = " << min_zeta << ";\n";
            rOStream << "P" << patch_id << "_params.max_zeta = " << 0.999999*max_zeta << ";\n";
            rOStream << "P" << patch_id << "_params.num_points3 = 100;\n";
        }
        rOStream << "plot_geom_hbsplines_" << TFESpaceType::Dim() << "d_cdb(";
        rOStream << "P" << patch_id << "_Xi";
        if (TFESpaceType::Dim() > 1) rOStream << ",P" << patch_id << "_Eta";
        if (TFESpaceType::Dim() > 2) rOStream << ",P" << patch_id << "_Zeta";
        rOStream << ",P" << patch_id << "_P";
        rOStream << ",P" << patch_id << "_W";
        rOStream << ",P" << patch_id << "_params);\n";
        rOStream << "\n";

        rOStream << "% visualize the control points\n";
        rOStream << "figure\n";
        rOStream << "plot_ctrl_points_hbsplines_" << TFESpaceType::Dim() << "d(P" << patch_id << "_P,P" << patch_id << "_W,P" << patch_id << "_EqId+1,P" << patch_id << "_params);\n";
        rOStream << "\n";

        rOStream << "% visualize the cells\n";
        rOStream << "P" << patch_id << "_cell_params.tol = 1.0e-6;\n";
        rOStream << "P" << patch_id << "_cell_params.p1 = " << pFESpace->Order(0) << ";\n";
        rOStream << "P" << patch_id << "_cell_params.max_xi = 1.0;\n";
        if (TFESpaceType::Dim() > 1)
        {
            rOStream << "P" << patch_id << "_cell_params.p2 = " << pFESpace->Order(1) << ";\n";
            rOStream << "P" << patch_id << "_cell_params.max_eta = 1.0;\n";
        }
        if (TFESpaceType::Dim() > 2)
        {
            rOStream << "P" << patch_id << "_cell_params.p3 = " << pFESpace->Order(2) << ";\n";
            rOStream << "P" << patch_id << "_cell_params.max_zeta = 1.0;\n";
        }
        rOStream << "P" << patch_id << "_cell_params.method = 'bezier';\n";
        rOStream << "P" << patch_id << "_cell_params.adjust = 1;\n";
        rOStream << "figure\n";
        rOStream << "plot_hbsplines_cells_" << TFESpaceType::Dim() << "d_with_id(P" << patch_id << "_Xi"
                 << ",P" << patch_id << "_Eta"
                 << ",P" << patch_id << "_P"
                 << ",P" << patch_id << "_W"
                 << ",P" << patch_id << "_EqId"
                 << ",P" << patch_id << "_S"
                 << ",P" << patch_id << "_C"
                 << ",P" << patch_id << "_N"
                 << ",P" << patch_id << "_CId"
                 << ",P" << patch_id << "_cell_params);\n";
        rOStream << "\n";

        rOStream << std::endl;
    }

    virtual void Export(typename MultiPatch<TFESpaceType::Dim()>::Pointer pMultiPatch, std::ostream& rOStream)
    {
        typedef typename MultiPatch<TFESpaceType::Dim()>::PatchContainerType PatchContainerType;

        for (typename PatchContainerType::ptr_iterator it = pMultiPatch->Patches().ptr_begin(); it != pMultiPatch->Patches().ptr_end(); ++it)
        {
            this->Export(*it, rOStream);
        }
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiPBSplinesPatchMatlabExporter";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};

/// output stream function
template<class TFESpaceType>
inline std::ostream& operator <<(std::ostream& rOStream, const MultiPBSplinesPatchMatlabExporter<TFESpaceType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}



} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_PBSPLINES_PATCH_MATLAB_EXPORTER_H_INCLUDED defined

