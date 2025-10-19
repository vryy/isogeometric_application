//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 9 Nov 2018 $
//
//

// Project includes
#include "custom_utilities/import_export/tsplines_patch_tsm_importer.h"
#include "custom_utilities/tsplines/tscell.h"
#include "custom_utilities/nurbs/pbbsplines_basis_function.h"
#include "custom_utilities/tsplines/tsplines_fespace.h"

// External includes
#include "rhbuilder.h"
#include "extractor.h"

namespace Kratos
{

Patch<2>::Pointer TSplinesPatchTSMImporter::ImportSingle(const std::string& filename) const
{
    std::cout << "Start reading TSM in " << filename << std::endl;
    RhBuilderPtr reader = TSPLINE::makePtr<RhBuilder>(filename);
    std::cout << "Start to find T-Splines in " << filename << std::endl;
    TSplinePtr tspline = reader->findTSpline();
    std::cout << "Reading T-Splines in " << filename << " completed" << std::endl;

    // create the PBSplinesFESpaceType
    typedef PBBSplinesBasisFunction<2, TsCell> PBBSplinesBasisFunctionType;
    typedef TSplinesFESpace<2, double, PBBSplinesBasisFunctionType, BCellManager<2, TsCell> > TSplinesFESpaceType;
    typedef typename Patch<2>::ControlPointType ControlPointType;
    typedef typename TSplinesFESpaceType::knot_t knot_t;
    typedef typename TSplinesFESpaceType::cell_t cell_t;

    typename TSplinesFESpaceType::Pointer pNewFESpace = TSplinesFESpaceType::Create();

    // set order information
    pNewFESpace->SetInfo(0, tspline->getSDegree());
    pNewFESpace->SetInfo(1, tspline->getTDegree());

    // extra step: determine xi_max and eta_max
    TPointsetPtr tpoints = tspline->getTPointset();
    // WATCH(tpoints->size())

    Real min_xi_min = 1.0e9, min_eta_min = 1.0e9;
    Real max_xi_max = -1.0e9, max_eta_max = -1.0e9;
    for (TObjVIterator pit = tpoints->iteratorBegin(); pit != tpoints->iteratorEnd(); ++pit)
    {
        TPointPtr ptpoint = std::dynamic_pointer_cast<TSPLINE::TPoint>(*pit);
        TNodeV4Ptr pnodev4 = std::dynamic_pointer_cast<TNodeV4>(ptpoint->getTNode());
        std::vector<std::vector<Real> > knot_vecs(2);
        TExtractor::extractUVKnotsFromTNodeV4(pnodev4, knot_vecs[0], knot_vecs[1]);

        for (std::size_t i1 = 0; i1 < knot_vecs[0].size(); ++i1)
        {
            if (knot_vecs[0][i1] > max_xi_max)
            {
                max_xi_max = knot_vecs[0][i1];
            }
            if (knot_vecs[0][i1] < min_xi_min)
            {
                min_xi_min = knot_vecs[0][i1];
            }
        }

        for (std::size_t j1 = 0; j1 < knot_vecs[1].size(); ++j1)
        {
            if (knot_vecs[1][j1] > max_eta_max)
            {
                max_eta_max = knot_vecs[1][j1];
            }
            if (knot_vecs[1][j1] < min_eta_min)
            {
                min_eta_min = knot_vecs[1][j1];
            }
        }
    }
    KRATOS_WATCH(min_xi_min)
    KRATOS_WATCH(max_xi_max)
    KRATOS_WATCH(min_eta_min)
    KRATOS_WATCH(max_eta_max)

    /** create basis functions **/
    std::size_t id = 0;
    double x, y, z, w;
    const double knot_tol = 1.0e-11; // tolerance to round off the knot to 0 or 1. We should parameterize it. TODO
    const double area_tol = 1.0e-6; // tolerance to accept the nonzero-area cell. We should parameterize it. TODO
    for (TObjVIterator pit = tpoints->iteratorBegin(); pit != tpoints->iteratorEnd(); ++pit)
    {
        ++id;
        TPointPtr ptpoint = std::dynamic_pointer_cast<TSPLINE::TPoint>(*pit);
        ptpoint->setId(id); // set id one time

        TNodeV4Ptr pnodev4 = std::dynamic_pointer_cast<TNodeV4>(ptpoint->getTNode());
        std::vector<std::vector<Real> > knot_vecs(2);
        TExtractor::extractUVKnotsFromTNodeV4(pnodev4, knot_vecs[0], knot_vecs[1]);

        // scale the knot vectors
        for (std::size_t i1 = 0; i1 < knot_vecs[0].size(); ++i1)
        {
            knot_vecs[0][i1] = (knot_vecs[0][i1] - min_xi_min) / (max_xi_max - min_xi_min);
            if (fabs(knot_vecs[0][i1]) < knot_tol)
            {
                knot_vecs[0][i1] = 0.0;
            }
            if (fabs(knot_vecs[0][i1] - 1.0) < knot_tol)
            {
                knot_vecs[0][i1] = 1.0;
            }
        }

        for (std::size_t j1 = 0; j1 < knot_vecs[1].size(); ++j1)
        {
            knot_vecs[1][j1] = (knot_vecs[1][j1] - min_eta_min) / (max_eta_max - min_eta_min);
            if (fabs(knot_vecs[1][j1]) < knot_tol)
            {
                knot_vecs[1][j1] = 0.0;
            }
            if (fabs(knot_vecs[1][j1] - 1.0) < knot_tol)
            {
                knot_vecs[1][j1] = 1.0;
            }
        }

        // create new basis function
        typename PBBSplinesBasisFunctionType::Pointer pnew_bf = pNewFESpace->CreateBf(id, knot_vecs);

        // create the cells for the basis function
        for (std::size_t i1 = 0; i1 < knot_vecs[0].size() - 1; ++i1)
        {
            Real xi_min = knot_vecs[0][i1];
            Real xi_max = knot_vecs[0][i1 + 1];

            for (std::size_t j1 = 0; j1 < knot_vecs[1].size() - 1; ++j1)
            {
                Real eta_min = knot_vecs[1][j1];
                Real eta_max = knot_vecs[1][j1 + 1];

                // check if the cell domain area is nonzero
                double area = (xi_max - xi_min) * (eta_max - eta_min);
                if (fabs(area) > area_tol)
                {
                    std::vector<Real> knots = {xi_min, xi_max, eta_min, eta_max};
                    cell_t p_cell = pNewFESpace->pCellManager()->CreateCell(knots);
                    pnew_bf->AddCell(p_cell);
                }
            }
        }

        // transfer the control point
        x = ptpoint->getX();
        y = ptpoint->getY();
        z = ptpoint->getZ();
        w = ptpoint->getW();
        ControlPointType c(w * x, w * y, w * z, w);
        pnew_bf->SetValue(CONTROL_POINT, c);
    }

    /** create faces **/
    TImagePtr timage = tspline->getTImage();
    // WATCH(timage->sizeFaces())
    double xi_min, xi_max, eta_min, eta_max;
    for (TFacVIterator fit = timage->faceIteratorBegin(); fit != timage->faceIteratorEnd(); ++fit)
    {
        xi_min = ((*fit)->northWest().s() - min_xi_min) / (max_xi_max - min_xi_min);
        xi_max = ((*fit)->southEast().s() - min_xi_min) / (max_xi_max - min_xi_min);
        eta_min = ((*fit)->southEast().t() - min_eta_min) / (max_eta_max - min_eta_min);
        eta_max = ((*fit)->northWest().t() - min_eta_min) / (max_eta_max - min_eta_min);
        std::vector<double> knots = {xi_min, xi_max, eta_min, eta_max};
        cell_t p_cell = pNewFESpace->pFaceManager()->CreateCell(knots);

        // std::cout << "checking cell " << xi_min << " " << xi_max << " " << eta_min << " " << eta_max << std::endl;

        // // adding the cell information to the basis functions
        // for (TNodVIterator nit = (*fit)->blendingNodeIteratorBegin(); nit != (*fit)->blendingNodeIteratorEnd(); ++nit)
        // {
        //     id = (*nit)->getTPoint()->getId();
        //     typename PBBSplinesBasisFunctionType::Pointer p_bf = pNewFESpace->operator()(id);
        //     p_bf->AddCell(p_cell);
        // }
    }

    // create and return patch
    std::size_t patch_id = 1;
    Patch<2>::Pointer pNewPatch = Patch<2>::Create(patch_id, pNewFESpace);

    // set control points grid function
    typename ControlGrid<ControlPointType>::Pointer pControlPointGrid
        = ControlGridUtility::CreatePointBasedControlGrid<ControlPointType, TSplinesFESpaceType>(CONTROL_POINT, pNewFESpace);
    pNewPatch->CreateControlPointGridFunction(pControlPointGrid);

    std::cout << "Import T-splines from " << filename << " successfully" << std::endl;

    return pNewPatch;
}

MultiPatch<2>::Pointer TSplinesPatchTSMImporter::Import(const std::string& filename) const
{
    KRATOS_ERROR << "TSM reader does not support multipatch";
}

} // namespace Kratos.
