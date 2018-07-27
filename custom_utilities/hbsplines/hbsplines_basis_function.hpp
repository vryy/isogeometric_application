//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 24 Nov 2017 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_BASIS_FUNCTION_HPP_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_BASIS_FUNCTION_HPP_INCLUDED

#include "custom_utilities/bezier_utils.h"
#include "custom_utilities/bspline_utils.h"

namespace Kratos
{

template<>
template<typename TVectorType, typename TIArrayType, typename TKnotContainerType, class TCellType>
inline void HBSplinesBasisFunction_Helper<1>::ComputeExtractionOperator(TVectorType& Crow,
        const TIArrayType& orders, const TKnotContainerType& local_knots, const TCellType& r_cell)
{
    std::stringstream ss;
    ss << __FUNCTION__ << " is not yet implemented for 1D";
    KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
}

template<>
template<typename TVectorType, typename TIArrayType, typename TKnotContainerType, class TCellType>
inline void HBSplinesBasisFunction_Helper<2>::ComputeExtractionOperator(TVectorType& Crow,
        const TIArrayType& orders, const TKnotContainerType& local_knots, const TCellType& r_cell)
{
    // a priori check
    for(std::size_t i = 0; i < local_knots[0].size(); ++i)
    {
        if(local_knots[0][i] > r_cell.LeftValue() && local_knots[0][i] < r_cell.RightValue())
        {
            std::stringstream ss;
            ss << "Error: the cell is not contained in one knot span in u-direction of the basis function" << std::endl;
            ss << "LeftValue: " << r_cell.LeftValue() << ", RightValue: " << r_cell.RightValue() << std::endl;
            ss << "local_knots[0]:";
            for(std::size_t j = 0; j < local_knots[0].size(); ++j)
                ss << " " << local_knots[0][j];
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
    }
    for(std::size_t i = 0; i < local_knots[1].size(); ++i)
    {
        if(local_knots[1][i] > r_cell.DownValue() && local_knots[1][i] < r_cell.UpValue())
        {
            std::stringstream ss;
            ss << "Error: the cell is not contained in one knot span in v-direction of the basis function" << std::endl;
            ss << "DownValue: " << r_cell.DownValue() << ", UpValue: " << r_cell.UpValue() << std::endl;
            ss << "local_knots[1]:";
            for(std::size_t j = 0; j < local_knots[1].size(); ++j)
                ss << " " << local_knots[1][j];
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
    }

    #ifdef DEBUG_BEZIER_EXTRACTION
    std::cout << "Bezier extraction debug at bf " << Id() << ":" << std::endl;
    std::cout << "local_knots[0]:";
    for(std::size_t i = 0; i < local_knots[0].size(); ++i)
        std::cout << " " << local_knots[0][i];
    std::cout << std::endl;
    std::cout << "local_knots[1]:";
    for(std::size_t i = 0; i < local_knots[1].size(); ++i)
        std::cout << " " << local_knots[1][i];
    std::cout << std::endl;
    KRATOS_WATCH(r_cell.LeftValue())
    KRATOS_WATCH(r_cell.RightValue())
    KRATOS_WATCH(r_cell.DownValue())
    KRATOS_WATCH(r_cell.UpValue())
    #endif

    // compute the inserted knot vector
    std::vector<double> ins_knots1;
//        std::vector<int> ins_span1;
    for(std::size_t i = 0; i < local_knots[0].size() - 1; ++i)
    {
        if(r_cell.LeftValue() > local_knots[0][i] && r_cell.LeftValue() < local_knots[0][i+1])
        {
            ins_knots1.push_back(r_cell.LeftValue());
//                ins_span1.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }
    }

    for(std::size_t i = 0; i < local_knots[0].size() - 1; ++i)
    {
        if(r_cell.RightValue() > local_knots[0][i] && r_cell.RightValue() < local_knots[0][i+1])
        {
            ins_knots1.push_back(r_cell.RightValue());
//                ins_span1.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }
    }

    std::vector<double> ins_knots2;
//        std::vector<int> ins_span2;
    for(std::size_t i = 0; i < local_knots[1].size() - 1; ++i)
    {
        if(r_cell.DownValue() > local_knots[1][i] && r_cell.DownValue() < local_knots[1][i+1])
        {
            ins_knots2.push_back(r_cell.DownValue());
//                ins_span2.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }
    }

    for(std::size_t i = 0; i < local_knots[1].size() - 1; ++i)
    {
        if(r_cell.UpValue() > local_knots[1][i] && r_cell.UpValue() < local_knots[1][i+1])
        {
            ins_knots2.push_back(r_cell.UpValue());
//                ins_span2.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }
    }

    #ifdef DEBUG_BEZIER_EXTRACTION
    std::cout << "ins_knots1:";
    for(std::size_t i = 0; i < ins_knots1.size(); ++i)
        std::cout << " " << ins_knots1[i];
    std::cout << std::endl;
//        std::cout << "ins_span1:";
//        for(std::size_t i = 0; i < ins_span1.size(); ++i)
//            std::cout << " " << ins_span1[i];
//        std::cout << std::endl;

    std::cout << "ins_knots2:";
    for(std::size_t i = 0; i < ins_knots2.size(); ++i)
        std::cout << " " << ins_knots2[i];
    std::cout << std::endl;
//        std::cout << "ins_span2:";
//        for(std::size_t i = 0; i < ins_span2.size(); ++i)
//            std::cout << " " << ins_span2[i];
//        std::cout << std::endl;
    #endif

    // compute the Bezier extraction operator
    std::vector<Vector> Crows;
    int nb_xi, nb_eta;
    Vector Ubar_xi, Ubar_eta;

//        BezierUtils::bezier_extraction_tsplines_2d(Crows,
//                                                   nb_xi,
//                                                   nb_eta,
//                                                   Ubar_xi,
//                                                   Ubar_eta,
//                                                   local_knots[0],
//                                                   local_knots[1],
//                                                   ins_knots1,
//                                                   ins_knots2,
//                                                   ins_span1,
//                                                   ins_span2,
//                                                   orders[0],
//                                                   orders[1]);

    BezierUtils::bezier_extraction_local_2d(Crows,
                                            nb_xi,
                                            nb_eta,
                                            Ubar_xi,
                                            Ubar_eta,
                                            local_knots[0],
                                            local_knots[1],
                                            ins_knots1,
                                            ins_knots2,
                                            orders[0],
                                            orders[1]);

    #ifdef DEBUG_BEZIER_EXTRACTION
    std::cout << "Crows:" << std::endl;
    for(std::size_t i = 0; i < Crows.size(); ++i)
        std::cout << Crows[i] << std::endl;
    #endif

    // extract the correct row
    std::size_t span1;// knot span in u-direction of the cell w.r.t basis function support
    std::size_t span2; // knot span in v-direction of the cell w.r.t basis function support
    #ifdef DEBUG_BEZIER_EXTRACTION
    KRATOS_WATCH(Ubar_xi)
    KRATOS_WATCH(Ubar_eta)
    #endif
    std::set<double> Ubar_xi_unique(Ubar_xi.begin(), Ubar_xi.end());
    std::set<double> Ubar_eta_unique(Ubar_eta.begin(), Ubar_eta.end());
    std::vector<double> Ubar_xi_unique_vector(Ubar_xi_unique.begin(), Ubar_xi_unique.end());
    std::vector<double> Ubar_eta_unique_vector(Ubar_eta_unique.begin(), Ubar_eta_unique.end());
    span1 = BSplineUtils::FindSpanLocal(r_cell.LeftValue(), Ubar_xi_unique_vector) - 1;
    span2 = BSplineUtils::FindSpanLocal(r_cell.DownValue(), Ubar_eta_unique_vector) - 1;

    #ifdef DEBUG_BEZIER_EXTRACTION
    KRATOS_WATCH(span1)
    KRATOS_WATCH(span2)
    #endif
//        std::size_t span = span2 * nb_xi + span1;
    std::size_t span = span1 * nb_eta + span2;
    if(Crow.size() != Crows[span].size())
        Crow.resize(Crows[span].size());
    #ifdef DEBUG_BEZIER_EXTRACTION
    KRATOS_WATCH(span)
    KRATOS_WATCH(Crows.size())
    #endif
    std::copy(Crows[span].begin(), Crows[span].end(), Crow.begin());
    #ifdef DEBUG_BEZIER_EXTRACTION
    std::cout << "----------------------" << std::endl;
    #endif
}


template<>
template<typename TVectorType, typename TIArrayType, typename TKnotContainerType, class TCellType>
inline void HBSplinesBasisFunction_Helper<3>::ComputeExtractionOperator(TVectorType& Crow,
        const TIArrayType& orders, const TKnotContainerType& local_knots, const TCellType& r_cell)
{
    // a priori check
    for(std::size_t i = 0; i < local_knots[0].size(); ++i)
    {
        if(local_knots[0][i] > r_cell.LeftValue() && local_knots[0][i] < r_cell.RightValue())
        {
            std::stringstream ss;
            ss << "Error: the cell is not contained in one knot span in u-direction of the basis function" << std::endl;
            ss << "LeftValue: " << r_cell.LeftValue() << ", RightValue: " << r_cell.RightValue() << std::endl;
            ss << "local_knots[0]:";
            for(std::size_t j = 0; j < local_knots[0].size(); ++j)
                ss << " " << local_knots[0][j];
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
    }

    for(std::size_t i = 0; i < local_knots[1].size(); ++i)
    {
        if(local_knots[1][i] > r_cell.DownValue() && local_knots[1][i] < r_cell.UpValue())
        {
            std::stringstream ss;
            ss << "Error: the cell is not contained in one knot span in v-direction of the basis function" << std::endl;
            ss << "DownValue: " << r_cell.DownValue() << ", UpValue: " << r_cell.UpValue() << std::endl;
            ss << "local_knots[1]:";
            for(std::size_t j = 0; j < local_knots[1].size(); ++j)
                ss << " " << local_knots[1][j];
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
    }

    for(std::size_t i = 0; i < local_knots[2].size(); ++i)
    {
        if(local_knots[2][i] > r_cell.BelowValue() && local_knots[2][i] < r_cell.AboveValue())
        {
            std::stringstream ss;
            ss << "Error: the cell is not contained in one knot span in w-direction of the basis function" << std::endl;
            ss << "BelowValue: " << r_cell.BelowValue() << ", AboveValue: " << r_cell.AboveValue() << std::endl;
            ss << "local_knots[2]:";
            for(std::size_t j = 0; j < local_knots[2].size(); ++j)
                ss << " " << local_knots[2][j];
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
    }

    #ifdef DEBUG_BEZIER_EXTRACTION
    std::cout << "Bezier extraction debug at bf " << Id() << ":" << std::endl;
    std::cout << "local_knots[0]:";
    for(std::size_t i = 0; i < local_knots[0].size(); ++i)
        std::cout << " " << local_knots[0][i];
    std::cout << std::endl;
    std::cout << "local_knots[1]:";
    for(std::size_t i = 0; i < local_knots[1].size(); ++i)
        std::cout << " " << local_knots[1][i];
    std::cout << std::endl;
    std::cout << "local_knots[2]:";
    for(std::size_t i = 0; i < local_knots[2].size(); ++i)
        std::cout << " " << local_knots[2][i];
    std::cout << std::endl;
    KRATOS_WATCH(r_cell.LeftValue())
    KRATOS_WATCH(r_cell.RightValue())
    KRATOS_WATCH(r_cell.DownValue())
    KRATOS_WATCH(r_cell.UpValue())
    KRATOS_WATCH(r_cell.BelowValue())
    KRATOS_WATCH(r_cell.AboveValue())
    #endif

    // compute the inserted knot vector
    std::vector<double> ins_knots1;
//        std::vector<int> ins_span1;
    for(std::size_t i = 0; i < local_knots[0].size() - 1; ++i)
    {
        if(r_cell.LeftValue() > local_knots[0][i] && r_cell.LeftValue() < local_knots[0][i+1])
        {
            ins_knots1.push_back(r_cell.LeftValue());
//                ins_span1.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }
    }
    for(std::size_t i = 0; i < local_knots[0].size() - 1; ++i)
    {
        if(r_cell.RightValue() > local_knots[0][i] && r_cell.RightValue() < local_knots[0][i+1])
        {
            ins_knots1.push_back(r_cell.RightValue());
//                ins_span1.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }
    }

    std::vector<double> ins_knots2;
//        std::vector<int> ins_span2;
    for(std::size_t i = 0; i < local_knots[1].size() - 1; ++i)
    {
        if(r_cell.DownValue() > local_knots[1][i] && r_cell.DownValue() < local_knots[1][i+1])
        {
            ins_knots2.push_back(r_cell.DownValue());
//                ins_span2.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }
    }
    for(std::size_t i = 0; i < local_knots[1].size() - 1; ++i)
    {
        if(r_cell.UpValue() > local_knots[1][i] && r_cell.UpValue() < local_knots[1][i+1])
        {
            ins_knots2.push_back(r_cell.UpValue());
//                ins_span2.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }
    }

    std::vector<double> ins_knots3;
//        std::vector<int> ins_span3;
    for(std::size_t i = 0; i < local_knots[2].size() - 1; ++i)
    {
        if(r_cell.BelowValue() > local_knots[2][i] && r_cell.BelowValue() < local_knots[2][i+1])
        {
            ins_knots3.push_back(r_cell.BelowValue());
//                ins_span3.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }
    }
    for(std::size_t i = 0; i < local_knots[2].size() - 1; ++i)
    {
        if(r_cell.AboveValue() > local_knots[2][i] && r_cell.AboveValue() < local_knots[2][i+1])
        {
            ins_knots3.push_back(r_cell.AboveValue());
//                ins_span3.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }
    }

    #ifdef DEBUG_BEZIER_EXTRACTION
    std::cout << "ins_knots1:";
    for(std::size_t i = 0; i < ins_knots1.size(); ++i)
        std::cout << " " << ins_knots1[i];
    std::cout << std::endl;
//        std::cout << "ins_span1:";
//        for(std::size_t i = 0; i < ins_span1.size(); ++i)
//            std::cout << " " << ins_span1[i];
//        std::cout << std::endl;

    std::cout << "ins_knots2:";
    for(std::size_t i = 0; i < ins_knots2.size(); ++i)
        std::cout << " " << ins_knots2[i];
    std::cout << std::endl;
//        std::cout << "ins_span2:";
//        for(std::size_t i = 0; i < ins_span2.size(); ++i)
//            std::cout << " " << ins_span2[i];
//        std::cout << std::endl;

    std::cout << "ins_knots3:";
    for(std::size_t i = 0; i < ins_knots3.size(); ++i)
        std::cout << " " << ins_knots3[i];
    std::cout << std::endl;
//        std::cout << "ins_span3:";
//        for(std::size_t i = 0; i < ins_span3.size(); ++i)
//            std::cout << " " << ins_span3[i];
//        std::cout << std::endl;
    #endif

    // compute the Bezier extraction operator
    std::vector<Vector> Crows;
    int nb_xi, nb_eta, nb_zeta;
    Vector Ubar_xi, Ubar_eta, Ubar_zeta;

//        BezierUtils::bezier_extraction_tsplines_3d(Crows,
//                                                   nb_xi,
//                                                   nb_eta,
//                                                   nb_zeta,
//                                                   Ubar_xi,
//                                                   Ubar_eta,
//                                                   Ubar_zeta,
//                                                   local_knots[0],
//                                                   local_knots[1],
//                                                   local_knots[2],
//                                                   ins_knots1,
//                                                   ins_knots2,
//                                                   ins_knots3,
//                                                   ins_span1,
//                                                   ins_span2,
//                                                   ins_span3,
//                                                   orders[0],
//                                                   orders[1]);

    BezierUtils::bezier_extraction_local_3d(Crows,
                                            nb_xi,
                                            nb_eta,
                                            nb_zeta,
                                            Ubar_xi,
                                            Ubar_eta,
                                            Ubar_zeta,
                                            local_knots[0],
                                            local_knots[1],
                                            local_knots[2],
                                            ins_knots1,
                                            ins_knots2,
                                            ins_knots3,
                                            orders[0],
                                            orders[1],
                                            orders[2]);

    #ifdef DEBUG_BEZIER_EXTRACTION
    std::cout << "Crows:" << std::endl;
    for(std::size_t i = 0; i < Crows.size(); ++i)
        std::cout << Crows[i] << std::endl;
    #endif

    // extract the correct row
    std::size_t span1; // knot span in u-direction of the cell w.r.t basis function support
    std::size_t span2; // knot span in v-direction of the cell w.r.t basis function support
    std::size_t span3; // knot span in w-direction of the cell w.r.t basis function support
    #ifdef DEBUG_BEZIER_EXTRACTION
    KRATOS_WATCH(Ubar_xi)
    KRATOS_WATCH(Ubar_eta)
    KRATOS_WATCH(Ubar_zeta)
    #endif
    std::set<double> Ubar_xi_unique(Ubar_xi.begin(), Ubar_xi.end());
    std::set<double> Ubar_eta_unique(Ubar_eta.begin(), Ubar_eta.end());
    std::set<double> Ubar_zeta_unique(Ubar_zeta.begin(), Ubar_zeta.end());
    std::vector<double> Ubar_xi_unique_vector(Ubar_xi_unique.begin(), Ubar_xi_unique.end());
    std::vector<double> Ubar_eta_unique_vector(Ubar_eta_unique.begin(), Ubar_eta_unique.end());
    std::vector<double> Ubar_zeta_unique_vector(Ubar_zeta_unique.begin(), Ubar_zeta_unique.end());
    span1 = BSplineUtils::FindSpanLocal(r_cell.LeftValue(), Ubar_xi_unique_vector) - 1;
    span2 = BSplineUtils::FindSpanLocal(r_cell.DownValue(), Ubar_eta_unique_vector) - 1;
    span3 = BSplineUtils::FindSpanLocal(r_cell.BelowValue(), Ubar_zeta_unique_vector) - 1;

    #ifdef DEBUG_BEZIER_EXTRACTION
    KRATOS_WATCH(span1)
    KRATOS_WATCH(span2)
    KRATOS_WATCH(span3)
    #endif
    std::size_t span = (span1 * nb_eta + span2) * nb_zeta + span3;
    if(Crow.size() != Crows[span].size())
        Crow.resize(Crows[span].size());
    #ifdef DEBUG_BEZIER_EXTRACTION
    KRATOS_WATCH(span)
    KRATOS_WATCH(Crows.size())
    #endif
    std::copy(Crows[span].begin(), Crows[span].end(), Crow.begin());
    #ifdef DEBUG_BEZIER_EXTRACTION
    std::cout << "----------------------" << std::endl;
    #endif
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_HBSPLINES_BASIS_FUNCTION_HPP_INCLUDED

