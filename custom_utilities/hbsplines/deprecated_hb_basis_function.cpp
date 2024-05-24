#include "custom_utilities/hbsplines/deprecated_hb_basis_function.h"

// #define DEBUG_BEZIER_EXTRACTION

namespace Kratos
{
void DeprecatedHBBasisFunction::ComputeExtractionOperator(cell_t p_cell, Vector& Crow, int p1, int p2)
{
    // a priori check
    for (std::size_t i = 0; i < mpLocalKnots1.size(); ++i)
        if (mpLocalKnots1[i]->Value() > p_cell->XiMinValue() && mpLocalKnots1[i]->Value() < p_cell->XiMaxValue())
        {
            std::stringstream ss;
            ss << "Error: the cell is not contained in one knot span in u-direction of the basis function" << std::endl;
            ss << "XiMinValue: " << p_cell->XiMinValue() << ", XiMaxValue: " << p_cell->XiMaxValue() << std::endl;
            ss << "mpLocalKnots1:";
            for (std::size_t j = 0; j < mpLocalKnots1.size(); ++j)
            {
                ss << " " << mpLocalKnots1[j];
            }
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
    for (std::size_t i = 0; i < mpLocalKnots2.size(); ++i)
        if (mpLocalKnots2[i]->Value() > p_cell->EtaMinValue() && mpLocalKnots2[i]->Value() < p_cell->EtaMaxValue())
        {
            std::stringstream ss;
            ss << "Error: the cell is not contained in one knot span in v-direction of the basis function" << std::endl;
            ss << "EtaMinValue: " << p_cell->EtaMinValue() << ", EtaMaxValue: " << p_cell->EtaMaxValue() << std::endl;
            ss << "mpLocalKnots2:";
            for (std::size_t j = 0; j < mpLocalKnots2.size(); ++j)
            {
                ss << " " << mpLocalKnots2[j];
            }
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }

#ifdef DEBUG_BEZIER_EXTRACTION
    std::cout << "Bezier extraction debug at bf " << Id() << ":" << std::endl;
    std::cout << "mpLocalKnots1:";
    for (std::size_t i = 0; i < mpLocalKnots1.size(); ++i)
    {
        std::cout << " " << mpLocalKnots1[i]->Value();
    }
    std::cout << std::endl;
    std::cout << "mpLocalKnots2:";
    for (std::size_t i = 0; i < mpLocalKnots2.size(); ++i)
    {
        std::cout << " " << mpLocalKnots2[i]->Value();
    }
    std::cout << std::endl;
    KRATOS_WATCH(p_cell->XiMinValue())
    KRATOS_WATCH(p_cell->XiMaxValue())
    KRATOS_WATCH(p_cell->EtaMinValue())
    KRATOS_WATCH(p_cell->EtaMaxValue())
#endif

    // compute the inserted knot vector
    std::vector<double> ins_knots1;
//        std::vector<int> ins_span1;
    for (std::size_t i = 0; i < mpLocalKnots1.size() - 1; ++i)
        if (p_cell->XiMinValue() > mpLocalKnots1[i]->Value() && p_cell->XiMinValue() < mpLocalKnots1[i + 1]->Value())
        {
            ins_knots1.push_back(p_cell->XiMinValue());
//                ins_span1.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }
    for (std::size_t i = 0; i < mpLocalKnots1.size() - 1; ++i)
        if (p_cell->XiMaxValue() > mpLocalKnots1[i]->Value() && p_cell->XiMaxValue() < mpLocalKnots1[i + 1]->Value())
        {
            ins_knots1.push_back(p_cell->XiMaxValue());
//                ins_span1.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }

    std::vector<double> ins_knots2;
//        std::vector<int> ins_span2;
    for (std::size_t i = 0; i < mpLocalKnots2.size() - 1; ++i)
        if (p_cell->EtaMinValue() > mpLocalKnots2[i]->Value() && p_cell->EtaMinValue() < mpLocalKnots2[i + 1]->Value())
        {
            ins_knots2.push_back(p_cell->EtaMinValue());
//                ins_span2.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }
    for (std::size_t i = 0; i < mpLocalKnots2.size() - 1; ++i)
        if (p_cell->EtaMaxValue() > mpLocalKnots2[i]->Value() && p_cell->EtaMaxValue() < mpLocalKnots2[i + 1]->Value())
        {
            ins_knots2.push_back(p_cell->EtaMaxValue());
//                ins_span2.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }

#ifdef DEBUG_BEZIER_EXTRACTION
    std::cout << "ins_knots1:";
    for (std::size_t i = 0; i < ins_knots1.size(); ++i)
    {
        std::cout << " " << ins_knots1[i];
    }
    std::cout << std::endl;
//        std::cout << "ins_span1:";
//        for(std::size_t i = 0; i < ins_span1.size(); ++i)
//            std::cout << " " << ins_span1[i];
//        std::cout << std::endl;

    std::cout << "ins_knots2:";
    for (std::size_t i = 0; i < ins_knots2.size(); ++i)
    {
        std::cout << " " << ins_knots2[i];
    }
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
    std::vector<double> LocalKnots1;
    std::vector<double> LocalKnots2;
    std::vector<double> LocalKnots3;
    this->GetLocalKnots(1, LocalKnots1);
    this->GetLocalKnots(2, LocalKnots2);

//        BezierUtils::bezier_extraction_tsplines_2d(Crows,
//                                                   nb_xi,
//                                                   nb_eta,
//                                                   Ubar_xi,
//                                                   Ubar_eta,
//                                                   LocalKnots1,
//                                                   LocalKnots2,
//                                                   ins_knots1,
//                                                   ins_knots2,
//                                                   ins_span1,
//                                                   ins_span2,
//                                                   p1,
//                                                   p2);

    BezierUtils::bezier_extraction_local_2d(Crows,
                                            nb_xi,
                                            nb_eta,
                                            Ubar_xi,
                                            Ubar_eta,
                                            LocalKnots1,
                                            LocalKnots2,
                                            ins_knots1,
                                            ins_knots2,
                                            p1,
                                            p2);

#ifdef DEBUG_BEZIER_EXTRACTION
    std::cout << "Crows:" << std::endl;
    for (std::size_t i = 0; i < Crows.size(); ++i)
    {
        std::cout << Crows[i] << std::endl;
    }
#endif

    // extract the correct row
    int span1;// knot span in u-direction of the cell w.r.t basis function support
    int span2; // knot span in v-direction of the cell w.r.t basis function support
#ifdef DEBUG_BEZIER_EXTRACTION
    KRATOS_WATCH(Ubar_xi)
    KRATOS_WATCH(Ubar_eta)
#endif
    std::set<double> Ubar_xi_unique(Ubar_xi.begin(), Ubar_xi.end());
    std::set<double> Ubar_eta_unique(Ubar_eta.begin(), Ubar_eta.end());
    std::vector<double> Ubar_xi_unique_vector(Ubar_xi_unique.begin(), Ubar_xi_unique.end());
    std::vector<double> Ubar_eta_unique_vector(Ubar_eta_unique.begin(), Ubar_eta_unique.end());
    span1 = this->FindSpanLocal(p_cell->XiMinValue(), Ubar_xi_unique_vector) - 1;
    span2 = this->FindSpanLocal(p_cell->EtaMinValue(), Ubar_eta_unique_vector) - 1;

#ifdef DEBUG_BEZIER_EXTRACTION
    KRATOS_WATCH(span1)
    KRATOS_WATCH(span2)
#endif
//        int span = span2 * nb_xi + span1;
    int span = span1 * nb_eta + span2;
    if (Crow.size() != Crows[span].size())
    {
        Crow.resize(Crows[span].size());
    }
#ifdef DEBUG_BEZIER_EXTRACTION
    KRATOS_WATCH(span)
    KRATOS_WATCH(Crows.size())
#endif
    std::copy(Crows[span].begin(), Crows[span].end(), Crow.begin());
#ifdef DEBUG_BEZIER_EXTRACTION
    std::cout << "----------------------" << std::endl;
#endif
}

void DeprecatedHBBasisFunction::ComputeExtractionOperator(cell_t p_cell, Vector& Crow, int p1, int p2, int p3)
{
    // a priori check
    for (std::size_t i = 0; i < mpLocalKnots1.size(); ++i)
        if (mpLocalKnots1[i]->Value() > p_cell->XiMinValue() && mpLocalKnots1[i]->Value() < p_cell->XiMaxValue())
        {
            std::stringstream ss;
            ss << "Error: the cell is not contained in one knot span in u-direction of the basis function" << std::endl;
            ss << "XiMinValue: " << p_cell->XiMinValue() << ", XiMaxValue: " << p_cell->XiMaxValue() << std::endl;
            ss << "mpLocalKnots1:";
            for (std::size_t j = 0; j < mpLocalKnots1.size(); ++j)
            {
                ss << " " << mpLocalKnots1[j];
            }
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
    for (std::size_t i = 0; i < mpLocalKnots2.size(); ++i)
        if (mpLocalKnots2[i]->Value() > p_cell->EtaMinValue() && mpLocalKnots2[i]->Value() < p_cell->EtaMaxValue())
        {
            std::stringstream ss;
            ss << "Error: the cell is not contained in one knot span in v-direction of the basis function" << std::endl;
            ss << "EtaMinValue: " << p_cell->EtaMinValue() << ", EtaMaxValue: " << p_cell->EtaMaxValue() << std::endl;
            ss << "mpLocalKnots2:";
            for (std::size_t j = 0; j < mpLocalKnots2.size(); ++j)
            {
                ss << " " << mpLocalKnots2[j];
            }
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
    for (std::size_t i = 0; i < mpLocalKnots3.size(); ++i)
        if (mpLocalKnots3[i]->Value() > p_cell->ZetaMinValue() && mpLocalKnots3[i]->Value() < p_cell->ZetaMaxValue())
        {
            std::stringstream ss;
            ss << "Error: the cell is not contained in one knot span in w-direction of the basis function" << std::endl;
            ss << "ZetaMinValue: " << p_cell->ZetaMinValue() << ", ZetaMaxValue: " << p_cell->ZetaMaxValue() << std::endl;
            ss << "mpLocalKnots3:";
            for (std::size_t j = 0; j < mpLocalKnots3.size(); ++j)
            {
                ss << " " << mpLocalKnots3[j];
            }
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }

#ifdef DEBUG_BEZIER_EXTRACTION
    std::cout << "Bezier extraction debug at bf " << Id() << ":" << std::endl;
    std::cout << "mpLocalKnots1:";
    for (std::size_t i = 0; i < mpLocalKnots1.size(); ++i)
    {
        std::cout << " " << mpLocalKnots1[i]->Value();
    }
    std::cout << std::endl;
    std::cout << "mpLocalKnots2:";
    for (std::size_t i = 0; i < mpLocalKnots2.size(); ++i)
    {
        std::cout << " " << mpLocalKnots2[i]->Value();
    }
    std::cout << std::endl;
    std::cout << "mpLocalKnots3:";
    for (std::size_t i = 0; i < mpLocalKnots3.size(); ++i)
    {
        std::cout << " " << mpLocalKnots3[i]->Value();
    }
    std::cout << std::endl;
    KRATOS_WATCH(p_cell->XiMinValue())
    KRATOS_WATCH(p_cell->XiMaxValue())
    KRATOS_WATCH(p_cell->EtaMinValue())
    KRATOS_WATCH(p_cell->EtaMaxValue())
    KRATOS_WATCH(p_cell->ZetaMinValue())
    KRATOS_WATCH(p_cell->ZetaMaxValue())
#endif

    // compute the inserted knot vector
    std::vector<double> ins_knots1;
//        std::vector<int> ins_span1;
    for (std::size_t i = 0; i < mpLocalKnots1.size() - 1; ++i)
        if (p_cell->XiMinValue() > mpLocalKnots1[i]->Value() && p_cell->XiMinValue() < mpLocalKnots1[i + 1]->Value())
        {
            ins_knots1.push_back(p_cell->XiMinValue());
//                ins_span1.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }
    for (std::size_t i = 0; i < mpLocalKnots1.size() - 1; ++i)
        if (p_cell->XiMaxValue() > mpLocalKnots1[i]->Value() && p_cell->XiMaxValue() < mpLocalKnots1[i + 1]->Value())
        {
            ins_knots1.push_back(p_cell->XiMaxValue());
//                ins_span1.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }

    std::vector<double> ins_knots2;
//        std::vector<int> ins_span2;
    for (std::size_t i = 0; i < mpLocalKnots2.size() - 1; ++i)
        if (p_cell->EtaMinValue() > mpLocalKnots2[i]->Value() && p_cell->EtaMinValue() < mpLocalKnots2[i + 1]->Value())
        {
            ins_knots2.push_back(p_cell->EtaMinValue());
//                ins_span2.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }
    for (std::size_t i = 0; i < mpLocalKnots2.size() - 1; ++i)
        if (p_cell->EtaMaxValue() > mpLocalKnots2[i]->Value() && p_cell->EtaMaxValue() < mpLocalKnots2[i + 1]->Value())
        {
            ins_knots2.push_back(p_cell->EtaMaxValue());
//                ins_span2.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }

    std::vector<double> ins_knots3;
//        std::vector<int> ins_span3;
    for (std::size_t i = 0; i < mpLocalKnots3.size() - 1; ++i)
        if (p_cell->ZetaMinValue() > mpLocalKnots3[i]->Value() && p_cell->ZetaMinValue() < mpLocalKnots3[i + 1]->Value())
        {
            ins_knots3.push_back(p_cell->ZetaMinValue());
//                ins_span3.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }
    for (std::size_t i = 0; i < mpLocalKnots3.size() - 1; ++i)
        if (p_cell->ZetaMaxValue() > mpLocalKnots3[i]->Value() && p_cell->ZetaMaxValue() < mpLocalKnots3[i + 1]->Value())
        {
            ins_knots3.push_back(p_cell->ZetaMaxValue());
//                ins_span3.push_back(i+1); // +1 because the bezier_extraction_tsplines_2d takes based-1 index
            break;
        }

#ifdef DEBUG_BEZIER_EXTRACTION
    std::cout << "ins_knots1:";
    for (std::size_t i = 0; i < ins_knots1.size(); ++i)
    {
        std::cout << " " << ins_knots1[i];
    }
    std::cout << std::endl;
//        std::cout << "ins_span1:";
//        for(std::size_t i = 0; i < ins_span1.size(); ++i)
//            std::cout << " " << ins_span1[i];
//        std::cout << std::endl;

    std::cout << "ins_knots2:";
    for (std::size_t i = 0; i < ins_knots2.size(); ++i)
    {
        std::cout << " " << ins_knots2[i];
    }
    std::cout << std::endl;
//        std::cout << "ins_span2:";
//        for(std::size_t i = 0; i < ins_span2.size(); ++i)
//            std::cout << " " << ins_span2[i];
//        std::cout << std::endl;

    std::cout << "ins_knots3:";
    for (std::size_t i = 0; i < ins_knots3.size(); ++i)
    {
        std::cout << " " << ins_knots3[i];
    }
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
    std::vector<double> LocalKnots1;
    std::vector<double> LocalKnots2;
    std::vector<double> LocalKnots3;
    this->GetLocalKnots(1, LocalKnots1);
    this->GetLocalKnots(2, LocalKnots2);
    this->GetLocalKnots(3, LocalKnots3);

//        BezierUtils::bezier_extraction_tsplines_2d(Crows,
//                                                   nb_xi,
//                                                   nb_eta,
//                                                   Ubar_xi,
//                                                   Ubar_eta,
//                                                   LocalKnots1,
//                                                   LocalKnots2,
//                                                   ins_knots1,
//                                                   ins_knots2,
//                                                   ins_span1,
//                                                   ins_span2,
//                                                   p1,
//                                                   p2);

    BezierUtils::bezier_extraction_local_3d(Crows,
                                            nb_xi,
                                            nb_eta,
                                            nb_zeta,
                                            Ubar_xi,
                                            Ubar_eta,
                                            Ubar_zeta,
                                            LocalKnots1,
                                            LocalKnots2,
                                            LocalKnots3,
                                            ins_knots1,
                                            ins_knots2,
                                            ins_knots3,
                                            p1,
                                            p2,
                                            p3);

#ifdef DEBUG_BEZIER_EXTRACTION
    std::cout << "Crows:" << std::endl;
    for (std::size_t i = 0; i < Crows.size(); ++i)
    {
        std::cout << Crows[i] << std::endl;
    }
#endif

    // extract the correct row
    int span1; // knot span in u-direction of the cell w.r.t basis function support
    int span2; // knot span in v-direction of the cell w.r.t basis function support
    int span3; // knot span in w-direction of the cell w.r.t basis function support
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
    span1 = this->FindSpanLocal(p_cell->XiMinValue(), Ubar_xi_unique_vector) - 1;
    span2 = this->FindSpanLocal(p_cell->EtaMinValue(), Ubar_eta_unique_vector) - 1;
    span3 = this->FindSpanLocal(p_cell->ZetaMinValue(), Ubar_zeta_unique_vector) - 1;

#ifdef DEBUG_BEZIER_EXTRACTION
    KRATOS_WATCH(span1)
    KRATOS_WATCH(span2)
    KRATOS_WATCH(span3)
#endif
    int span = (span1 * nb_eta + span2) * nb_zeta + span3;
    if (Crow.size() != Crows[span].size())
    {
        Crow.resize(Crows[span].size());
    }
#ifdef DEBUG_BEZIER_EXTRACTION
    KRATOS_WATCH(span)
    KRATOS_WATCH(Crows.size())
#endif
    std::copy(Crows[span].begin(), Crows[span].end(), Crow.begin());
#ifdef DEBUG_BEZIER_EXTRACTION
    std::cout << "----------------------" << std::endl;
#endif
}
}

#undef DEBUG_BEZIER_EXTRACTION
