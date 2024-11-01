#include "includes/define.h"
#include "custom_utilities/bezier_utils.h"
#include "custom_utilities/isogeometric_math_utils.h"

using namespace Kratos;

int main(int argc, char** argv)
{
    // test the extended knot vector computation
    Vector Ubar_xi;
    int nt_xi;
    std::vector<double> Xi;
    int p = 3;

    Xi.resize(5);
    Xi[0] = 0.0;
    Xi[1] = 1.0;
    Xi[2] = 2.0;
    Xi[3] = 3.0;
    Xi[4] = 4.0;

    IsogeometricMathUtils<double>::compute_extended_knot_vector(Ubar_xi, nt_xi, Xi, p);
    KRATOS_WATCH(Ubar_xi)
    KRATOS_WATCH(nt_xi)

    // test the bezier extraction in 1D
    std::vector<double> Uxi;
    Uxi.resize(4);
    Uxi[0] = 0.5;
    Uxi[1] = 1.5;
    Uxi[2] = 2.5;
    Uxi[3] = 3.5;

    std::vector<int> spans_xi;
    spans_xi.resize(4);
    spans_xi[0] = 1;
    spans_xi[1] = 2;
    spans_xi[2] = 3;
    spans_xi[3] = 4;

    std::vector<Vector> Crows;
    int nb_xi;
    BezierUtils::bezier_extraction_tsplines_1d(Crows, nb_xi, Ubar_xi, Xi, Uxi, spans_xi, p);

    KRATOS_WATCH(nb_xi)
    for (std::size_t i = 0; i < Crows.size(); ++i)
    {
        std::cout << "Crows 1d " << i << ": " << Crows[i] << std::endl;
    }
    KRATOS_WATCH(Ubar_xi)

    // test the bezier extraction in 2D
    std::vector<double> Eta = Xi;
    Vector Ubar_eta;

    std::vector<double> Ueta;
    Ueta.resize(3);
    Ueta[0] = 0.5;
    Ueta[1] = 1.5;
    Ueta[2] = 2.5;

    std::vector<int> spans_eta;
    spans_eta.resize(3);
    spans_eta[0] = 1;
    spans_eta[1] = 2;
    spans_eta[2] = 3;

    int q = 3;
    int nb_eta;

    BezierUtils::bezier_extraction_tsplines_2d(Crows, nb_xi, nb_eta, Ubar_xi, Ubar_eta, Xi, Eta, Uxi, Ueta, spans_xi, spans_eta, p, q);
    KRATOS_WATCH(nb_xi)
    KRATOS_WATCH(nb_eta)
    for (std::size_t i = 0; i < Crows.size(); ++i)
    {
        std::cout << "Crows 2d " << i << ": " << Crows[i] << std::endl;
    }
    KRATOS_WATCH(Ubar_xi)
    KRATOS_WATCH(Ubar_eta)
}
