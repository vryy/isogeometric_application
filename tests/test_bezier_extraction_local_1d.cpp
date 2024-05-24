#include "includes/define.h"
#include "custom_utilities/bezier_utils.h"
#include "custom_utilities/isogeometric_math_utils.h"

using namespace Kratos;

int main(int argc, char** argv)
{
    std::vector<Vector> Crows;
    int nb;
    Vector Ubar;
    std::vector<double> Xi;
    std::vector<double> U;

//    int p = 3;
//    Xi.push_back(0.0);
//    Xi.push_back(0.0);
//    Xi.push_back(1.0);
//    Xi.push_back(1.0);
//    Xi.push_back(2.0);
//    U.push_back(0.1);
//    U.push_back(1.1);
//    U.push_back(0.05);

    int p = 1;
    Xi.push_back(0.0);
    Xi.push_back(0.0);
    Xi.push_back(0.5);

    BezierUtils::bezier_extraction_local_1d(Crows, nb, Ubar, Xi, U, p);

    KRATOS_WATCH(Ubar)
    KRATOS_WATCH(nb)
    for (int i = 0; i < Crows.size(); ++i)
        KRATOS_WATCH(Crows[i])

        return 0;
}
