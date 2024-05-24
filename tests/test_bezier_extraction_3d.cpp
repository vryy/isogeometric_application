#include "includes/define.h"
#include "custom_utilities/bezier_utils.h"
#include "custom_utilities/isogeometric_math_utils.h"

using namespace Kratos;

int main(int argc, char** argv)
{
    std::vector<Matrix> C;
    int nb1, nb2, nb3;
    std::vector<double> U = {0,     0,    0,  1.0000,  1.0000,  1.0000};
    std::vector<double> V = {0,     0,    0,  0.5,  1.0000,  1.0000,  1.0000};
    std::vector<double> W = {0,     0,  1.0000,  1.0000};
    int p = 2;
    int q = 2;
    int r = 1;

    BezierUtils::bezier_extraction_3d(C, nb1, nb2, nb3, U, V, W, p, q, r);

    KRATOS_WATCH(nb1)
    KRATOS_WATCH(nb2)
    KRATOS_WATCH(nb3)
    KRATOS_WATCH(C.size())
    for (int i = 0; i < C.size(); ++i)
        KRATOS_WATCH(C[i])

        return 0;
}
