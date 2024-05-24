#include "includes/define.h"
#include "custom_utilities/bezier_utils.h"
#include "custom_utilities/isogeometric_math_utils.h"

using namespace Kratos;

int main(int argc, char** argv)
{
    std::vector<Matrix> C;
    int nb1, nb2;
    std::vector<double> U = {0,     0,    0,  1.0 / 3,  2.0 / 3,  1.0000,  1.0000,  1.0000};
    std::vector<double> V = {0,     0,    0,  1.0 / 3,  2.0 / 3,  1.0000,  1.0000,  1.0000};
    int p = 2;
    int q = 2;

    BezierUtils::bezier_extraction_2d(C, nb1, nb2, U, V, p, q);

    KRATOS_WATCH(nb1)
    KRATOS_WATCH(nb2)
    KRATOS_WATCH(C.size())
    for (int i = 0; i < C.size(); ++i)
        KRATOS_WATCH(C[i])

        return 0;
}
