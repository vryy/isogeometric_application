#include "includes/define.h"
#include "custom_utilities/bezier_utils.h"
#include "custom_utilities/isogeometric_math_utils.h"

using namespace Kratos;

int main(int argc, char** argv)
{
    std::vector<Matrix> C;
    int nb;
    std::vector<double> U = {0.00, 0.00, 0.00, 0.00, 0.25, 0.50, 0.75, 1.00, 1.00, 1.00, 1.00};
    int p = 3;

    BezierUtils::bezier_extraction_1d(C, nb, U, p);

    KRATOS_WATCH(nb)
    KRATOS_WATCH(C.size())
    for (int i = 0; i < C.size(); ++i)
        KRATOS_WATCH(C[i])

        return 0;
}
