#include "includes/define.h"
#include "custom_utilities/bspline_utils.h"

using namespace Kratos;

void LocalBasisFun(int p, double xi, const std::vector<double>& Xi)
{
//    int span = BSplineUtils::FindSpan2(1, p, xi, Xi) + 1;
//    std::vector<double> v(1);
//    BSplineUtils::BasisFuns(v, span, xi, p, Xi);
//    std::cout << "local knots:";
//    for (std::size_t i = 0; i < Xi.size(); ++i)
//        std::cout << " " << Xi[i];
//    std::cout << "; p = " << p;
//    std::cout << "; xi = " << xi;
//    std::cout << "; span = " << span;
//    std::cout << "; v = " << v[0];
//    std::cout << std::endl;

    double v = BSplineUtils::CoxDeBoor(xi, 0, p, Xi);
    std::cout << "local knots:";
    for (std::size_t i = 0; i < Xi.size(); ++i)
        std::cout << " " << Xi[i];
    std::cout << "; p = " << p;
    std::cout << "; xi = " << xi;
    std::cout << "; v = " << v;
    std::cout << std::endl;
}


int main(int argc, char** argv)
{
    std::vector<double> Xi1 = {0.0, 0.0, 0.5};
    LocalBasisFun(1, 0.0, Xi1);
    LocalBasisFun(1, 0.2, Xi1);
    LocalBasisFun(1, 0.25, Xi1);
    LocalBasisFun(1, 0.5, Xi1);

    std::cout << "-----------------" << std::endl;

    std::vector<double> Xi2 = {0.0, 0.5, 1.0};
    LocalBasisFun(1, 0.0, Xi2);
    LocalBasisFun(1, 0.1, Xi2);
    LocalBasisFun(1, 0.2, Xi2);
    LocalBasisFun(1, 0.3, Xi2);
    LocalBasisFun(1, 0.4, Xi2);
    LocalBasisFun(1, 0.5, Xi2);
    LocalBasisFun(1, 0.6, Xi2);
    LocalBasisFun(1, 0.7, Xi2);
    LocalBasisFun(1, 0.8, Xi2);
    LocalBasisFun(1, 0.9, Xi2);
    LocalBasisFun(1, 1.0, Xi2);

    std::cout << "-----------------" << std::endl;

    std::vector<double> Xi3 = {0.5, 1.0, 1.0};
    LocalBasisFun(1, 0.5, Xi3);
    LocalBasisFun(1, 0.75, Xi3);
    LocalBasisFun(1, 0.8, Xi3);
    LocalBasisFun(1, 0.99, Xi3);
    LocalBasisFun(1, 1.0, Xi3);
}

