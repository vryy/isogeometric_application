//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 20 Aug 2015 $
//   Revision:            $Revision: 1.0 $
//
//

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "bezier_utils.h"
#include "bspline_utils.h"

namespace Kratos
{

void BSplineUtils::test_ComputeBsplinesKnotInsertionCoefficients1DLocal()
{
    std::vector<double> local_knots;
    std::vector<double> ins_knots;

    local_knots.push_back(0);
    local_knots.push_back(0);
    local_knots.push_back(1);
    local_knots.push_back(2);
    local_knots.push_back(3);

    ins_knots.push_back(0.5);
    ins_knots.push_back(1.5);
    ins_knots.push_back(2.5);

    int p = 3;
    Vector D;
    std::vector<double> new_knots;
    ComputeBsplinesKnotInsertionCoefficients1DLocal(D, new_knots, p, local_knots, ins_knots);

    std::cout << "new_knots:";
    for (unsigned int i = 0; i < new_knots.size(); ++i)
    {
        std::cout << " " << new_knots[i];
    }
    std::cout << std::endl;

    KRATOS_WATCH(D)
}

}// namespace Kratos.
