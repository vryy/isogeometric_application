#include "includes/define.h"
#include "custom_utilities/nurbs/knot_array_1d.h"

using namespace Kratos;

int main(int argc, char** argv)
{
    typedef KnotArray1D<double> knot_container_t;
    typedef knot_container_t::knot_t knot_t;

    knot_container_t knot_vector;

    knot_vector.pCreateKnot(0.0);
    knot_vector.pCreateKnot(0.0);
    knot_vector.pCreateKnot(0.0);
    knot_vector.pCreateKnot(0.0);
    knot_vector.pCreateKnot(0.25);
    knot_vector.pCreateKnot(0.5);
    knot_vector.pCreateKnot(0.5);
    knot_vector.pCreateKnot(0.75);
    knot_vector.pCreateKnot(1.0);
    knot_vector.pCreateKnot(1.0);
    knot_vector.pCreateKnot(1.0);
    knot_vector.pCreateKnot(1.0);

    KRATOS_WATCH(knot_vector)
    KRATOS_WATCH(knot_vector.nspans())

    for (std::size_t i = 0; i < knot_vector.nspans(); ++i)
    {
        std::cout << "span " << i + 1 << ":" << std::endl;
        std::tuple<knot_t, knot_t> span = knot_vector.span(i + 1);
        KRATOS_WATCH(std::get<0>(span)->Index())
        KRATOS_WATCH(std::get<0>(span)->Value())
        KRATOS_WATCH(std::get<1>(span)->Index())
        KRATOS_WATCH(std::get<1>(span)->Value())
        std::cout << "---------------" << std::endl;
    }

    return 0;
}
