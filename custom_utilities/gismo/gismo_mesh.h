//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 15 Feb 2016 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_GISMO_MESH_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_GISMO_MESH_H_INCLUDED

// System includes
#include <ctime>
#include <cmath>
#include <string>
#include <vector>
#include <deque>
#include <iostream>
#include <fstream>
#include <set>
#include <list>

// External includes
#include <omp.h>
#include "boost/progress.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/foreach.hpp"
#include "boost/python.hpp"
#include "boost/python/stl_iterator.hpp"
#include "gismo.h"

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "custom_utilities/iga_define.h"

namespace Kratos
{

/**
    Class represents a Hierarchical NURBS mesh by G+Smo
*/
class GismoMesh : public IsogeometricEcho
{
public:
    /// Const definition
    static const int NO_READ          = 0;
    static const int READ_PATCH       = 1;
    static const int READ_ORDER       = 2;
    static const int READ_NUMBER      = 3;
    static const int READ_KNOTS       = 4;
    static const int READ_COORDINATES = 5;
    static const int READ_WEIGHTS     = 6;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(GismoMesh);

    /// Type definition

    /// Default constructor
    GismoMesh(std::string Name) : mName(Name)
    {}

    /// Destructor
    ~GismoMesh()
    {}

    /**************************************************************************
                           MODIFICATION SUBROUTINES
    **************************************************************************/

    /// Get the name of this hierarchical mesh
    std::string Name() const {return mName;}

    /**************************************************************************
                            INPUT SUBROUTINES
    **************************************************************************/

    /// Read and construct the first level from file
    void ReadMesh(std::string fn)
    {
        std::ifstream infile(fn.c_str());
        if (!infile)
            KRATOS_ERROR << "Error open file " << fn;

        std::string line;
        std::vector<std::string> words;
        int read_mode = READ_PATCH;
        int npatches, dim_index = 0;
        std::vector<std::vector<int> > orders;
        std::vector<std::vector<int> > numbers;
        std::vector<std::vector<double> > knots_1;
        std::vector<std::vector<double> > knots_2;
        std::vector<std::vector<double> > knots_3;
        std::vector<std::vector<double> > x_coords;
        std::vector<std::vector<double> > y_coords;
        std::vector<std::vector<double> > z_coords;
        std::vector<std::vector<double> > weights;
        while (!infile.eof())
        {
            std::getline(infile, line);
            boost::trim_if(line, boost::is_any_of("\t ")); // ignore trailing spaces
            boost::split(words, line, boost::is_any_of(" \t"), boost::token_compress_on);

            if (words.size() != 0)
            {
                if (words[0] == std::string("#") || words[0][0] == '#')
                {
                    continue;
                }

                if (read_mode == READ_PATCH)
                {
                    // bound check
                    if (words.size() < 2)
                    {
                        KRATOS_ERROR << "The Patch section need to contain information about dimension and number of patches, current number of information = " << words.size();
                    }

                    // read info
                    mDim = atoi(words[0].c_str());
                    npatches = atoi(words[1].c_str());
                    read_mode = READ_ORDER;
                    continue;
                }

                if (read_mode == READ_ORDER)
                {
                    // bound check
                    if (words.size() != mDim)
                        KRATOS_ERROR << "The Order section must contained number of information equal to dimension, current number of information = " << words.size();

                    // read info
                    std::vector<int> dummy;
                    orders.push_back(dummy);
                    for (std::size_t i = 0; i < mDim; ++i)
                    {
                        orders.back().push_back(atoi(words[i].c_str()));
                    }
                    read_mode = READ_NUMBER;
                    continue;
                }

                if (read_mode == READ_NUMBER)
                {
                    // bound check
                    if (words.size() != mDim)
                        KRATOS_ERROR << "The Number section must contained number of information equal to dimension, current number of information = " << words.size();

                    // read info
                    std::vector<int> dummy;
                    numbers.push_back(dummy);
                    for (std::size_t i = 0; i < mDim; ++i)
                    {
                        numbers.back().push_back(atoi(words[i].c_str()));
                    }
                    read_mode = READ_KNOTS;
                    continue;
                }

                if (read_mode == READ_KNOTS)
                {
                    // bound check
                    int knot_len = numbers.back()[dim_index] + orders.back()[dim_index] + 1;
                    if (words.size() != knot_len)
                        KRATOS_ERROR << "The Knots section must contained number of information equal to n+p+1, current number of information = " << words.size();

                    // read info
                    if (dim_index == 0)
                    {
                        std::vector<double> dummy;
                        knots_1.push_back(dummy);
                    }
                    else if (dim_index == 1)
                    {
                        std::vector<double> dummy;
                        knots_2.push_back(dummy);
                    }
                    else if (dim_index == 2)
                    {
                        std::vector<double> dummy;
                        knots_3.push_back(dummy);
                    }

                    for (std::size_t i = 0; i < knot_len; ++i)
                    {
                        double k = atof(words[i].c_str());
                        if (dim_index == 0)
                        {
                            knots_1.back().push_back(k);
                        }
                        else if (dim_index == 1)
                        {
                            knots_2.back().push_back(k);
                        }
                        else if (dim_index == 2)
                        {
                            knots_3.back().push_back(k);
                        }
                        else
                            KRATOS_ERROR << "Wrong knot dimension index. Something must be wrong";
                    }

                    ++dim_index;
                    if (dim_index == mDim)
                    {
                        dim_index = 0;
                        read_mode = READ_COORDINATES;
                    }
                    continue;
                }

                if (read_mode == READ_COORDINATES)
                {
                    // bound check
                    int num_basis = 1;
                    for (std::size_t i = 0; i < mDim; ++i)
                    {
                        num_basis *= numbers.back()[i];
                    }
                    if (words.size() != num_basis)
                        KRATOS_ERROR << "The Coordinates section must contained number of information equal to prod(ni), current number of information = ", words.size();

                    if (dim_index == 0)
                    {
                        std::vector<double> dummy;
                        x_coords.push_back(dummy);
                        for (std::size_t i = 0; i < num_basis; ++i)
                        {
                            x_coords.back().push_back(atof(words[i].c_str()));
                        }
                    }
                    else if (dim_index == 1)
                    {
                        std::vector<double> dummy;
                        y_coords.push_back(dummy);
                        for (std::size_t i = 0; i < num_basis; ++i)
                        {
                            y_coords.back().push_back(atof(words[i].c_str()));
                        }
                    }
                    else if (dim_index == 2)
                    {
                        std::vector<double> dummy;
                        z_coords.push_back(dummy);
                        for (std::size_t i = 0; i < num_basis; ++i)
                        {
                            z_coords.back().push_back(atof(words[i].c_str()));
                        }
                    }

                    ++dim_index;
                    if (dim_index == mDim)
                    {
                        dim_index = 0;
                        read_mode = READ_WEIGHTS;
                    }
                    continue;
                }

                if (read_mode == READ_WEIGHTS)
                {
                    // bound check
                    int num_basis = 1;
                    for (std::size_t i = 0; i < mDim; ++i)
                    {
                        num_basis *= numbers.back()[i];
                    }
                    if (words.size() != num_basis)
                        KRATOS_ERROR << "The Weights section must contained number of information equal to prod(ni), current number of information = " << words.size();

                    // read info
                    std::vector<double> dummy;
                    weights.push_back(dummy);
                    for (std::size_t i = 0; i < num_basis; ++i)
                    {
                        weights.back().push_back(atof(words[i].c_str()));
                    }

                    if (npatches > 1)
                    {
                        read_mode = READ_ORDER;
                        npatches--;
                    }
                    else
                    {
                        read_mode = NO_READ;
                    }
                    continue;
                }
            }
        }

        // close the file
        infile.close();

        // create Gismo geometries
        npatches = orders.size();
        if (mDim == 2)
        {
            for (unsigned int ipatch = 0; ipatch < npatches; ++ipatch)
            {
                gismo::gsKnotVector<real_t> kv1, kv2;
                FillKnotVector(knots_1[ipatch], kv1, 1.0e-6);
                FillKnotVector(knots_2[ipatch], kv2, 1.0e-6);
                KRATOS_WATCH(kv1)
                KRATOS_WATCH(kv2)

                unsigned int num_basis = numbers[ipatch][0] * numbers[ipatch][1];
                KRATOS_WATCH(num_basis)

                gismo::gsMatrix<real_t> coefs(num_basis, 3);
                for (unsigned int i = 0; i < num_basis; ++i)
                {
                    coefs << x_coords[ipatch][i] / weights[ipatch][i], y_coords[ipatch][i] / weights[ipatch][i], weights[ipatch][i];
                }

                gismo::gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);

                gismo::gsTensorBSpline<2, real_t> patch(basis, coefs);
                KRATOS_WATCH(patch)
                gismo::gsWriteParaview(patch, mName.c_str());
            }
        }
        else if (mDim == 3)
        {
            for (unsigned int ipatch = 0; ipatch < npatches; ++ipatch)
            {
                gismo::gsKnotVector<real_t> kv1, kv2, kv3;
                FillKnotVector(knots_1[ipatch], kv1, 1.0e-6);
                FillKnotVector(knots_2[ipatch], kv2, 1.0e-6);
                FillKnotVector(knots_3[ipatch], kv3, 1.0e-6);
                KRATOS_WATCH(kv1)
                KRATOS_WATCH(kv2)
                KRATOS_WATCH(kv3)

                unsigned int num_basis = numbers[ipatch][0] * numbers[ipatch][1] * numbers[ipatch][2];

                gismo::gsMatrix<real_t> coefs(num_basis, 4);
                for (unsigned int i = 0; i < num_basis; ++i)
                {
                    coefs << x_coords[ipatch][i] / weights[ipatch][i], y_coords[ipatch][i] / weights[ipatch][i], z_coords[ipatch][i] / weights[ipatch][i], weights[ipatch][i];
                }

                gismo::gsTensorBSplineBasis<3, real_t> basis(kv1, kv2, kv3);

                gismo::gsTensorBSpline<3, real_t> patch(basis, coefs);
                KRATOS_WATCH(patch)
            }
        }

        if (GetEchoLevel() > 0)
        {
            std::cout << __FUNCTION__ << ": Traverse file completed" << std::endl;
        }
    }

    /**************************************************************************
                            REFINEMENT SUBROUTINES
    **************************************************************************/

    /**************************************************************************
                            INFORMATION SUBROUTINES
    **************************************************************************/

    /// Print information of this hierarchical mesh
    void PrintInfo(std::ostream& rOStream) const
    {}

    /// Print data of this hierarchical mesh
    void PrintData(std::ostream& rOStream) const
    {}

    /**************************************************************************
                            MATLAB INTERFACE
    **************************************************************************/

    /**************************************************************************
                            KRATOS INTERFACE
    **************************************************************************/

    /**************************************************************************
                            DEBUG INTERFACE
    **************************************************************************/

private:
    std::string mName;
    unsigned int mDim, mEchoLevel;

    template<class knot_vector_t>
    void FillKnotVector(const std::vector<double>& knots, knot_vector_t& knotvec, double tol)
    {
        std::cout << "knots:";
        for (unsigned int i = 0; i < knots.size(); ++i)
        {
            std::cout << " " << knots[i];
        }
        std::cout << std::endl;
        double k;
        int cnt = 0, mult;
        k = knots[cnt];
        mult = 1;
        while (cnt < knots.size() - 1)
        {
            ++cnt;
            if (fabs(knots[cnt] - k) < tol)
            {
                ++mult;
            }
            else
            {
                knotvec.insert(k, mult);
                k = knots[cnt];
                mult = 1;
            }
        }
        knotvec.insert(k, mult);
    }

};

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const GismoMesh& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_GISMO_MESH_H_INCLUDED
