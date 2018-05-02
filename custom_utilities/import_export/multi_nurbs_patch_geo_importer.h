//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_NURBS_PATCH_GEO_IMPORTER_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_NURBS_PATCH_GEO_IMPORTER_H_INCLUDED

// System includes
#include <vector>
#include <fstream>
#include <iomanip>

// External includes
#include <boost/algorithm/string.hpp>

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"
#include "custom_utilities/import_export/multipatch_importer.h"

namespace Kratos
{

enum ReadMode
{
    _NO_READ          = 0,
    _READ_PATCH       = 1,
    _CHECK_PATCH      = 7,
    _READ_ORDER       = 2,
    _READ_NUMBER      = 3,
    _READ_KNOTS       = 4,
    _READ_COORDINATES = 5,
    _READ_WEIGHTS     = 6
};


/// Get the dimension of underlying NURBS in geo file
static int GetDimensionOfGeoHelper(const std::string& fn)
{
    std::ifstream infile(fn.c_str());
    if(!infile)
        KRATOS_THROW_ERROR(std::logic_error, "Error open file", fn)

    std::string line;
    std::vector<std::string> words;
    int read_mode = _READ_PATCH;
    while(!infile.eof())
    {
        std::getline(infile, line);
        boost::trim_if(line, boost::is_any_of("\t ")); // ignore trailing spaces
        boost::split(words, line, boost::is_any_of(" \t"), boost::token_compress_on);

        if(words.size() != 0)
        {
            if(words[0] == std::string("#") || words[0][0] == '#')
                continue;

            if(read_mode == _READ_PATCH)
            {
                // bound check
                if(words.size() < 2)
                {
                    std::cout << "Error at line: " << line << std::endl;
                    KRATOS_THROW_ERROR(std::logic_error, "The Patch section need to contain information about dimension and number of patches, current number of information =", words.size())
                }

                // read info
                int Dim = atoi(words[0].c_str());
                return Dim;
            }
        }
    }

    infile.close();

    return 0;
}


template<int TDim>
class MultiNURBSPatchGeoImporter : public MultiPatchImporter<TDim>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MultiNURBSPatchGeoImporter);

    virtual typename Patch<TDim>::Pointer ImportSingle(const std::string& filename)
    {
        std::ifstream infile(filename.c_str());
        if(!infile)
            KRATOS_THROW_ERROR(std::logic_error, "Error open file", filename)

        std::vector<std::size_t> orders;
        std::vector<std::size_t> numbers;
        std::vector<std::vector<double> > knots(3);
        std::vector<std::vector<double> > wcoords(3);
        std::vector<double> weights;

        // firstly check the version
        std::string firstline;
        std::vector<std::string> words;
        std::getline(infile, firstline);
        boost::trim_if(firstline, boost::is_any_of("\t ")); // ignore trailing spaces
        boost::split(words, firstline, boost::is_any_of(" \t"), boost::token_compress_on);

        if(words[3] == std::string("v.0.6"))
        {
            ReadV06Single(infile, orders, numbers, knots, wcoords, weights);
        }
        else if(words[3] == std::string("v.0.7"))
        {
            ReadV07Single(infile, orders, numbers, knots, wcoords, weights);
        }
        else if(words[3] == std::string("v.2.1"))
        {
            ReadV21Single(infile, orders, numbers, knots, wcoords, weights);
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Unknown NURBS file format", words[3])

        infile.close();

        // create the FESpace
        typename BSplinesFESpace<TDim>::Pointer pNewFESpace = BSplinesFESpace<TDim>::Create();
        for (int dim = 0; dim < TDim; ++dim)
        {
            pNewFESpace->SetKnotVector(dim, knots[dim]);
            pNewFESpace->SetInfo(dim, numbers[dim], orders[dim]);
        }

        // reset function indices and enumerate it first time to give each function in the FESpace a different id
        pNewFESpace->ResetFunctionIndices();
        std::size_t start = 0;
        pNewFESpace->Enumerate(start);

        // create new patch
        typename Patch<TDim>::Pointer pNewPatch = Patch<TDim>::Create(0, pNewFESpace);

        // create control grid and assign to new patch
        typedef ControlPoint<double> ControlPointType;
        typename StructuredControlGrid<TDim, ControlPointType>::Pointer pControlPointGrid = StructuredControlGrid<TDim, ControlPointType>::Create(numbers);
        std::size_t total_number = 1;
        for (int dim = 0; dim < TDim; ++dim)
            total_number *= numbers[dim];

        for (std::size_t i = 0; i < total_number; ++i)
        {
            ControlPointType c;
            if (TDim == 2)
                c.SetCoordinates(wcoords[0][i]/weights[i], wcoords[1][i]/weights[i], 0.0, weights[i]);
            else if (TDim == 3)
                c.SetCoordinates(wcoords[0][i]/weights[i], wcoords[1][i]/weights[i], wcoords[2][i]/weights[i], weights[i]);
            pControlPointGrid->SetData(i, c);
        }

        pControlPointGrid->SetName("CONTROL_POINT");
        pNewPatch->CreateControlPointGridFunction(pControlPointGrid);

        std::cout << __FUNCTION__ << ": Read NURBS from " << filename << " completed" << std::endl;
        return pNewPatch;
    }

    virtual typename MultiPatch<TDim>::Pointer Import(const std::string& filename)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not implemented")
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiNURBSPatchGeoImporter";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:

    void ReadV06Single(std::ifstream& infile, std::vector<std::size_t>& orders,
        std::vector<std::size_t>& numbers,
        std::vector<std::vector<double> >& knots,
        std::vector<std::vector<double> >& wcoords,
        std::vector<double>& weights)
    {
        // REMARK: i think the v.0.6 is the same as v.0.7, but we have to bear in mind that it's different
        ReadV07Single(infile, orders, numbers, knots, wcoords, weights);
    }

    void ReadV07Single(std::ifstream& infile, std::vector<std::size_t>& orders,
        std::vector<std::size_t>& numbers,
        std::vector<std::vector<double> >& knots,
        std::vector<std::vector<double> >& wcoords,
        std::vector<double>& weights)
    {
        std::string line;
        std::vector<std::string> words;
        int read_mode = _READ_PATCH;
        int npatches, dim_index = 0;

        while(!infile.eof())
        {
            std::getline(infile, line);
            boost::trim_if(line, boost::is_any_of("\t ")); // ignore trailing spaces
            boost::split(words, line, boost::is_any_of(" \t"), boost::token_compress_on);

            if(words.size() != 0)
            {
                if(words[0] == std::string("#") || words[0][0] == '#')
                    continue;

                if(read_mode == _READ_PATCH)
                {
                    // bound check
                    if(words.size() < 2)
                    {
                        std::cout << "Error at line: " << line << std::endl;
                        KRATOS_THROW_ERROR(std::logic_error, "The Patch section need to contain information about dimension and number of patches, current number of information =", words.size())
                    }

                    // read info
                    int Dim = atoi(words[0].c_str());
                    if (Dim != TDim)
                        KRATOS_THROW_ERROR(std::logic_error, "The input dimension is invalid", "")
                    npatches = atoi(words[1].c_str());
                    if(npatches > 1)
                    {
                        KRATOS_WATCH(line)
                        KRATOS_WATCH(words[0])
                        KRATOS_WATCH(words[1])
                        KRATOS_THROW_ERROR(std::logic_error, "At present, the number of patches > 1 is not supported, npatches =", npatches)
                    }
                    read_mode = _READ_ORDER;
                    continue;
                }

                if(read_mode == _READ_ORDER)
                {
                    // bound check
                    if(words.size() != TDim)
                        KRATOS_THROW_ERROR(std::logic_error, "The Order section must contained number of information equal to dimension, current number of information =", words.size())

                    // read info
                    for(std::size_t i = 0; i < TDim; ++i)
                        orders.push_back(static_cast<std::size_t>(atoi(words[i].c_str())));
                    read_mode = _READ_NUMBER;
                    continue;
                }

                if(read_mode == _READ_NUMBER)
                {
                    // bound check
                    if(words.size() != TDim)
                        KRATOS_THROW_ERROR(std::logic_error, "The Number section must contained number of information equal to dimension, current number of information =", words.size())

                    for(std::size_t i = 0; i < TDim; ++i)
                        numbers.push_back(static_cast<std::size_t>(atoi(words[i].c_str())));
                    read_mode = _READ_KNOTS;
                    continue;
                }

                if(read_mode == _READ_KNOTS)
                {
                    // bound check
                    int knot_len = numbers[dim_index] + orders[dim_index] + 1;
                    if(words.size() != knot_len)
                        KRATOS_THROW_ERROR(std::logic_error, "The Knots section must contained number of information equal to n+p+1, current number of information =", words.size())

                    for(std::size_t i = 0; i < knot_len; ++i)
                    {
                        double k = atof(words[i].c_str());
                        knots[dim_index].push_back(k);
                    }

                    ++dim_index;
                    if(dim_index == TDim)
                    {
                        dim_index = 0;
                        read_mode = _READ_COORDINATES;
                    }
                    continue;
                }

                if(read_mode == _READ_COORDINATES)
                {
                    // bound check
                    int num_basis = 1;
                    for(std::size_t i = 0; i < TDim; ++i)
                        num_basis *= numbers[i];
                    if(words.size() != num_basis)
                        KRATOS_THROW_ERROR(std::logic_error, "The Coordinates section must contained number of information equal to prod(ni), current number of information =", words.size())

                    for(std::size_t i = 0; i < num_basis; ++i)
                        wcoords[dim_index].push_back(atof(words[i].c_str()));

                    ++dim_index;
                    if(dim_index == TDim)
                    {
                        dim_index = 0;
                        read_mode = _READ_WEIGHTS;
                    }
                    continue;
                }

                if(read_mode == _READ_WEIGHTS)
                {
                    // bound check
                    int num_basis = 1;
                    for(std::size_t i = 0; i < TDim; ++i)
                        num_basis *= numbers[i];
                    if(words.size() != num_basis)
                        KRATOS_THROW_ERROR(std::logic_error, "The Weights section must contained number of information equal to prod(ni), current number of information =", words.size())

                    for(std::size_t i = 0; i < num_basis; ++i)
                        weights.push_back(atof(words[i].c_str()));

                    read_mode = _NO_READ;
                    continue;
                }
            }
        }
    }

    void ReadV21Single(std::ifstream& infile, std::vector<std::size_t>& orders,
        std::vector<std::size_t>& numbers,
        std::vector<std::vector<double> >& knots,
        std::vector<std::vector<double> >& wcoords,
        std::vector<double>& weights)
    {
        std::string line;
        std::vector<std::string> words;
        int read_mode = _READ_PATCH;
        int ipatch = 0, npatches, rdim, dim_index = 0;

        while(!infile.eof())
        {
            std::getline(infile, line);
            boost::trim_if(line, boost::is_any_of("\t ")); // ignore trailing spaces
            boost::split(words, line, boost::is_any_of(" \t"), boost::token_compress_on);

            if(words.size() != 0)
            {
                if(words[0] == std::string("#") || words[0][0] == '#')
                    continue;

                if(read_mode == _READ_PATCH)
                {
                    // bound check
                    if(words.size() < 2)
                    {
                        std::cout << "Error at line: " << line << std::endl;
                        KRATOS_THROW_ERROR(std::logic_error, "The Patch section need to contain information about dimension and number of patches, current number of information =", words.size())
                    }

                    // read info
                    int Dim = atoi(words[0].c_str());
                    if (Dim != TDim)
                        KRATOS_THROW_ERROR(std::logic_error, "The input dimension is invalid", "")
                    rdim = atoi(words[1].c_str());
                    npatches = atoi(words[2].c_str());
                    KRATOS_WATCH(rdim)
                    KRATOS_WATCH(npatches)
                    if(npatches > 1)
                    {
                        KRATOS_WATCH(line)
                        KRATOS_WATCH(words[0])
                        KRATOS_WATCH(words[1])
                        KRATOS_THROW_ERROR(std::logic_error, "At present, the number of patches > 1 is not supported, npatches =", npatches)
                    }
                    read_mode = _CHECK_PATCH;
                    continue;
                }

                if(read_mode == _CHECK_PATCH)
                {
                    // bound check
                    if(words.size() < 2)
                    {
                        std::cout << "Error at line: " << line << std::endl;
                        KRATOS_THROW_ERROR(std::logic_error, "The Patch section need to contain PATCH and the patch index, current number of information =", words.size())
                    }

                    if(words[0] == "PATCH")
                    {
                        if(ipatch == npatches)
                            break;
                        ++ipatch;
                    }
                    else
                        KRATOS_THROW_ERROR(std::logic_error, "The patch section has wrong keyword", words[0])
                    read_mode = _READ_ORDER;
                    continue;
                }

                if(read_mode == _READ_ORDER)
                {
                    // bound check
                    if(words.size() != TDim)
                        KRATOS_THROW_ERROR(std::logic_error, "The Order section must contained number of information equal to dimension, current number of information =", words.size())

                    // read info
                    for(std::size_t i = 0; i < TDim; ++i)
                        orders.push_back(static_cast<std::size_t>(atoi(words[i].c_str())));
                    read_mode = _READ_NUMBER;
                    continue;
                }

                if(read_mode == _READ_NUMBER)
                {
                    // bound check
                    if(words.size() != TDim)
                        KRATOS_THROW_ERROR(std::logic_error, "The Number section must contained number of information equal to dimension, current number of information =", words.size())

                    for(std::size_t i = 0; i < TDim; ++i)
                        numbers.push_back(static_cast<std::size_t>(atoi(words[i].c_str())));
                    read_mode = _READ_KNOTS;
                    continue;
                }

                if(read_mode == _READ_KNOTS)
                {
                    // bound check
                    int knot_len = numbers[dim_index] + orders[dim_index] + 1;
                    if(words.size() != knot_len)
                        KRATOS_THROW_ERROR(std::logic_error, "The Knots section must contained number of information equal to n+p+1, current number of information =", words.size())

                    for(std::size_t i = 0; i < knot_len; ++i)
                    {
                        double k = atof(words[i].c_str());
                        knots[dim_index].push_back(k);
                    }

                    ++dim_index;
                    if(dim_index == TDim)
                    {
                        dim_index = 0;
                        read_mode = _READ_COORDINATES;
                    }
                    continue;
                }

                if(read_mode == _READ_COORDINATES)
                {
                    // bound check
                    int num_basis = 1;
                    for(std::size_t i = 0; i < TDim; ++i)
                        num_basis *= numbers[i];
                    if(words.size() != num_basis)
                        KRATOS_THROW_ERROR(std::logic_error, "The Coordinates section must contained number of information equal to prod(ni), current number of information =", words.size())

                    for(std::size_t i = 0; i < num_basis; ++i)
                        wcoords[dim_index].push_back(atof(words[i].c_str()));

                    ++dim_index;
                    if(dim_index == rdim)
                    {
                        dim_index = 0;
                        read_mode = _READ_WEIGHTS;
                    }
                    continue;
                }

                if(read_mode == _READ_WEIGHTS)
                {
                    // bound check
                    int num_basis = 1;
                    for(std::size_t i = 0; i < TDim; ++i)
                        num_basis *= numbers[i];
                    if(words.size() != num_basis)
                        KRATOS_THROW_ERROR(std::logic_error, "The Weights section must contained number of information equal to prod(ni), current number of information =", words.size())

                    for(std::size_t i = 0; i < num_basis; ++i)
                        weights.push_back(atof(words[i].c_str()));

                    read_mode = _NO_READ;
                    continue;
                }
            }
        }
    }

};


/// output stream function
template<int TDim>
inline std::ostream& operator <<(std::ostream& rOStream, const MultiNURBSPatchGeoImporter<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}


} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_NURBS_PATCH_GEO_IMPORTER_H_INCLUDED defined

