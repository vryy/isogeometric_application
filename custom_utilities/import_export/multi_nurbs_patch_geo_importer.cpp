//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

// System includes
#include <vector>
#include <fstream>
#include <iomanip>

// External includes

// Project includes
#include "custom_utilities/nurbs/bsplines_fespace.h"
#include "custom_utilities/nurbs/bsplines_patch_utility.h"
#include "custom_utilities/import_export/multi_nurbs_patch_geo_importer.h"

namespace Kratos
{

template<int TDim>
typename Patch<TDim>::Pointer MultiNURBSPatchGeoImporter<TDim>::ImportSingle(const std::string& filename) const
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

    // create the patch
    typename Patch<TDim>::Pointer pNewPatch = this->CreateNewPatch(1, orders, numbers, knots, wcoords, weights);

    // first enumeration
    std::size_t start = 0;
    pNewPatch->pFESpace()->Enumerate(start);

    std::cout << __FUNCTION__ << ": Read NURBS from " << filename << " completed" << std::endl;
    return pNewPatch;
}

template<int TDim>
typename MultiPatch<TDim>::Pointer MultiNURBSPatchGeoImporter<TDim>::Import(const std::string& filename) const
{
    std::ifstream infile(filename.c_str());
    if(!infile)
        KRATOS_THROW_ERROR(std::logic_error, "Error open file", filename)

    std::vector<std::vector<std::size_t> > orders;
    std::vector<std::vector<std::size_t> > numbers;
    std::vector<std::vector<std::vector<double> > > knots;
    std::vector<std::vector<std::vector<double> > > wcoords;
    std::vector<std::vector<double> > weights;
    std::vector<GeoInterface> interfaces;

    // firstly check the version
    std::string firstline;
    std::vector<std::string> words;
    std::getline(infile, firstline);
    boost::trim_if(firstline, boost::is_any_of("\t ")); // ignore trailing spaces
    boost::split(words, firstline, boost::is_any_of(" \t"), boost::token_compress_on);

    if(words[3] == std::string("v.2.1"))
    {
        ReadV21Multi(infile, orders, numbers, knots, wcoords, weights, interfaces);
    }
    else
        KRATOS_THROW_ERROR(std::logic_error, "Unknown NURBS file format", words[3])

    infile.close();

    typename MultiPatch<TDim>::Pointer pNewMultiPatch = typename MultiPatch<TDim>::Pointer(new MultiPatch<TDim>());

    std::size_t npatches = orders.size();

    // create new patch and insert to multipatch
    for (std::size_t i = 0; i < npatches; ++i)
    {
        // create the patch
        typename Patch<TDim>::Pointer pNewPatch = this->CreateNewPatch(i+1,
            orders[i], numbers[i], knots[i], wcoords[i], weights[i]);

        pNewMultiPatch->AddPatch(pNewPatch);
    }

    // incorporate neighbor information
    std::size_t ninterfaces = interfaces.size();
    for (std::size_t i = 0; i < ninterfaces; ++i)
    {
        const GeoInterface& interface = interfaces[i];

        std::size_t patch1_id = static_cast<std::size_t>(interface.patch1);
        BoundarySide side1 = BoundarySideHelper<TDim>::Get(interface.side1);
        std::size_t patch2_id = static_cast<std::size_t>(interface.patch2);
        BoundarySide side2 = BoundarySideHelper<TDim>::Get(interface.side2);

        typename Patch<TDim>::Pointer pPatch1 = pNewMultiPatch->pGetPatch(patch1_id);
        typename Patch<TDim>::Pointer pPatch2 = pNewMultiPatch->pGetPatch(patch2_id);

        if (TDim == 2)
        {
            BoundaryDirection dir;
            if (interface.ornt1 == 1)
                dir = _FORWARD_;
            else if (interface.ornt1 == -1)
                dir = _REVERSED_;

            BSplinesPatchUtility::MakeInterface2D(pPatch1, side1, pPatch2, side2, dir);
        }
        else if (TDim == 3)
        {
            bool uv_or_vu;
            BoundaryDirection dir1, dir2;

            if (interface.flag == 1)
                uv_or_vu = true;
            else
                uv_or_vu = false;

            if (interface.ornt1 == 1)
                dir1 = _FORWARD_;
            else if (interface.ornt1 == -1)
                dir1 = _REVERSED_;

            if (interface.ornt2 == 1)
                dir2 = _FORWARD_;
            else if (interface.ornt2 == -1)
                dir2 = _REVERSED_;

            BSplinesPatchUtility::MakeInterface3D(pPatch1, side1, pPatch2, side2, uv_or_vu, dir1, dir2);
        }
    }

    // first enumeration
    std::size_t start = 0;
    pNewMultiPatch->Enumerate(start);

    std::cout << __FUNCTION__ << ": Read multipatch NURBS from " << filename << " completed" << std::endl;
    return pNewMultiPatch;
}

/// REF: https://github.com/rafavzqz/geopdes/blob/master/geopdes/doc/geo_specs_v06.txt
template<int TDim>
void MultiNURBSPatchGeoImporter<TDim>::ReadV06Single(std::ifstream& infile,
    std::vector<std::size_t>& orders,
    std::vector<std::size_t>& numbers,
    std::vector<std::vector<double> >& knots,
    std::vector<std::vector<double> >& wcoords,
    std::vector<double>& weights) const
{
    // REMARK: i think the v.0.6 is the same as v.0.7, but we have to bear in mind that it's different
    ReadV07Single(infile, orders, numbers, knots, wcoords, weights);
}

/// REF: https://github.com/rafavzqz/geopdes/blob/master/geopdes/doc/geo_specs_v07.txt
template<int TDim>
void MultiNURBSPatchGeoImporter<TDim>::ReadV07Single(std::ifstream& infile,
    std::vector<std::size_t>& orders,
    std::vector<std::size_t>& numbers,
    std::vector<std::vector<double> >& knots,
    std::vector<std::vector<double> >& wcoords,
    std::vector<double>& weights) const
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

/// REF: https://github.com/rafavzqz/geopdes/blob/master/geopdes/doc/geo_specs_v21.txt
template<int TDim>
void MultiNURBSPatchGeoImporter<TDim>::ReadV21Single(std::ifstream& infile,
    std::vector<std::size_t>& orders,
    std::vector<std::size_t>& numbers,
    std::vector<std::vector<double> >& knots,
    std::vector<std::vector<double> >& wcoords,
    std::vector<double>& weights) const
{
    std::string line;
    std::vector<std::string> words;
    int read_mode = _READ_PATCH;
    int ipatch = 0, npatches, rdim;

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

                this->ReadPatchData(infile, rdim, orders, numbers, knots, wcoords, weights);

                break;
            }
        }
    }
}

/// REF: https://github.com/rafavzqz/geopdes/blob/master/geopdes/doc/mp_geo_specs_v21.txt
template<int TDim>
void MultiNURBSPatchGeoImporter<TDim>::ReadV21Multi(std::ifstream& infile,
    std::vector<std::vector<std::size_t> >& orders,
    std::vector<std::vector<std::size_t> >& numbers,
    std::vector<std::vector<std::vector<double> > >& knots,
    std::vector<std::vector<std::vector<double> > >& wcoords,
    std::vector<std::vector<double> >& weights,
    std::vector<GeoInterface>& interfaces) const
{
    std::string line;
    std::vector<std::string> words;
    int read_mode = _READ_PATCH;
    int ipatch = 0, npatches, rdim, number_of_patches_read = 0, number_of_interfaces_read = 0;
    int iinterface = 0, ninterfaces;

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
                ninterfaces = atoi(words[3].c_str());
                KRATOS_WATCH(rdim)
                KRATOS_WATCH(npatches)
                KRATOS_WATCH(ninterfaces)
                orders.resize(npatches);
                numbers.resize(npatches);
                knots.resize(npatches);
                wcoords.resize(npatches);
                weights.resize(npatches);
                interfaces.resize(ninterfaces);
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
                    ipatch = atoi(words[1].c_str()) - 1;
                }
                else
                {
                    KRATOS_THROW_ERROR(std::logic_error, "The patch section has wrong keyword", words[0])
                }

                knots[ipatch].resize(3);
                wcoords[ipatch].resize(3);
                this->ReadPatchData(infile, rdim, orders[ipatch], numbers[ipatch], knots[ipatch], wcoords[ipatch], weights[ipatch]);
                ++number_of_patches_read;

                if (number_of_patches_read < npatches)
                    read_mode = _CHECK_PATCH;
                else
                    read_mode = _READ_INTERFACE;
                continue;
            }

            if(read_mode == _READ_INTERFACE)
            {
                // bound check
                if(words.size() < 2)
                {
                    std::cout << "Error at line: " << line << std::endl;
                    KRATOS_THROW_ERROR(std::logic_error, "The Interface section need to contain INTERFACE and the interface index, current number of information =", words.size())
                }

                if(words[0] == "INTERFACE")
                {
                    iinterface = atoi(words[1].c_str()) - 1;
                }
                else
                {
                    KRATOS_WATCH(words[0])
                    KRATOS_WATCH(words[1])
                    KRATOS_THROW_ERROR(std::logic_error, "The interface section has wrong keyword", words[0])
                }

                std::getline(infile, line);
                boost::trim_if(line, boost::is_any_of("\t ")); // ignore trailing spaces
                boost::split(words, line, boost::is_any_of(" \t"), boost::token_compress_on);

                interfaces[iinterface].patch1 = atoi(words[0].c_str());
                interfaces[iinterface].side1 = atoi(words[1].c_str());

                std::getline(infile, line);
                boost::trim_if(line, boost::is_any_of("\t ")); // ignore trailing spaces
                boost::split(words, line, boost::is_any_of(" \t"), boost::token_compress_on);

                interfaces[iinterface].patch2 = atoi(words[0].c_str());
                interfaces[iinterface].side2 = atoi(words[1].c_str());

                std::getline(infile, line);
                boost::trim_if(line, boost::is_any_of("\t ")); // ignore trailing spaces
                boost::split(words, line, boost::is_any_of(" \t"), boost::token_compress_on);

                if (TDim == 2)
                {
                    interfaces[iinterface].ornt1 = atoi(words[0].c_str());
                }
                else if (TDim == 3)
                {
                    interfaces[iinterface].flag = atoi(words[0].c_str());
                    interfaces[iinterface].ornt1 = atoi(words[1].c_str());
                    interfaces[iinterface].ornt2 = atoi(words[2].c_str());
                }
                ++number_of_interfaces_read;

                if (number_of_interfaces_read == ninterfaces)
                    break;
            }
        }
    }
}

/// Read the data for one patch
template<int TDim>
void MultiNURBSPatchGeoImporter<TDim>::ReadPatchData(std::ifstream& infile,
    const int& rdim,
    std::vector<std::size_t>& orders,
    std::vector<std::size_t>& numbers,
    std::vector<std::vector<double> >& knots,
    std::vector<std::vector<double> >& wcoords,
    std::vector<double>& weights) const
{
    int read_mode = _READ_ORDER, dim_index = 0;
    std::string line;
    std::vector<std::string> words;

    while (read_mode != _NO_READ)
    {
        std::getline(infile, line);
        boost::trim_if(line, boost::is_any_of("\t ")); // ignore trailing spaces
        boost::split(words, line, boost::is_any_of(" \t"), boost::token_compress_on);

        if(words.size() != 0)
        {
            if(words[0] == std::string("#") || words[0][0] == '#')
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

/// Create a new patch based on patch data
template<int TDim>
typename Patch<TDim>::Pointer MultiNURBSPatchGeoImporter<TDim>::CreateNewPatch(const std::size_t& Id,
    const std::vector<std::size_t>& orders,
    const std::vector<std::size_t>& numbers,
    const std::vector<std::vector<double> >& knots,
    const std::vector<std::vector<double> >& wcoords,
    const std::vector<double>& weights) const
{
    // create the FESpace
    typename BSplinesFESpace<TDim>::Pointer pNewFESpace = BSplinesFESpace<TDim>::Create();
    for (int dim = 0; dim < TDim; ++dim)
    {
        pNewFESpace->SetKnotVector(dim, knots[dim]);
        pNewFESpace->SetInfo(dim, numbers[dim], orders[dim]);
    }

    // reset function indices and enumerate it first time to give each function in the FESpace a different id
    pNewFESpace->ResetFunctionIndices();

    // create new patch
    typename Patch<TDim>::Pointer pNewPatch = Patch<TDim>::Create(Id, pNewFESpace);

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

    return pNewPatch;
}

// template instantiation
template class MultiNURBSPatchGeoImporter<1>;
template class MultiNURBSPatchGeoImporter<2>;
template class MultiNURBSPatchGeoImporter<3>;

} // namespace Kratos.

