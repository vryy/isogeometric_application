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
#include "custom_utilities/import_export/multipatch_importer.h"

namespace Kratos
{

enum ReadMode
{
    _NO_READ          = 0,
    _READ_PATCH       = 1,
    _READ_ORDER       = 2,
    _READ_NUMBER      = 3,
    _READ_KNOTS       = 4,
    _READ_COORDINATES = 5,
    _READ_WEIGHTS     = 6,
    _CHECK_PATCH      = 7,
    _READ_INTERFACE   = 8,
    _READ_BOUNDARY    = 9
};

struct GeoInterface
{
    int patch1;
    int side1;
    int patch2;
    int side2;
    int flag;
    int ornt1;
    int ornt2;
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
struct BoundarySideHelper
{
    static BoundarySide Get(const int& i)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not implemented")
    }
};

template<int TDim>
class MultiNURBSPatchGeoImporter : public MultiPatchImporter<TDim>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MultiNURBSPatchGeoImporter);

    virtual typename Patch<TDim>::Pointer ImportSingle(const std::string& filename) const;

    virtual typename MultiPatch<TDim>::Pointer Import(const std::string& filename) const;

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiNURBSPatchGeoImporter";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:

    void ReadV06Single(std::ifstream& infile,
        std::vector<std::size_t>& orders,
        std::vector<std::size_t>& numbers,
        std::vector<std::vector<double> >& knots,
        std::vector<std::vector<double> >& wcoords,
        std::vector<double>& weights) const;

    void ReadV07Single(std::ifstream& infile,
        std::vector<std::size_t>& orders,
        std::vector<std::size_t>& numbers,
        std::vector<std::vector<double> >& knots,
        std::vector<std::vector<double> >& wcoords,
        std::vector<double>& weights) const;

    void ReadV21Single(std::ifstream& infile,
        std::vector<std::size_t>& orders,
        std::vector<std::size_t>& numbers,
        std::vector<std::vector<double> >& knots,
        std::vector<std::vector<double> >& wcoords,
        std::vector<double>& weights) const;

    void ReadV21Multi(std::ifstream& infile,
        std::vector<std::vector<std::size_t> >& orders,
        std::vector<std::vector<std::size_t> >& numbers,
        std::vector<std::vector<std::vector<double> > >& knots,
        std::vector<std::vector<std::vector<double> > >& wcoords,
        std::vector<std::vector<double> >& weights,
        std::vector<GeoInterface>& interfaces) const;

    void ReadPatchData(std::ifstream& infile,
        const int& rdim,
        std::vector<std::size_t>& orders,
        std::vector<std::size_t>& numbers,
        std::vector<std::vector<double> >& knots,
        std::vector<std::vector<double> >& wcoords,
        std::vector<double>& weights) const;

    typename Patch<TDim>::Pointer CreateNewPatch(const std::size_t& Id,
        const std::vector<std::size_t>& orders,
        const std::vector<std::size_t>& numbers,
        const std::vector<std::vector<double> >& knots,
        const std::vector<std::vector<double> >& wcoords,
        const std::vector<double>& weights) const;
};

template<>
struct BoundarySideHelper<2>
{
    static BoundarySide Get(const int& i)
    {
        switch(i)
        {
            case 1: return _BLEFT_;
            case 2: return _BRIGHT_;
            case 3: return _BBOTTOM_;
            case 4: return _BTOP_;
            default: KRATOS_THROW_ERROR(std::logic_error, i, "is not a valid side");
        }
        return _NUMBER_OF_BOUNDARY_SIDE;
    }
};

template<>
struct BoundarySideHelper<3>
{
    static BoundarySide Get(const int& i)
    {
        switch(i)
        {
            case 1: return _BLEFT_;
            case 2: return _BRIGHT_;
            case 3: return _BFRONT_;
            case 4: return _BBACK_;
            case 5: return _BBOTTOM_;
            case 6: return _BTOP_;
            default: KRATOS_THROW_ERROR(std::logic_error, i, "is not a valid side");
        }
        return _NUMBER_OF_BOUNDARY_SIDE;
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

