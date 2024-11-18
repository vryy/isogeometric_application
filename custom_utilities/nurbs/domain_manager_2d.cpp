#include "custom_utilities/nurbs/domain_manager_2d.h"

namespace Kratos
{

void DomainManager2D::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "DomainManager in 2D";
}

void DomainManager2D::PrintData(std::ostream& rOStream) const
{
    rOStream << "X-coordinates:";
    for (coords_container_t::const_iterator it = BaseType::mXcoords.begin(); it != BaseType::mXcoords.end(); ++it)
    {
        rOStream << " " << *it;
    }
    rOStream << std::endl;

    rOStream << "Y-coordinates:";
    for (coords_container_t::const_iterator it = BaseType::mYcoords.begin(); it != BaseType::mYcoords.end(); ++it)
    {
        rOStream << " " << *it;
    }
    rOStream << std::endl;

    rOStream << "Cells:" << std::endl;
    for (map_t::const_iterator it1 = mActiveCells.begin(); it1 != mActiveCells.end(); ++it1)
    {
        rOStream << " column " << it1->first << ":";
        for (index_container_t::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
        {
            rOStream << " " << *it2;
        }
        rOStream << std::endl;
    }
}

void DomainManager2D::AddCell(const std::vector<double>& box)
{
    coords_container_t::iterator it_x1 = BaseType::mXcoords.find(box[0]); //Xmin
    coords_container_t::iterator it_x2 = BaseType::mXcoords.find(box[1]); //Xmax
    coords_container_t::iterator it_y1 = BaseType::mYcoords.find(box[2]); //Ymin
    coords_container_t::iterator it_y2 = BaseType::mYcoords.find(box[3]); //Ymax

    if (it_x1 == BaseType::mXcoords.end() || it_x2 == BaseType::mXcoords.end())
        KRATOS_ERROR << "Cell does not align with x-coordinates";

    if (it_y1 == BaseType::mYcoords.end() || it_y2 == BaseType::mYcoords.end())
        KRATOS_ERROR << "Cell does not align with y-coordinates";

    std::size_t i1 = std::distance(BaseType::mXcoords.begin(), it_x1);
    std::size_t i2 = std::distance(BaseType::mXcoords.begin(), it_x2);
    std::size_t j1 = std::distance(BaseType::mYcoords.begin(), it_y1);
    std::size_t j2 = std::distance(BaseType::mYcoords.begin(), it_y2);

    for (std::size_t i = i1; i < i2; ++i)
    {
        for (std::size_t j = j1; j < j2; ++j)
        {
            mActiveCells[i].insert(j);
        }
    }
}

bool DomainManager2D::IsInside(const std::vector<double>& bounding_box) const
{
    double tol = 1.0e-10;

    // find the lower bound for the Xmin and upper bound for Xmax
    std::size_t i1 = 0, i2 = 0;
    for (coords_container_t::iterator it = BaseType::mXcoords.begin(); it != BaseType::mXcoords.end(); ++it)
        if (bounding_box[0] > *it - tol)
        {
            ++i1;
        }
    if (i1 == 0 || i1 == BaseType::mXcoords.size())
    {
        return false;
    }
    else
    {
        --i1;
    }

    //
    for (coords_container_t::iterator it = BaseType::mXcoords.begin(); it != BaseType::mXcoords.end(); ++it)
        if (bounding_box[1] > *it + tol)
        {
            ++i2;
        }
    if (i2 == 0 || i2 == BaseType::mXcoords.size())
    {
        return false;
    }

    // find the lower bound for the Ymin and upper bound for Ymax
    std::size_t j1 = 0, j2 = 0;
    for (coords_container_t::iterator it = BaseType::mYcoords.begin(); it != BaseType::mYcoords.end(); ++it)
        if (bounding_box[2] > *it - tol)
        {
            ++j1;
        }
    if (j1 == 0 || j1 == BaseType::mYcoords.size())
    {
        return false;
    }
    else
    {
        --j1;
    }

    //
    for (coords_container_t::iterator it = BaseType::mYcoords.begin(); it != BaseType::mYcoords.end(); ++it)
        if (bounding_box[3] > *it + tol)
        {
            ++j2;
        }
    if (j2 == 0 || j2 == BaseType::mYcoords.size())
    {
        return false;
    }

//        KRATOS_WATCH(i1)
//        KRATOS_WATCH(i2)
//        KRATOS_WATCH(j1)
//        KRATOS_WATCH(j2)

    // check if in the bound if all the cells are active
    for (std::size_t i = i1; i < i2; ++i)
    {
        map_t::const_iterator it1 = mActiveCells.find(i);
        if (it1 == mActiveCells.end())
        {
            return false;
        }

        for (std::size_t j = j1; j < j2; ++j)
        {
            index_container_t::iterator it2 = it1->second.find(j);
            if (it2 == it1->second.end())
            {
                return false;
            }
        }
    }

    return true;
}

void DomainManager2D::ExportDomain(const std::string& fn, const std::string& color, double distance) const
{
    std::ofstream outfile(fn.c_str(), std::ios_base::app);

    outfile << "verts = zeros(" << BaseType::mXcoords.size() * BaseType::mYcoords.size() << ",3);\n";
    std::size_t i = 0, j = 0, row;
    for (coords_container_t::const_iterator it_y = BaseType::mYcoords.begin(); it_y != BaseType::mYcoords.end(); ++it_y)
    {
        for (coords_container_t::const_iterator it_x = BaseType::mXcoords.begin(); it_x != BaseType::mXcoords.end(); ++it_x)
        {
            row = j * BaseType::mXcoords.size() + i + 1;
            outfile << "verts(" << row << ",1) = " << *it_x << ";\n";
            outfile << "verts(" << row << ",2) = " << *it_y << ";\n";
            outfile << "verts(" << row << ",3) = " << distance << ";\n";
            ++i;
        }
        ++j;
        i = 0;
    }

//        this->PrintData(std::cout);

    std::size_t cnt = 0;
    std::size_t start, number_of_patches = 0;
    for (map_t::const_iterator it = mActiveCells.begin(); it != mActiveCells.end(); ++it)
    {
        number_of_patches += it->second.size();
    }
    outfile << "faces = zeros(" << number_of_patches << ",4);\n";
    for (map_t::const_iterator it = mActiveCells.begin(); it != mActiveCells.end(); ++it)
    {
        for (index_container_t::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
            start = *it2 * mXcoords.size() + it->first + 1;
            outfile << "faces(" << ++cnt << ",:) = [" << start << " " << (start + 1) << " " << (start + 1 + mXcoords.size()) << " " << (start + mXcoords.size()) << "];\n";
        }
    }

    outfile << "patch('Faces',faces,'Vertices',verts,'FaceColor'," << color << ");\n\n";

    outfile.close();
}

}
