#include "custom_utilities/nurbs/domain_manager_3d.h"

namespace Kratos
{

    void DomainManager3D::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DomainManager in 3D";
    }

    void DomainManager3D::PrintData(std::ostream& rOStream) const
    {
        rOStream << "X-coordinates:";
        for(coords_container_t::const_iterator it = BaseType::mXcoords.begin(); it != BaseType::mXcoords.end(); ++it)
            rOStream << " " << *it;
        rOStream << std::endl;

        rOStream << "Y-coordinates:";
        for(coords_container_t::const_iterator it = BaseType::mYcoords.begin(); it != BaseType::mYcoords.end(); ++it)
            rOStream << " " << *it;
        rOStream << std::endl;

        rOStream << "Z-coordinates:";
        for(coords_container_t::const_iterator it = BaseType::mZcoords.begin(); it != BaseType::mZcoords.end(); ++it)
            rOStream << " " << *it;
        rOStream << std::endl;

        rOStream << "Cells:" << std::endl;
        for(map_t::const_iterator it1 = mActiveCells.begin(); it1 != mActiveCells.end(); ++it1)
        {
            rOStream << " face " << it1->first.first << "," << it1->first.second << ":";
            for(std::set<std::size_t>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
                rOStream << " " << *it2;
            rOStream << std::endl;
        }
    }

    std::size_t DomainManager3D::GetIndex(std::size_t X, std::size_t Y, std::size_t Z, std::size_t numX, std::size_t numY, std::size_t numZ) const
    {
        return (Z * numY + Y) * numX + X;
    }

    void DomainManager3D::AddCell(const std::vector<double>& box)
    {
        coords_container_t::iterator it_x1 = BaseType::mXcoords.find(box[0]);
        coords_container_t::iterator it_x2 = BaseType::mXcoords.find(box[1]);
        coords_container_t::iterator it_y1 = BaseType::mYcoords.find(box[2]);
        coords_container_t::iterator it_y2 = BaseType::mYcoords.find(box[3]);
        coords_container_t::iterator it_z1 = BaseType::mZcoords.find(box[4]);
        coords_container_t::iterator it_z2 = BaseType::mZcoords.find(box[5]);

        if(it_x1 == BaseType::mXcoords.end() || it_x2 == BaseType::mXcoords.end())
            KRATOS_THROW_ERROR(std::runtime_error, "Cell does not align with x-coordinates", "")

        if(it_y1 == BaseType::mYcoords.end() || it_y2 == BaseType::mYcoords.end())
            KRATOS_THROW_ERROR(std::runtime_error, "Cell does not align with y-coordinates", "")

        if(it_z1 == BaseType::mZcoords.end() || it_z2 == BaseType::mZcoords.end())
            KRATOS_THROW_ERROR(std::runtime_error, "Cell does not align with z-coordinates", "")

        std::size_t i1 = std::distance(BaseType::mXcoords.begin(), it_x1);
        std::size_t i2 = std::distance(BaseType::mXcoords.begin(), it_x2);
        std::size_t j1 = std::distance(BaseType::mYcoords.begin(), it_y1);
        std::size_t j2 = std::distance(BaseType::mYcoords.begin(), it_y2);
        std::size_t k1 = std::distance(BaseType::mZcoords.begin(), it_z1);
        std::size_t k2 = std::distance(BaseType::mZcoords.begin(), it_z2);

        for(std::size_t i = i1; i < i2; ++i)
        {
            for(std::size_t j = j1; j < j2; ++j)
            {
                for(std::size_t k = k1; k < k2; ++k)
                {
                    key_t key(i, j);
                    mActiveCells[key].insert(k);
                }
            }
        }
    }

    bool DomainManager3D::IsInside(const std::vector<double>& bounding_box) const
    {
        double tol = 1.0e-10;

        // find the lower bound for the Xmin and upper bound for Xmax
        std::size_t i1 = 0, i2 = 0;
        for(coords_container_t::iterator it = BaseType::mXcoords.begin(); it != BaseType::mXcoords.end(); ++it)
            if(bounding_box[0] > *it - tol)
                ++i1;
        if(i1 == 0 || i1 == BaseType::mXcoords.size())
            return false;
        else
            --i1;

        //
        for(coords_container_t::iterator it = BaseType::mXcoords.begin(); it != BaseType::mXcoords.end(); ++it)
            if(bounding_box[1] > *it + tol)
                ++i2;
        if(i2 == 0 || i2 == BaseType::mXcoords.size())
            return false;

        // find the lower bound for the Ymin and upper bound for Ymax
        std::size_t j1 = 0, j2 = 0;
        for(coords_container_t::iterator it = BaseType::mYcoords.begin(); it != BaseType::mYcoords.end(); ++it)
            if(bounding_box[2] > *it - tol)
                ++j1;
        if(j1 == 0 || j1 == BaseType::mYcoords.size())
            return false;
        else
            --j1;

        //
        for(coords_container_t::iterator it = BaseType::mYcoords.begin(); it != BaseType::mYcoords.end(); ++it)
            if(bounding_box[3] > *it + tol)
                ++j2;
        if(j2 == 0 || j2 == BaseType::mYcoords.size())
            return false;

        // find the lower bound for the Zmin and upper bound for Zmax
        std::size_t k1 = 0, k2 = 0;
        for(coords_container_t::iterator it = BaseType::mZcoords.begin(); it != BaseType::mZcoords.end(); ++it)
            if(bounding_box[4] > *it - tol)
                ++k1;
        if(k1 == 0 || k1 == BaseType::mZcoords.size())
            return false;
        else
            --k1;

        //
        for(coords_container_t::iterator it = BaseType::mZcoords.begin(); it != BaseType::mZcoords.end(); ++it)
            if(bounding_box[5] > *it + tol)
                ++k2;
        if(k2 == 0 || k2 == BaseType::mZcoords.size())
            return false;


        // check if in the bound if all the cells are active
        for(std::size_t i = i1; i < i2; ++i)
        {
            for(std::size_t j = j1; j < j2; ++j)
            {
                key_t key(i, j);
                map_t::const_iterator it1 = mActiveCells.find(key);
                if(it1 == mActiveCells.end())
                    return false;

                for(std::size_t k = k1; k < k2; ++k)
                {
                    index_container_t::iterator it2 = it1->second.find(k);
                    if(it2 == it1->second.end())
                        return false;
                }
            }
        }

        return true;
    }

    void DomainManager3D::ExportDomain(const std::string& fn, const std::string& color, const double& distance) const
    {
        std::ofstream outfile(fn.c_str(), std::ios_base::app);

        outfile << "verts = zeros(" << BaseType::mXcoords.size() * BaseType::mYcoords.size() * BaseType::mZcoords.size() << ",3);\n";
        std::size_t i, j, k = 0, row;
        std::size_t num_X_coords = BaseType::mXcoords.size();
        std::size_t num_Y_coords = BaseType::mYcoords.size();
        std::size_t num_Z_coords = BaseType::mZcoords.size();
        for(coords_container_t::const_iterator it_z = BaseType::mZcoords.begin(); it_z != BaseType::mZcoords.end(); ++it_z)
        {
            j = 0;
            for(coords_container_t::const_iterator it_y = BaseType::mYcoords.begin(); it_y != BaseType::mYcoords.end(); ++it_y)
            {
                i = 0;
                for(coords_container_t::const_iterator it_x = BaseType::mXcoords.begin(); it_x != BaseType::mXcoords.end(); ++it_x)
                {
                    row = GetIndex(i, j, k, num_X_coords, num_Y_coords, num_Z_coords) + 1;
                    outfile << "verts(" << row << ",1) = " << *it_x + distance << ";\n";
                    outfile << "verts(" << row << ",2) = " << *it_y << ";\n";
                    outfile << "verts(" << row << ",3) = " << *it_z << ";\n";
                    ++i;
                }
                ++j;
            }
            ++k;
        }

//        this->PrintData(std::cout);

        std::size_t Xi, Yi, Zi, start;
        outfile << "faces = [";
        for(map_t::const_iterator it = mActiveCells.begin(); it != mActiveCells.end(); ++it)
        {
            Xi = it->first.first;
            Yi = it->first.second;
            for(index_container_t::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
            {
                Zi = *it2;
                std::size_t faces[][4][3] = {
                            { {Xi, Yi, Zi}, {Xi+1, Yi, Zi}, {Xi+1, Yi+1, Zi}, {Xi, Yi+1, Zi} },
                            { {Xi, Yi, Zi+1}, {Xi+1, Yi, Zi+1}, {Xi+1, Yi+1, Zi+1}, {Xi, Yi+1, Zi+1} },
                            { {Xi, Yi, Zi}, {Xi+1, Yi, Zi}, {Xi+1, Yi, Zi+1}, {Xi, Yi, Zi+1} },
                            { {Xi+1, Yi, Zi}, {Xi+1, Yi+1, Zi}, {Xi+1, Yi+1, Zi+1}, {Xi+1, Yi, Zi+1} },
                            { {Xi+1, Yi+1, Zi}, {Xi, Yi+1, Zi}, {Xi, Yi+1, Zi+1}, {Xi+1, Yi+1, Zi+1} },
                            { {Xi, Yi+1, Zi}, {Xi, Yi, Zi}, {Xi, Yi, Zi+1}, {Xi, Yi+1, Zi+1} }
                            };

                for(unsigned int i = 0; i < 6; ++i)
                {
                    for(unsigned int j = 0; j < 4; ++j)
                        outfile << " " << GetIndex(faces[i][j][0], faces[i][j][1], faces[i][j][2], num_X_coords, num_Y_coords, num_Z_coords) + 1;
                    outfile << ";\n";
                }
            }
        }
        outfile << "];\n";

        outfile << "patch('Faces',faces,'Vertices',verts,'FaceColor'," << color << ",'FaceAlpha',0.0);\n\n";

        outfile.close();
    }

}

