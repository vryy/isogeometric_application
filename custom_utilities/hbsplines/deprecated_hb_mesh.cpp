#include <fstream>
#include <algorithm>

#include "custom_utilities/hbsplines/deprecated_hb_mesh.h"

#include "utilities/openmp_utils.h"
#include "custom_utilities/bspline_utils.h"
#include "custom_utilities/bezier_utils.h"
#include "custom_utilities/nurbs/cell_manager_2d.h"
#include "custom_utilities/nurbs/cell_manager_3d.h"
#include "custom_utilities/triangulation_utils.h"
#include "custom_utilities/isogeometric_math_utils.h"
#include "utilities/auto_collapse_spatial_binning.h"

#ifdef ISOGEOMETRIC_USE_TETGEN
#define TETLIBRARY
#include "custom_external_libraries/tetgen1.5.0/tetgen.h"
#endif

#define ENABLE_PROFILING
#define DEBUG_REFINE

namespace Kratos
{

    template<int TDim>
    DeprecatedHBMesh<TDim>::DeprecatedHBMesh(const std::size_t& Id, const std::string& Name)
    : BaseType(Id), mName(Name), mEchoLevel(0), mLastLevel(0), mMaxLevels(10)
    {
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "hierarchical B-Splines mesh " << Name() << ", Id = " << BaseType::Id() << ", number of levels = " << mLastLevel;
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::PrintData(std::ostream& rOStream) const
    {
        BaseType::PrintData(rOStream);
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::PrintKnotVectors() const
    {
        std::cout << "###############Begin knot vectors################" << std::endl;
        std::cout << "knot vector 1:" << mKnots1 << std::endl;
        std::cout << "knot vector 2:" << mKnots2 << std::endl;
        std::cout << "knot vector 3:" << mKnots3 << std::endl;
        std::cout << "###############End knot vectors##################" << std::endl;
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::PrintCells(int level) const
    {
        if(level > 0)
        {
            std::cout << "###############Begin cells at level " << level << "################" << std::endl;
            std::size_t n = 0;
            for(typename cell_container_t::iterator it = mpCellManager->begin(); it != mpCellManager->end(); ++it)
                if((*it)->Level() == level)
                    std::cout << "(" << ++n << ") " << *(*it) << std::endl;
            std::cout << "###############End cells at level " << level << "################" << std::endl;
        }
        else
        {
            std::cout << "###############Begin cells at all levels" << "################" << std::endl;
            std::size_t n = 0;
            for(typename cell_container_t::iterator it = mpCellManager->begin(); it != mpCellManager->end(); ++it)
                std::cout << "(" << ++n << ") " << *(*it) << std::endl;
            std::cout << "###############End cells at levels" << "################" << std::endl;
        }
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::PrintBasisFuncs(int level) const
    {
        if(level > 0)
        {
            std::cout << "###############Begin basis functions at level " << level << "################" << std::endl;
            std::size_t n = 0;
            for(typename bf_container_t::iterator it = mBasisFuncs.begin(); it != mBasisFuncs.end(); ++it)
                if((*it)->Level() == level)
                    std::cout << "(" << ++n << ") " << *(*it) << std::endl;
            std::cout << "###############End basis functions at level " << level << "##################" << std::endl;
        }
        else
        {
            std::cout << "###############Begin basis functions at all levels" << "################" << std::endl;
            std::size_t n = 0;
            for(typename bf_container_t::iterator it = mBasisFuncs.begin(); it != mBasisFuncs.end(); ++it)
                std::cout << "(" << ++n << ") " << *(*it) << std::endl;
            std::cout << "###############End basis functions at all levels" << "##################" << std::endl;
        }
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::PrintRefinementHistory() const
    {
        std::cout << "Refinement history:";
        for(std::size_t i = 0; i < mRefinementHistory.size(); ++i)
            std::cout << ", " << mRefinementHistory[i];
        std::cout << std::endl;
    }

    // template<int TDim>
//    void DeprecatedHBMesh<TDim>::CheckNestedSpace()
//    {
//        std::cout << __FUNCTION__ << " starts" << std::endl;
//        double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;
//        for(bf_container_t::iterator it_bf = mBasisFuncs.begin(); it_bf != mBasisFuncs.end(); ++it_bf)
//        {
//            // extract the support domain of the next level
//            std::size_t level = (*it_bf)->Level();
//            domain_t p_domain = GetSupportDomain(level + 1);
//
//            // get the support domain of the bf
//            (*it_bf)->GetBoundingBox(Xmin, Xmax, Ymin, Ymax, Zmin, Zmax);
//
//            // check if the bf support domain contained in the refined domain managed by the domain manager
//            std::size_t subdomain_id = p_domain->IsInside(Xmin, Xmax, Ymin, Ymax, Zmin, Zmax);
//            if(subdomain_id != 0)
//                std::cout << "Bf " << (*it_bf)->Id() << " of level " << level << " has the support lie in the support domain of next level " << (level + 1) << std::endl;
//        }
//        std::cout << __FUNCTION__ << " completed" << std::endl;
//    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::BuildNestedSpace(std::size_t level, std::map<std::size_t, std::set<std::size_t> >& rK)
    {
        // check the first criteria
        if(level == 1)
        {
            for(bf_container_t::iterator it_bf = mBasisFuncs.begin(); it_bf != mBasisFuncs.end(); ++it_bf)
            {
                if((*it_bf)->Level() == level)
                    rK[level].insert((*it_bf)->Id());
            }
        }
        else
        {
            double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;
            for(std::set<std::size_t>::iterator it = rK[level - 1].begin(); it != rK[level - 1].end(); ++it)
            {
                mBasisFuncs[*it]->GetBoundingBox(Xmin, Xmax, Ymin, Ymax, Zmin, Zmax);

                bool is_inside;
                if(TDim == 2)
                {
                    std::vector<double> bounding_box = {Xmin, Xmax, Ymin, Ymax};
                    is_inside = GetSupportDomain(level)->IsInside(bounding_box);
                }
                else if(TDim == 3)
                {
                    std::vector<double> bounding_box = {Xmin, Xmax, Ymin, Ymax, Zmin, Zmax};
                    is_inside = GetSupportDomain(level)->IsInside(bounding_box);
                }

                if(is_inside)
                    rK[level].insert(*it);
            }
//            KRATOS_WATCH(aux)
        }

        // check the second criteria
        if(level > 1)
            for(bf_container_t::iterator it_bf = mBasisFuncs.begin(); it_bf != mBasisFuncs.end(); ++it_bf)
            {
                if((*it_bf)->Level() == level)
                    rK[level].insert((*it_bf)->Id());
            }
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::CheckNestedSpace()
    {
        std::cout << __FUNCTION__ << " starts" << std::endl;

        std::map<std::size_t, std::set<std::size_t> > K;
        for(std::size_t level = 1; level <= mLastLevel; ++level)
            this->BuildNestedSpace(level, K);

        for(std::map<std::size_t, std::set<std::size_t> >::iterator it = K.begin(); it != K.end(); ++it)
        {
            std::cout << "K[" << it->first << "] =";
            for(std::set<std::size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
                std::cout << " " << *it2;
            std::cout << std::endl;
        }

        KRATOS_WATCH(mBasisFuncs.size())
        KRATOS_WATCH(K[mLastLevel].size())

//        std::cout << "mBasisFuncs:";
//        for(bf_container_t::iterator it_bf = mBasisFuncs.begin(); it_bf != mBasisFuncs.end(); ++it_bf)
//            std::cout << " " << (*it_bf)->Id();
//        std::cout << std::endl;

        // check if any bf in level 3 lie completely in level 5
        unsigned int level_to_compare = 3;
        unsigned int base_level = 4;
        domain_t p_domain = GetSupportDomain(base_level);
        double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;
        for(bf_container_t::iterator it_bf = mBasisFuncs.begin(); it_bf != mBasisFuncs.end(); ++it_bf)
        {
            if((*it_bf)->Level() == level_to_compare)
            {
                (*it_bf)->GetBoundingBox(Xmin, Xmax, Ymin, Ymax, Zmin, Zmax);

                bool is_inside;
                if(TDim == 2)
                {
                    std::vector<double> bounding_box = {Xmin, Xmax, Ymin, Ymax};
                    is_inside = p_domain->IsInside(bounding_box);
                }
                else if(TDim == 3)
                {
                    std::vector<double> bounding_box = {Xmin, Xmax, Ymin, Ymax, Zmin, Zmax};
                    is_inside = p_domain->IsInside(bounding_box);
                }

                if(is_inside)
                    std::cout << "bf " << (*it_bf)->Id() << " of level " << level_to_compare << " lie completely in level " << base_level << std::endl;
            }
        }

        unsigned int sample_id = 1;
        KRATOS_WATCH(mBasisFuncs[sample_id]->Level())
        mBasisFuncs[sample_id]->GetBoundingBox(Xmin, Xmax, Ymin, Ymax, Zmin, Zmax);
        std::cout << "bf " << sample_id << " support domain: " << Xmin << " " << Xmax << " " << Ymin << " " << Ymax << " " << Zmin << " " << Zmax << std::endl;
        KRATOS_WATCH(GetSupportDomain(2)->IsInside(std::vector<double>{Xmin, Xmax, Ymin, Ymax}))
        KRATOS_WATCH(*GetSupportDomain(2))

//        unsigned int bfs[] = {201, 202, 205, 206, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 308};
//        for(unsigned int i = 0; i < 17; ++i)
//            std::cout << "bf " << bfs[i] << " level: " << mBasisFuncs[bfs[i]]->Level() << std::endl;

        std::cout << __FUNCTION__ << " completed" << std::endl;
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::ReadMesh(const std::string& fn)
    {
        std::ifstream infile(fn.c_str());
        if(!infile)
            KRATOS_THROW_ERROR(std::logic_error, "Error open file", fn)

        std::string line;
        std::vector<std::string> words;
        int read_mode = READ_PATCH;
        int npatches, dim_index = 0;
        std::vector<int> orders;
        std::vector<int> numbers;
        std::vector<double> x_coords;
        std::vector<double> y_coords;
        std::vector<double> z_coords;
        std::vector<double> weights;
        while(!infile.eof())
        {
            std::getline(infile, line);
            boost::trim_if(line, boost::is_any_of("\t ")); // ignore trailing spaces
            boost::split(words, line, boost::is_any_of(" \t"), boost::token_compress_on);

            if(words.size() != 0)
            {
                if(words[0] == std::string("#") || words[0][0] == '#')
                    continue;

                if(read_mode == READ_PATCH)
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
                        KRATOS_THROW_ERROR(std::logic_error, "Currently number of patches > 1 is not supported, npatches =", npatches)
                    }
                    read_mode = READ_ORDER;
                    continue;
                }

                if(read_mode == READ_ORDER)
                {
                    // bound check
                    if(words.size() != TDim)
                        KRATOS_THROW_ERROR(std::logic_error, "The Order section must contained number of information equal to dimension, current number of information =", words.size())

                    // read info
                    for(std::size_t i = 0; i < TDim; ++i)
                        orders.push_back(atoi(words[i].c_str()));
                    read_mode = READ_NUMBER;
                    continue;
                }

                if(read_mode == READ_NUMBER)
                {
                    // bound check
                    if(words.size() != TDim)
                        KRATOS_THROW_ERROR(std::logic_error, "The Number section must contained number of information equal to dimension, current number of information =", words.size())

                    for(std::size_t i = 0; i < TDim; ++i)
                        numbers.push_back(atoi(words[i].c_str()));
                    read_mode = READ_KNOTS;
                    continue;
                }

                if(read_mode == READ_KNOTS)
                {
                    // bound check
                    int knot_len = numbers[dim_index] + orders[dim_index] + 1;
                    if(words.size() != knot_len)
                        KRATOS_THROW_ERROR(std::logic_error, "The Knots section must contained number of information equal to n+p+1, current number of information =", words.size())

                    for(std::size_t i = 0; i < knot_len; ++i)
                    {
                        double k = atof(words[i].c_str());
                        if(dim_index == 0)
                            mKnots1.pCreateKnot(k);
                        else if(dim_index == 1)
                            mKnots2.pCreateKnot(k);
                        else if(dim_index == 2)
                            mKnots3.pCreateKnot(k);
                        else
                            KRATOS_THROW_ERROR(std::logic_error, "Wrong knot dimension index. Something must be wrong", "")
                    }

                    ++dim_index;
                    if(dim_index == TDim)
                    {
                        dim_index = 0;
                        read_mode = READ_COORDINATES;
                    }
                    continue;
                }

                if(read_mode == READ_COORDINATES)
                {
                    // bound check
                    int num_basis = 1;
                    for(std::size_t i = 0; i < TDim; ++i)
                        num_basis *= numbers[i];
                    if(words.size() != num_basis)
                        KRATOS_THROW_ERROR(std::logic_error, "The Coordinates section must contained number of information equal to prod(ni), current number of information =", words.size())

                    if(dim_index == 0)
                        for(std::size_t i = 0; i < num_basis; ++i)
                            x_coords.push_back(atof(words[i].c_str()));
                    else if(dim_index == 1)
                        for(std::size_t i = 0; i < num_basis; ++i)
                            y_coords.push_back(atof(words[i].c_str()));
                    else if(dim_index == 2)
                        for(std::size_t i = 0; i < num_basis; ++i)
                            z_coords.push_back(atof(words[i].c_str()));

                    ++dim_index;
                    if(dim_index == TDim)
                    {
                        dim_index = 0;
                        read_mode = READ_WEIGHTS;
                    }
                    continue;
                }

                if(read_mode == READ_WEIGHTS)
                {
                    // bound check
                    int num_basis = 1;
                    for(std::size_t i = 0; i < TDim; ++i)
                        num_basis *= numbers[i];
                    if(words.size() != num_basis)
                        KRATOS_THROW_ERROR(std::logic_error, "The Weights section must contained number of information equal to prod(ni), current number of information =", words.size())

                    for(std::size_t i = 0; i < num_basis; ++i)
                        weights.push_back(atof(words[i].c_str()));

                    read_mode = NO_READ;
                    continue;
                }
            }
        }

        // close the file
        infile.close();

        if((GetEchoLevel() & ECHO_REFINEMENT) == ECHO_REFINEMENT)
            std::cout << __FUNCTION__ << ": Traverse file completed" << std::endl;

        // update the order
        if(orders.size() > 0)
            mOrder1 = orders[0];
        if(orders.size() > 1)
            mOrder2 = orders[1];
        if(orders.size() > 2)
            mOrder3 = orders[2];

        // rescale the coordinates
        for(std::size_t i = 0; i < weights.size(); ++i)
        {
            x_coords[i] /= weights[i];
            y_coords[i] /= weights[i];
            if(TDim > 2)
                z_coords[i] /= weights[i];
        }

        // initialize the cell container
        if(TDim == 2)
            mpCellManager = cell_container_t::Pointer(new CellManager2D<CellType>());
        else if(TDim == 3)
            mpCellManager = cell_container_t::Pointer(new CellManager3D<CellType>());

        // create bfs for the first level
        unsigned int lastID = 0;
        unsigned int level = 1;
        mLastLevel = 1;
        double area_tol = 1.0e-6; // tolerance to accept the nonzero-area cell. We should parameterize it. TODO
        if(TDim == 2)
        {
            int num1 = mKnots1.size() - mOrder1 - 1;
            int num2 = mKnots2.size() - mOrder2 - 1;
            for(std::size_t j = 0; j < num2; ++j)
            {
                // create and fill the local knot vector
                std::vector<knot_t> pLocalKnots2;
                for(std::size_t k = 0; k < mOrder2 + 2; ++k)
                    pLocalKnots2.push_back(mKnots2.pKnotAt(j + k));

                for(std::size_t i = 0; i < num1; ++i)
                {
                    // create and fill the local knot vector
                    std::vector<knot_t> pLocalKnots1;
                    for(std::size_t k = 0; k < mOrder1 + 2; ++k)
                        pLocalKnots1.push_back(mKnots1.pKnotAt(i + k));

                    // create the basis function object
                    std::vector<knot_t> pLocalKnots3;
                    DeprecatedHBBasisFunction::Pointer p_bf = mBasisFuncs.CreateBf(level, pLocalKnots1, pLocalKnots2, pLocalKnots3);

                    // assign the coordinates and weight
                    int ifunc = j * num1 + i;
                    p_bf->GetControlPoint().SetCoordinates(x_coords[ifunc], y_coords[ifunc], 0.0, weights[ifunc]);

                    // create the cells for the basis function
                    for(std::size_t i1 = 0; i1 < mOrder1 + 1; ++i1)
                    {
                        knot_t pLeft = mKnots1.pKnotAt(i + i1);
                        knot_t pRight = mKnots1.pKnotAt(i + i1 + 1);
                        for(std::size_t j1 = 0; j1 < mOrder2 + 1; ++j1)
                        {
                            knot_t pDown = mKnots2.pKnotAt(j + j1);
                            knot_t pUp = mKnots2.pKnotAt(j + j1 + 1);

                            // check if the cell domain area is nonzero
                            double area = (pRight->Value() - pLeft->Value()) * (pUp->Value() - pDown->Value());
                            if(fabs(area) > area_tol)
                            {
                                std::vector<knot_t> pKnots = {pLeft, pRight, pDown, pUp};
                                cell_t p_cell = mpCellManager->CreateCell(pKnots);
                                p_cell->SetLevel(level);
                                p_bf->AddCell(p_cell);
                                p_cell->AddBf(p_bf);
                            }
                        }
                    }
                }
            }
        }
        else if(TDim == 3)
        {
            int num1 = mKnots1.size() - mOrder1 - 1;
            int num2 = mKnots2.size() - mOrder2 - 1;
            int num3 = mKnots3.size() - mOrder3 - 1;
            for(std::size_t l = 0; l < num3; ++l)
            {
                // create and fill the local knot vector
                std::vector<knot_t> pLocalKnots3;
                for(std::size_t k = 0; k < mOrder3 + 2; ++k)
                    pLocalKnots3.push_back(mKnots3.pKnotAt(l + k));

                for(std::size_t j = 0; j < num2; ++j)
                {
                    // create and fill the local knot vector
                    std::vector<knot_t> pLocalKnots2;
                    for(std::size_t k = 0; k < mOrder2 + 2; ++k)
                        pLocalKnots2.push_back(mKnots2.pKnotAt(j + k));

                    for(std::size_t i = 0; i < num1; ++i)
                    {
                        // create and fill the local knot vector
                        std::vector<knot_t> pLocalKnots1;
                        for(std::size_t k = 0; k < mOrder1 + 2; ++k)
                            pLocalKnots1.push_back(mKnots1.pKnotAt(i + k));

                        // create the basis function object
                        DeprecatedHBBasisFunction::Pointer p_bf = mBasisFuncs.CreateBf(level, pLocalKnots1, pLocalKnots2, pLocalKnots3);

                        // assign the coordinates and weight
                        int ifunc = (l * num2 + j) * num1 + i;
                        p_bf->GetControlPoint().SetCoordinates(x_coords[ifunc], y_coords[ifunc], z_coords[ifunc], weights[ifunc]);

                        // create the cells for the basis function
                        for(std::size_t i1 = 0; i1 < mOrder1 + 1; ++i1)
                        {
                            knot_t pLeft = mKnots1.pKnotAt(i + i1);
                            knot_t pRight = mKnots1.pKnotAt(i + i1 + 1);
                            for(std::size_t j1 = 0; j1 < mOrder2 + 1; ++j1)
                            {
                                knot_t pDown = mKnots2.pKnotAt(j + j1);
                                knot_t pUp = mKnots2.pKnotAt(j + j1 + 1);
                                for(std::size_t l1 = 0; l1 < mOrder3 + 1; ++l1)
                                {
                                    knot_t pBelow = mKnots3.pKnotAt(l + l1);
                                    knot_t pAbove = mKnots3.pKnotAt(l + l1 + 1);

                                    // check if the cell domain area is nonzero
                                    double area = (pRight->Value() - pLeft->Value()) * (pUp->Value() - pDown->Value()) * (pAbove->Value() - pBelow->Value());
                                    if(fabs(area) > area_tol)
                                    {
                                        std::vector<knot_t> pKnots = {pLeft, pRight, pDown, pUp, pBelow, pAbove};
                                        cell_t p_cell = mpCellManager->CreateCell(pKnots);
                                        p_cell->SetLevel(level);
                                        p_bf->AddCell(p_cell);
                                        p_cell->AddBf(p_bf);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if((GetEchoLevel() & ECHO_REFINEMENT) == ECHO_REFINEMENT)
            std::cout << "ReadMesh and build level 1 completed" << std::endl;
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::Refine(const std::size_t& Id)
    {
        #ifdef ENABLE_PROFILING
        double start = OpenMPUtils::GetCurrentTime();
        #endif

        bf_t p_bf;
        bool found = false;
        for(typename bf_container_t::iterator it = mBasisFuncs.begin(); it != mBasisFuncs.end(); ++it)
            if((*it)->Id() == Id)
            {
                p_bf = *it;
                found = true;
            }

        if(!found)
            return;

        // does not refine if maximum level is reached
        if(p_bf->Level() == mMaxLevels)
        {
            std::cout << "Maximum level is reached, basis function " << p_bf->Id() << " is skipped" << std::endl;
            return;
        }

        /* create a list of basis function in the next level representing this basis function */
        // create a list of new knots
        double tol = 1.0e-10;
        std::vector<std::vector<knot_t> > pnew_local_knots(3);
        std::vector<std::vector<double> > ins_knots(3);
        for(unsigned int i = 0; i < 3; ++i)
        {
            const std::vector<knot_t>& pLocalKnots = p_bf->LocalKnots(i + 1);

            for(std::vector<knot_t>::const_iterator it = pLocalKnots.begin(); it != pLocalKnots.end(); ++it)
            {
                pnew_local_knots[i].push_back(*it);

                std::vector<knot_t>::const_iterator it2 = it + 1;
                if(it2 != pLocalKnots.end())
                {
                    if(fabs((*it2)->Value() - (*it)->Value()) > tol)
                    {
                        // now we just add the middle one, but in the general we can add arbitrary values
                        // TODO: generalize this
                        double ins_knot = 0.5 * ((*it)->Value() + (*it2)->Value());
                        knot_t p_new_knot;
                        if(i == 0)
                            p_new_knot = mKnots1.pCreateUniqueKnot(ins_knot, tol);
                        else if(i == 1)
                            p_new_knot = mKnots2.pCreateUniqueKnot(ins_knot, tol);
                        else if(i == 2)
                            p_new_knot = mKnots3.pCreateUniqueKnot(ins_knot, tol);

                        pnew_local_knots[i].push_back(p_new_knot);
                        ins_knots[i].push_back(p_new_knot->Value());
                    }
                }
            }
        }

        /* compute the refinement coefficients */
        Vector RefinedCoeffs;
        std::vector<std::vector<double> > local_knots(3);
        for(unsigned int i = 0; i < 3; ++i)
            p_bf->GetLocalKnots(i + 1, local_knots[i]);

        std::vector<std::vector<double> > new_knots(TDim);
        if(TDim == 2)
            BSplineUtils::ComputeBsplinesKnotInsertionCoefficients2DLocal(RefinedCoeffs,
                                                                          new_knots[0],
                                                                          new_knots[1],
                                                                          mOrder1,
                                                                          mOrder2,
                                                                          local_knots[0],
                                                                          local_knots[1],
                                                                          ins_knots[0],
                                                                          ins_knots[1]);
        else if(TDim == 3)
            BSplineUtils::ComputeBsplinesKnotInsertionCoefficients3DLocal(RefinedCoeffs,
                                                                          new_knots[0],
                                                                          new_knots[1],
                                                                          new_knots[2],
                                                                          mOrder1,
                                                                          mOrder2,
                                                                          mOrder3,
                                                                          local_knots[0],
                                                                          local_knots[1],
                                                                          local_knots[2],
                                                                          ins_knots[0],
                                                                          ins_knots[1],
                                                                          ins_knots[2]);

        #ifdef ENABLE_PROFILING
        double time_1 = OpenMPUtils::GetCurrentTime() - start;
        start = OpenMPUtils::GetCurrentTime();
        #endif

        /* create new basis functions */
        unsigned int next_level = p_bf->Level() + 1;
        if(next_level > mLastLevel)
            mLastLevel = next_level;
        double area_tol = 1.0e-6; // tolerance to accept the nonzero-area cell. We should parameterize it.
        typename cell_container_t::Pointer pnew_cells;
        if(TDim == 2)
            pnew_cells = cell_container_t::Pointer(new CellManager2D<CellType>());
        else if(TDim == 3)
            pnew_cells = cell_container_t::Pointer(new CellManager3D<CellType>());
        double father_weight = p_bf->GetControlPoint().W();
        double father_X = p_bf->GetControlPoint().X();
        double father_Y = p_bf->GetControlPoint().Y();
        double father_Z = p_bf->GetControlPoint().Z();
        if(TDim == 2)
        {
            int num1 = pnew_local_knots[0].size() - mOrder1 - 1;
            int num2 = pnew_local_knots[1].size() - mOrder2 - 1;
            for(std::size_t j = 0; j < num2; ++j)
            {
                // create and fill the local knot vector
                std::vector<knot_t> pLocalKnots2;
                for(std::size_t k = 0; k < mOrder2 + 2; ++k)
                    pLocalKnots2.push_back(pnew_local_knots[1][j + k]);

                for(std::size_t i = 0; i < num1; ++i)
                {
                    // create and fill the local knot vector
                    std::vector<knot_t> pLocalKnots1;
                    for(std::size_t k = 0; k < mOrder1 + 2; ++k)
                        pLocalKnots1.push_back(pnew_local_knots[0][i + k]);

                    // create the basis function object
                    std::vector<knot_t> pLocalKnots3;
                    DeprecatedHBBasisFunction::Pointer pnew_bf = mBasisFuncs.CreateBf(next_level, pLocalKnots1, pLocalKnots2, pLocalKnots3);

                    // update the coordinates
                    int ifunc = j * num1 + i;
                    double add_weight = father_weight * RefinedCoeffs[ifunc];
                    pnew_bf->GetControlPoint().AddCoordinates(father_X, father_Y, father_Z, add_weight);

                    // create the cells for the basis function
                    for(std::size_t i1 = 0; i1 < mOrder1 + 1; ++i1)
                    {
                        knot_t pLeft = pnew_local_knots[0][i + i1];
                        knot_t pRight = pnew_local_knots[0][i + i1 + 1];
                        for(std::size_t j1 = 0; j1 < mOrder2 + 1; ++j1)
                        {
                            knot_t pDown = pnew_local_knots[1][j + j1];
                            knot_t pUp = pnew_local_knots[1][j + j1 + 1];

                            // check if the cell domain area is nonzero
                            double area = (pRight->Value() - pLeft->Value()) * (pUp->Value() - pDown->Value());
                            if(fabs(area) > area_tol)
                            {
                                std::vector<knot_t> pKnots = {pLeft, pRight, pDown, pUp};
                                cell_t pnew_cell = mpCellManager->CreateCell(pKnots);
                                pnew_cell->SetLevel(next_level);
                                pnew_bf->AddCell(pnew_cell);
                                pnew_cell->AddBf(pnew_bf);
                                pnew_cells->insert(pnew_cell);
                            }
                        }
                    }
                }
            }
        }
        else if(TDim == 3)
        {
            int num1 = pnew_local_knots[0].size() - mOrder1 - 1;
            int num2 = pnew_local_knots[1].size() - mOrder2 - 1;
            int num3 = pnew_local_knots[2].size() - mOrder3 - 1;
            for(std::size_t l = 0; l < num3; ++l)
            {
                // create and fill the local knot vector
                std::vector<knot_t> pLocalKnots3;
                for(std::size_t k = 0; k < mOrder3 + 2; ++k)
                    pLocalKnots3.push_back(pnew_local_knots[2][l + k]);

                for(std::size_t j = 0; j < num2; ++j)
                {
                    // create and fill the local knot vector
                    std::vector<knot_t> pLocalKnots2;
                    for(std::size_t k = 0; k < mOrder2 + 2; ++k)
                        pLocalKnots2.push_back(pnew_local_knots[1][j + k]);

                    for(std::size_t i = 0; i < num1; ++i)
                    {
                        // create and fill the local knot vector
                        std::vector<knot_t> pLocalKnots1;
                        for(std::size_t k = 0; k < mOrder1 + 2; ++k)
                            pLocalKnots1.push_back(pnew_local_knots[0][i + k]);

                        // create the basis function object
                        DeprecatedHBBasisFunction::Pointer pnew_bf = mBasisFuncs.CreateBf(next_level, pLocalKnots1, pLocalKnots2, pLocalKnots3); //#

                        // update the coordinates
                        int ifunc = (l * num2 + j) * num1 + i;
                        double add_weight = father_weight * RefinedCoeffs[ifunc];
                        pnew_bf->GetControlPoint().AddCoordinates(father_X, father_Y, father_Z, add_weight);

                        // create the cells for the basis function
                        for(std::size_t i1 = 0; i1 < mOrder1 + 1; ++i1)
                        {
                            knot_t pLeft = pnew_local_knots[0][i + i1];
                            knot_t pRight = pnew_local_knots[0][i + i1 + 1];
                            for(std::size_t j1 = 0; j1 < mOrder2 + 1; ++j1)
                            {
                                knot_t pDown = pnew_local_knots[1][j + j1];
                                knot_t pUp = pnew_local_knots[1][j + j1 + 1];
                                for(std::size_t l1 = 0; l1 < mOrder3 + 1; ++l1)
                                {
                                    knot_t pBelow = pnew_local_knots[2][l + l1];
                                    knot_t pAbove = pnew_local_knots[2][l + l1 + 1];

                                    // check if the cell domain area is nonzero
                                    double area = (pRight->Value() - pLeft->Value()) * (pUp->Value() - pDown->Value()) * (pAbove->Value() - pBelow->Value());
                                    if(fabs(area) > area_tol)
                                    {
                                        std::vector<knot_t> pKnots = {pLeft, pRight, pDown, pUp, pBelow, pAbove};
                                        cell_t pnew_cell = mpCellManager->CreateCell(pKnots);
                                        pnew_cell->SetLevel(next_level);
                                        pnew_bf->AddCell(pnew_cell);
                                        pnew_cell->AddBf(pnew_bf);
                                        pnew_cells->insert(pnew_cell);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        #ifdef ENABLE_PROFILING
        double time_2 = OpenMPUtils::GetCurrentTime() - start;
        start = OpenMPUtils::GetCurrentTime();
        #endif

        /* remove the cells of the old basis function (remove only the cell in the current level) */
        typename cell_container_t::Pointer pcells_to_remove;
        if(TDim == 2)
            pcells_to_remove = cell_container_t::Pointer(new CellManager2D<CellType>());
        else if(TDim == 3)
            pcells_to_remove = cell_container_t::Pointer(new CellManager3D<CellType>());

        // firstly we check if the cell c of the current bf in the current level cover any sub-cells. Then the sub-cell includes all bfs of the cell c.
        for(DeprecatedHBBasisFunction::cell_iterator it_cell = p_bf->cell_begin(); it_cell != p_bf->cell_end(); ++it_cell)
        {
            if((*it_cell)->Level() == p_bf->Level())
            {
                for(cell_container_t::iterator it_subcell = pnew_cells->begin(); it_subcell != pnew_cells->end(); ++it_subcell)
                {
                    if((*it_subcell)->template IsCovered<TDim>(*it_cell))
                    {
                        for(typename CellType::bf_iterator it_bf = (*it_cell)->bf_begin(); it_bf != (*it_cell)->bf_end(); ++it_bf)
                        {
                            (*it_subcell)->AddBf(*it_bf);
                            (*it_bf)->AddCell(*it_subcell);
                        }
                    }
                }

                // mark to remove the old cell
                pcells_to_remove->insert(*it_cell);
            }
        }

        // secondly, it happens that new cell c cover several existing cells. In this case cell c must be removed, and its bfs will be transferred to sub-cells.
        for(cell_container_t::iterator it_cell = pnew_cells->begin(); it_cell != pnew_cells->end(); ++it_cell)
        {
            std::vector<cell_t> p_cells = mpCellManager->GetCells(*it_cell);
            if(p_cells.size() > 0)
            {
                if((GetEchoLevel() & ECHO_REFINEMENT) == ECHO_REFINEMENT)
                {
                    std::cout << "cell " << (*it_cell)->Id() << " is detected to contain some smaller cells:";
                    for(std::size_t i = 0; i < p_cells.size(); ++i)
                        std::cout << " " << p_cells[i]->Id();
                    std::cout << std::endl;
                }

                pcells_to_remove->insert(*it_cell);
                for(std::size_t i = 0; i < p_cells.size(); ++i)
                {
                    for(typename CellType::bf_iterator it_bf = (*it_cell)->bf_begin(); it_bf != (*it_cell)->bf_end(); ++it_bf)
                    {
                        p_cells[i]->AddBf(*it_bf);
                        (*it_bf)->AddCell(p_cells[i]);
                    }
                }
            }
        }

        /* remove the cells */
        for(cell_container_t::iterator it_cell = pcells_to_remove->begin(); it_cell != pcells_to_remove->end(); ++it_cell)
        {
            mpCellManager->erase(*it_cell);
            for(bf_container_t::iterator it_bf = mBasisFuncs.begin(); it_bf != mBasisFuncs.end(); ++it_bf)
                (*it_bf)->RemoveCell(*it_cell);
        }

        /* remove the basis function from all the cells */
        for(cell_container_t::iterator it_cell = mpCellManager->begin(); it_cell != mpCellManager->end(); ++it_cell)
            (*it_cell)->RemoveBf(p_bf);

        /* remove the old basis function */
        mBasisFuncs.erase(p_bf);

        #ifdef ENABLE_PROFILING
        double time_3 = OpenMPUtils::GetCurrentTime() - start;
        start = OpenMPUtils::GetCurrentTime();
        #endif

        /* Debug part: to be removed when the code is stable
        for(cell_container_t::iterator it_cell = mpCellManager->begin(); it_cell != mpCellManager->end(); ++it_cell)
        {
            std::cout << "cell " << (*it_cell)->Id() << " supports:";
            for(typename CellType::bf_iterator it_bf = (*it_cell)->bf_begin(); it_bf != (*it_cell)->bf_end(); ++it_bf)
                std::cout << " " << (*it_bf)->Id();
            std::cout << std::endl;
        }

        for(bf_container_t::iterator it_bf = mBasisFuncs.begin(); it_bf != mBasisFuncs.end(); ++it_bf)
        {
            std::cout << "bf " << (*it_bf)->Id() << " contains cell:";
            for(DeprecatedHBBasisFunction::cell_iterator it_cell = (*it_bf)->cell_begin(); it_cell != (*it_bf)->cell_end(); ++it_cell)
                std::cout << " " << (*it_cell)->Id();
            std::cout << std::endl;
        }
        */

        mRefinementHistory.push_back(p_bf->Id());
        if((GetEchoLevel() & ECHO_REFINEMENT) == ECHO_REFINEMENT)
        {
        #ifdef ENABLE_PROFILING
            std::cout << "Refine bf " << p_bf->Id() << " completed" << std::endl;
            std::cout << " Time to compute the refinement coefficients: " << time_1 << " s" << std::endl;
            std::cout << " Time to create new cells and new bfs: " << time_2 << " s" << std::endl;
            std::cout << " Time to clean up: " << time_3 << " s" << std::endl;
        #else
            std::cout << "Refine bf " << p_bf->Id() << " completed" << std::endl;
        #endif
        }
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::RefineNodes(boost::python::list& pyList)
    {
        // extract the python list to std::set
        std::set<std::size_t> NodeIds;
        typedef boost::python::stl_input_iterator<int> iterator_value_type;
        BOOST_FOREACH(const iterator_value_type::value_type& node_id,
                      std::make_pair(iterator_value_type(pyList), // begin
                      iterator_value_type() ) ) // end
            NodeIds.insert(node_id);

        // cluster the nodes based on the level
        std::map<unsigned int, std::set<std::size_t> > NodesToRefineAtLevel;
        for(std::set<std::size_t>::iterator it = NodeIds.begin(); it != NodeIds.end(); ++it)
        {
            unsigned int level = mBasisFuncs[*it]->Level();
            NodesToRefineAtLevel[level].insert(*it);
        }

        // refine at each level from low to high
        for(std::map<unsigned int, std::set<std::size_t> >::iterator it = NodesToRefineAtLevel.begin(); it != NodesToRefineAtLevel.end(); ++it)
        {
            std::cout << "Level " << it->first << " has " << it->second.size() << " bfs to refine" << std::endl;
            for(std::set<std::size_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
                this->Refine(*it2);
        }
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::RefineWindow(const double& Xi_min, const double& Xi_max,
            const double& Eta_min, const double& Eta_max, const double& Zeta_min, const double& Zeta_max)
    {
        // search and mark all basis functions need to refine on all level (starting from the last level) which support is contained in the refining domain
        double bounding_box_xi_min;
        double bounding_box_xi_max;
        double bounding_box_eta_min;
        double bounding_box_eta_max;
        double bounding_box_zeta_min;
        double bounding_box_zeta_max;
        std::map<int, std::vector<bf_t> > ToBeRefinedBfs;
        for(bf_container_t::iterator it_bf = mBasisFuncs.begin(); it_bf != mBasisFuncs.end(); ++it_bf)
        {
            // get the bounding box (support domain of the basis function)
            (*it_bf)->GetBoundingBox( bounding_box_xi_min,
                                      bounding_box_xi_max,
                                      bounding_box_eta_min,
                                      bounding_box_eta_max,
                                      bounding_box_zeta_min,
                                      bounding_box_zeta_max );

            // check if the bounding box lie in the refined domain
            // Remarks: this can be changed by a refinement indicator (i.e from error estimator)
            if(    bounding_box_xi_min >= Xi_min
                && bounding_box_xi_max <= Xi_max
                && bounding_box_eta_min >= Eta_min
                && bounding_box_eta_max <= Eta_max
                && bounding_box_zeta_min >= Zeta_min
                && bounding_box_zeta_max <= Zeta_max )
            {
                this->Refine((*it_bf)->Id());
            }
        }
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::LinearDependencyRefine(const std::size_t& refine_cycle)
    {
        if(mLastLevel < 1)
            return;

        #ifdef ENABLE_PROFILING
        double start = OpenMPUtils::GetCurrentTime();
        #endif

        // rebuild support domain
        mSupportDomains.clear();
        for(std::size_t level = 1; level <= mLastLevel; ++level)
        {
            domain_t p_domain = GetSupportDomain(level);

            // add the knots to the domain manager
            for(std::size_t next_level = level; next_level <= mLastLevel; ++next_level)
                for(bf_container_t::iterator it_bf = mBasisFuncs.begin(); it_bf != mBasisFuncs.end(); ++it_bf)
                    if((*it_bf)->Level() == next_level)
                        for(DeprecatedHBBasisFunction::cell_iterator it_cell = (*it_bf)->cell_begin(); it_cell != (*it_bf)->cell_end(); ++it_cell)
                            if(TDim == 2)
                            {
                                p_domain->AddXcoord((*it_cell)->LeftValue());
                                p_domain->AddXcoord((*it_cell)->RightValue());
                                p_domain->AddYcoord((*it_cell)->DownValue());
                                p_domain->AddYcoord((*it_cell)->UpValue());
                            }
                            else if(TDim == 3)
                            {
                                p_domain->AddXcoord((*it_cell)->LeftValue());
                                p_domain->AddXcoord((*it_cell)->RightValue());
                                p_domain->AddYcoord((*it_cell)->DownValue());
                                p_domain->AddYcoord((*it_cell)->UpValue());
                                p_domain->AddZcoord((*it_cell)->BelowValue());
                                p_domain->AddZcoord((*it_cell)->AboveValue());
                            }

            // add the cells to the domain manager
            for(std::size_t next_level = level; next_level <= mLastLevel; ++next_level)
                for(bf_container_t::iterator it_bf = mBasisFuncs.begin(); it_bf != mBasisFuncs.end(); ++it_bf)
                    if((*it_bf)->Level() == next_level)
                        for(DeprecatedHBBasisFunction::cell_iterator it_cell = (*it_bf)->cell_begin(); it_cell != (*it_bf)->cell_end(); ++it_cell)
                            if(TDim == 2)
                            {
                                std::vector<double> box = {(*it_cell)->LeftValue(), (*it_cell)->RightValue(), (*it_cell)->DownValue(), (*it_cell)->UpValue()};
                                p_domain->AddCell(box);
                            }
                            else if(TDim == 3)
                            {
                                std::vector<double> box = {(*it_cell)->LeftValue(), (*it_cell)->RightValue(), (*it_cell)->DownValue(), (*it_cell)->UpValue(), (*it_cell)->BelowValue(), (*it_cell)->AboveValue()};
                                p_domain->AddCell(box);
                            }

//            std::cout << "support domain level " << level << *p_domain << std::endl;
        }

        // refine based on the rule that if a bf has support domain contained in next level, it must be refined
        for(std::size_t level = 1; level <= mLastLevel - 1; ++level)
        {
            std::vector<std::size_t> refined_bfs;
            for(bf_container_t::iterator it_bf = mBasisFuncs.begin(); it_bf != mBasisFuncs.end(); ++it_bf)
            {
                // extract the support domain of the next level
                if((*it_bf)->Level() != level)
                    continue;
                domain_t p_domain = GetSupportDomain(level + 1);

                // get the support domain of the bf
                double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;
                (*it_bf)->GetBoundingBox(Xmin, Xmax, Ymin, Ymax, Zmin, Zmax);

                // check if the bf support domain contained in the refined domain managed by the domain manager
                bool is_inside;
                if(TDim == 2)
                {
                    std::vector<double> bounding_box = {Xmin, Xmax, Ymin, Ymax};
                    is_inside = p_domain->IsInside(bounding_box);
                }
                else if(TDim == 3)
                {
                    std::vector<double> bounding_box = {Xmin, Xmax, Ymin, Ymax, Zmin, Zmax};
                    is_inside = p_domain->IsInside(bounding_box);
                }

                if(is_inside)
                    refined_bfs.push_back((*it_bf)->Id());
            }

            if(refined_bfs.size() > 0)
            {
                if((GetEchoLevel() & ECHO_REFINEMENT) == ECHO_REFINEMENT)
                {
                    std::cout << "Additional Bf";
                    for(std::size_t i = 0; i < refined_bfs.size(); ++i)
                        std::cout << " " << refined_bfs[i];
                    std::cout << " of level " << level << " will be refined to maintain the linear independence..." << std::endl;
                }

                for(std::size_t i = 0; i < refined_bfs.size(); ++i)
                    this->Refine(refined_bfs[i]);

                // perform another round to make sure all bfs has support domain in the domain manager of each level
                this->LinearDependencyRefine(refine_cycle + 1);
            }
        }

        #ifdef ENABLE_PROFILING
        std::cout << "LinearDependencyRefine cycle " << refine_cycle << " completed: " << OpenMPUtils::GetCurrentTime() - start << " s" << std::endl;
        #else
        std::cout << "LinearDependencyRefine cycle " << refine_cycle << " completed" << std::endl;
        #endif
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::BuildMesh()
    {
        Vector Crow;
        for(cell_container_t::iterator it_cell = mpCellManager->begin(); it_cell != mpCellManager->end(); ++it_cell)
        {
            (*it_cell)->Reset();
            for(typename CellType::bf_iterator it_bf = (*it_cell)->bf_begin(); it_bf != (*it_cell)->bf_end(); ++it_bf)
            {
                if(TDim == 2)
                    (*it_bf)->ComputeExtractionOperator(*it_cell, Crow, mOrder1, mOrder2);
                else if(TDim == 3)
                    (*it_bf)->ComputeExtractionOperator(*it_cell, Crow, mOrder1, mOrder2, mOrder3);
                (*it_cell)->AddAnchor((*it_bf)->Id(), (*it_bf)->GetControlPoint().W(), Crow);
            }
        }
    }

    /// Validate the patch
    template<int TDim>
    bool DeprecatedHBMesh<TDim>::Validate() const
    {
        // TODO
        return false;
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::BuildBoundaryMesh(DeprecatedHBMesh<TDim>& rMesh, std::string boundary_mesh_type) const
    {
       // TODO
       return;
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::ExportCellTopology(std::string fn, bool cell_numbering) const
    {
        std::ofstream outfile(fn.c_str());
        outfile << "%% Cell topology generated from hierarchical B-Splines mesh, (c) Hoang Giang Bui, 2018\n";
        outfile << "clc\n";
        outfile << "close all\n";
        outfile << "hold on\n";
        outfile << "axis equal\n";

        for(typename cell_container_t::iterator it = mpCellManager->begin(); it != mpCellManager->end(); ++it)
        {
            cell_t p_cell = (*it);
            if(TDim == 2)
            {
//                outfile << "line([" << p_cell->LeftValue() << " " << p_cell->RightValue() << "],[" << p_cell->DownValue() << " " << p_cell->DownValue() << "]);\n";
//                outfile << "line([" << p_cell->LeftValue() << " " << p_cell->RightValue() << "],[" << p_cell->UpValue() << " " << p_cell->UpValue() << "]);\n";
//                outfile << "line([" << p_cell->LeftValue() << " " << p_cell->LeftValue() << "],[" << p_cell->DownValue() << " " << p_cell->UpValue() << "]);\n";
//                outfile << "line([" << p_cell->RightValue() << " " << p_cell->RightValue() << "],[" << p_cell->DownValue() << " " << p_cell->UpValue() << "]);\n";

                outfile << "verts = [" << p_cell->LeftValue() << " " << p_cell->DownValue() << ";"
                                       << p_cell->RightValue() << " " << p_cell->DownValue() << ";"
                                       << p_cell->RightValue() << " " << p_cell->UpValue() << ";"
                                       << p_cell->LeftValue() << " " << p_cell->UpValue() << "];\n";
                outfile << "faces = [1 2 3 4];\n";
                outfile << "patch('Faces',faces,'Vertices',verts,'FaceColor','white');\n";

//                outfile << "rectangle('Position',[" << p_cell->LeftValue() << "," << p_cell->DownValue()
//                        << (p_cell->RightValue() - p_cell->LeftValue()) << ","
//                        << (p_cell->UpValue() - p_cell-DownValue()) << "]);\n";

                if(cell_numbering)
                    outfile << "text(" << 0.5*(p_cell->LeftValue() + p_cell->RightValue()) << ","
                                       << 0.5*(p_cell->DownValue() + p_cell->UpValue()) << ",'"
                                       << p_cell->Id() << "');\n";
            }
            else if(TDim == 3)
            {
                outfile << "verts = [" << p_cell->LeftValue() << " " << p_cell->DownValue() << " " << p_cell->BelowValue() << ";"
                                       << p_cell->RightValue() << " " << p_cell->DownValue() << " " << p_cell->BelowValue() << ";"
                                       << p_cell->RightValue() << " " << p_cell->UpValue() << " " << p_cell->BelowValue() << ";"
                                       << p_cell->LeftValue() << " " << p_cell->UpValue() << " " << p_cell->BelowValue() << ";"
                                       << p_cell->LeftValue() << " " << p_cell->DownValue() << " " << p_cell->AboveValue() << ";"
                                       << p_cell->RightValue() << " " << p_cell->DownValue() << " " << p_cell->AboveValue() << ";"
                                       << p_cell->RightValue() << " " << p_cell->UpValue() << " " << p_cell->AboveValue() << ";"
                                       << p_cell->LeftValue() << " " << p_cell->UpValue() << " " << p_cell->AboveValue() << "];\n";
                outfile << "faces = [1 2 3 4;2 6 7 3;4 3 7 8;1 5 8 4;1 2 6 5;5 6 7 8];\n";
                outfile << "patch('Faces',faces,'Vertices',verts,'FaceColor','white','FaceAlpha',0.5);\n";

                if(cell_numbering)
                    outfile << "text(" << 0.5*(p_cell->LeftValue() + p_cell->RightValue()) << ","
                                       << 0.5*(p_cell->DownValue() + p_cell->UpValue()) << ","
                                       << 0.5*(p_cell->BelowValue() + p_cell->AboveValue()) << ",'"
                                       << p_cell->Id() << "');\n";
            }
        }

        outfile.close();

        std::cout << "Export cell topology to " << fn << " completed" << std::endl;
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::ExportCellGeology(std::string fn)
    {
        // generate cell geology
        std::vector<unsigned int> point_list;
        std::map<unsigned int, double> X_list;
        std::map<unsigned int, double> Y_list;
        std::map<unsigned int, double> Z_list;
        std::map<unsigned int, double> xi_list;
        std::map<unsigned int, double> eta_list;
        std::map<unsigned int, double> zeta_list;
        std::map<unsigned int, unsigned int> cell_list;
        std::map<unsigned int, std::vector<std::vector<unsigned int> > > Connectivities;
        GenerateCellGeology(point_list, X_list, Y_list, Z_list, xi_list, eta_list, zeta_list, cell_list, Connectivities);

        // create file handler
        std::ofstream outfile(fn.c_str());
        outfile << "%% hierarchical B-Splines mesh cell geology, (c) Hoang Giang Bui, 2018\n";
        outfile << "clc\n";
        outfile << "close all\n";
        outfile << "hold on\n";
        outfile << "axis equal\n";
        if(TDim == 2)
        {
//            for(std::set<unordered_pair<unsigned int> >::iterator it = Lines.begin(); it != Lines.end(); ++it)
//                outfile << "line([" << X_list[it->first()] << " " << X_list[it->second()] << "],[" << Y_list[it->first()] << " " << Y_list[it->second()] << "]);\n";

            outfile << "verts = [";
            for(std::size_t i = 0; i < point_list.size(); ++i)
                outfile << X_list[point_list[i]] << " " << Y_list[point_list[i]] << ";\n";
            outfile << "];\n";

            outfile << "faces3 = [";
            for(std::map<unsigned int, std::vector<std::vector<unsigned int> > >::iterator it = Connectivities.begin(); it != Connectivities.end(); ++it) // iterate through cells
                for(std::size_t i = 0; i < it->second.size(); ++i) // iterate through elements in cell
                    if(it->second[i].size() == 3)
                    {
                        for(std::size_t j = 0; j < 3; ++j)
                            outfile << " " << it->second[i][j];
                        outfile << ";\n";
                    }
            outfile << "];\n";
            outfile << "patch('Faces',faces3,'Vertices',verts,'FaceColor','white');\n";

            outfile << "faces4 = [";
            for(std::map<unsigned int, std::vector<std::vector<unsigned int> > >::iterator it = Connectivities.begin(); it != Connectivities.end(); ++it) // iterate through cells
                for(std::size_t i = 0; i < it->second.size(); ++i) // iterate through elements in cell
                    if(it->second[i].size() == 4)
                    {
                        for(std::size_t j = 0; j < 4; ++j)
                            outfile << " " << it->second[i][j];
                        outfile << ";\n";
                    }
            outfile << "];\n";
            outfile << "patch('Faces',faces4,'Vertices',verts,'FaceColor','white');\n";
        }
        else if(TDim == 3)
        {
            outfile << "verts = [";
            for(std::size_t i = 0; i < point_list.size(); ++i)
                outfile << X_list[point_list[i]] << " " << Y_list[point_list[i]] << " " << Z_list[point_list[i]] << ";\n";
            outfile << "];\n";

            outfile << "faces3 = [";
            for(std::map<unsigned int, std::vector<std::vector<unsigned int> > >::iterator it = Connectivities.begin(); it != Connectivities.end(); ++it) // iterate through cells
                for(std::size_t i = 0; i < it->second.size(); ++i) // iterate through elements in cell
                    if(it->second[i].size() == 4)
                    {
                        outfile << it->second[i][0] << " " << it->second[i][1] << " " << it->second[i][2] << ";\n";
                        outfile << it->second[i][0] << " " << it->second[i][1] << " " << it->second[i][3] << ";\n";
                        outfile << it->second[i][0] << " " << it->second[i][2] << " " << it->second[i][3] << ";\n";
                        outfile << it->second[i][1] << " " << it->second[i][2] << " " << it->second[i][3] << ";\n";
                    }
            outfile << "];\n";
            outfile << "patch('Faces',faces3,'Vertices',verts,'FaceColor','white','FaceAlpha',0.0);\n";

            outfile << "faces4 = [";
            for(std::map<unsigned int, std::vector<std::vector<unsigned int> > >::iterator it = Connectivities.begin(); it != Connectivities.end(); ++it) // iterate through cells
                for(std::size_t i = 0; i < it->second.size(); ++i) // iterate through elements in cell
                    if(it->second[i].size() == 8)
                    {
                        outfile << it->second[i][0] << " " << it->second[i][1] << " " << it->second[i][2] << " " << it->second[i][3] << ";\n";
                        outfile << it->second[i][4] << " " << it->second[i][5] << " " << it->second[i][6] << " " << it->second[i][7] << ";\n";
                        outfile << it->second[i][0] << " " << it->second[i][1] << " " << it->second[i][5] << " " << it->second[i][4] << ";\n";
                        outfile << it->second[i][1] << " " << it->second[i][2] << " " << it->second[i][6] << " " << it->second[i][5] << ";\n";
                        outfile << it->second[i][2] << " " << it->second[i][3] << " " << it->second[i][7] << " " << it->second[i][6] << ";\n";
                        outfile << it->second[i][3] << " " << it->second[i][0] << " " << it->second[i][4] << " " << it->second[i][7] << ";\n";
                    }
            outfile << "];\n";
            outfile << "patch('Faces',faces4,'Vertices',verts,'FaceColor','white','FaceAlpha',0.0);\n";
        }

        outfile.close();

        std::cout << "Export cell geology to " << fn << " completed" << std::endl;
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::ExportSupportDomain(std::string fn, double distance)
    {
        std::ofstream outfile(fn.c_str());
        outfile << "%% Support domain topology generated from hierarchical B-Splines mesh, (c) Hoang Giang Bui, 2018\n";
        outfile << "clc\n";
//        outfile << "close all\n";
        outfile << "hold on\n";
        outfile << "axis equal\n";
        outfile.close();

//        std::string color_library[] = {"'r'", "'g'", "'b'", "'y'", "'m'", "'c'", "'w'", "'k'"};

        for(unsigned int level = 1; level <= mLastLevel; ++level)
        {
            double level_distance = level * distance;
            std::stringstream Color;
            double aux = (double)(level - 1) / mLastLevel;
            Color << "[" << 0.5 << "," << 0.5 << "," << aux << "]";
            GetSupportDomain(level)->ExportDomain(fn, Color.str(), level_distance);
        }

        std::cout << "Export support domain to " << fn << " completed" << std::endl;
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::ExportMatlab(std::string fn) const
    {
        std::ofstream outfile(fn.c_str());
        outfile << "%% hierarchical B-Splines mesh information, (c) Hoang Giang Bui, 2018\n";
        outfile << "clc\n";
        outfile << "clear\n";
        outfile << "close all\n";
        outfile << "hold on\n";
        outfile << "axis equal\n\n";

        // export the basis function information
        int cnt = 0;
        std::map<int, int> bfMap; // map from bf Id to the local id
        std::vector<double> LocalKnots1, LocalKnots2, LocalKnots3;
        double min_xi = static_cast<double>(INT_MAX);
        double max_xi = -min_xi;
        double min_eta = min_xi;
        double max_eta = -min_eta;
        double min_zeta = min_xi;
        double max_zeta = -min_zeta;
        for(bf_container_t::iterator it_bf = mBasisFuncs.begin(); it_bf != mBasisFuncs.end(); ++it_bf)
        {
            (*it_bf)->GetLocalKnots(1, LocalKnots1);
            (*it_bf)->GetLocalKnots(2, LocalKnots2);
            if(TDim == 3)
                (*it_bf)->GetLocalKnots(3, LocalKnots3);

            double min_xi_bf = *std::min_element(LocalKnots1.begin(), LocalKnots1.end());
            double max_xi_bf = *std::max_element(LocalKnots1.begin(), LocalKnots1.end());
            if(min_xi_bf < min_xi)
                min_xi = min_xi_bf;
            if(max_xi_bf > max_xi)
                max_xi = max_xi_bf;

            double min_eta_bf = *std::min_element(LocalKnots2.begin(), LocalKnots2.end());
            double max_eta_bf = *std::max_element(LocalKnots2.begin(), LocalKnots2.end());
            if(min_eta_bf < min_eta)
                min_eta = min_eta_bf;
            if(max_eta_bf > max_eta)
                max_eta = max_eta_bf;

            if(TDim == 3)
            {
                double min_zeta_bf = *std::min_element(LocalKnots3.begin(), LocalKnots3.end());
                double max_zeta_bf = *std::max_element(LocalKnots3.begin(), LocalKnots3.end());
                if(min_zeta_bf < min_zeta)
                    min_zeta = min_zeta_bf;
                if(max_zeta_bf > max_zeta)
                    max_zeta = max_zeta_bf;
            }

            ++cnt;
            bfMap[(*it_bf)->Id()] = cnt;

            outfile << "Xi{" << cnt << "} = [";
            for(std::size_t i = 0; i < LocalKnots1.size(); ++i)
                outfile << " " << LocalKnots1[i];
            outfile << "];\n";

            outfile << "Eta{" << cnt << "} = [";
            for(std::size_t i = 0; i < LocalKnots2.size(); ++i)
                outfile << " " << LocalKnots2[i];
            outfile << "];\n";

            if(TDim == 3)
            {
                outfile << "Zeta{" << cnt << "} = [";
                for(std::size_t i = 0; i < LocalKnots3.size(); ++i)
                    outfile << " " << LocalKnots3[i];
                outfile << "];\n";
            }

            outfile << "P(" << cnt << ",:) = [" << (*it_bf)->GetControlPoint().X() << " " << (*it_bf)->GetControlPoint().Y() << " " << (*it_bf)->GetControlPoint().Z() << "];\n";
            outfile << "W(" << cnt << ") = " << (*it_bf)->GetControlPoint().W() << ";\n";
            outfile << "Id(" << cnt << ") = " << (*it_bf)->Id() << ";\n";
            outfile << std::endl;
        }

        // export the cell information
        cnt = 0;
        for(cell_container_t::iterator it_cell = mpCellManager->begin(); it_cell != mpCellManager->end(); ++it_cell)
        {
            ++cnt;

            // write the boundary of the cell
            outfile << "S{" << cnt << "} = [" << (*it_cell)->LeftValue() << " " << (*it_cell)->RightValue() << ";";
            outfile << (*it_cell)->DownValue() << " " << (*it_cell)->UpValue() << "];\n";

            // write the extraction operator
            Matrix C = (*it_cell)->GetExtractionOperator();
            outfile << "C{" << cnt << "} = [";
            for(std::size_t i = 0; i < C.size1(); ++i)
            {
                for(std::size_t j = 0;  j < C.size2(); ++ j)
                    outfile << " " << C(i, j);
                outfile << ";";
            }
            outfile << "];\n";

            // write the supported basis functions
            const std::vector<std::size_t>& bfs = (*it_cell)->GetSupportedAnchors();
            outfile << "N{" << cnt << "} = [";
            for(std::size_t i = 0; i < bfs.size(); ++i)
                outfile << " " << bfs[i];
            outfile << "];\n" << std::endl;
        }

        outfile << std::endl;

        outfile.close();
        std::cout << "Export mesh information to " << fn << " completed" << std::endl;
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::ExportMDPA(std::string fn) const
    {
        std::ofstream outfile(fn.c_str());

        // write header
        outfile << "//KRATOS isogeometric application data file for hierarchical B-Splines\n";
        outfile << "//(c) 2018 Hoang Giang Bui, Ruhr-University Bochum\n";

        IsogeometricMathUtils::timestamp(outfile);

        // write model_part data section
        outfile << "Begin ModelPartData\n";
        outfile << "End ModelPartData\n\n";

        // write properties
        outfile << "Begin Properties 1\n";
        outfile << "End Properties\n\n";

        // write nodes
        outfile << "Begin Nodes\n";
        int BfId = 0;
        #ifdef MDPA_NODE_RENUMBERING
        std::map<int, int> OldToNewBfId; // map from old basis function id to new basis function id
        #endif
        for(typename bf_container_t::iterator it_bf = mBasisFuncs.begin(); it_bf != mBasisFuncs.end(); ++it_bf)
        {
            #ifdef MDPA_NODE_RENUMBERING
            outfile << ++BfId << "\t" << (*it_bf)->GetControlPoint().X() << "\t" << (*it_bf)->GetControlPoint().Y() << "\t" << (*it_bf)->GetControlPoint().Z() << std::endl; // assign new basis function id
            OldToNewBfId[(*it_bf)->Id()] = BfId;
            #else
            outfile << (*it_bf)->Id() << "\t" << (*it_bf)->GetControlPoint().X() << "\t" << (*it_bf)->GetControlPoint().Y() << "\t" << (*it_bf)->GetControlPoint().Z() << std::endl; // assign new basis function id
            #endif
        }
        outfile << "End Nodes\n\n";

        // write elements
        if(TDim == 2)
            outfile << "Begin Elements KinematicLinearGeo2dBezier\n";
        else if(TDim == 3)
            outfile << "Begin Elements KinematicLinearGeo3dBezier\n";
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid Dimension", "")
        #ifdef MDPA_CELL_RENUMBERING
        int ElemId = 0;
        std::map<int, int> OldToNewElemId; // map from old element id to new element id
        #endif
        for(typename cell_container_t::iterator it = mpCellManager->begin(); it != mpCellManager->end(); ++it)
        {
            const std::vector<std::size_t>& supported_bfs = (*it)->GetSupportedAnchors();
            #ifdef MDPA_CELL_RENUMBERING
            outfile << ++ElemId << " 1"; // assign new element Id
            OldToNewElemId[(*it)->Id()] = ElemId;
            #else
            outfile << (*it)->Id() << " 1"; // assign new element Id
            #endif
            for(std::size_t i = 0; i < supported_bfs.size(); ++i)
                #ifdef MDPA_NODE_RENUMBERING
                outfile << " " << OldToNewBfId[supported_bfs[i]];
                #else
                outfile << " " << supported_bfs[i];
                #endif
            outfile << std::endl;
        }
        outfile << "End Elements\n\n";

        // write weights
        outfile << "Begin ElementalData NURBS_WEIGHT\n";
        for(typename cell_container_t::iterator it = mpCellManager->begin(); it != mpCellManager->end(); ++it)
        {
            const std::vector<double>& weights = (*it)->GetAnchorWeights();
            #ifdef MDPA_CELL_RENUMBERING
            outfile << OldToNewElemId[(*it)->Id()] << " [" << weights.size() << "] (";
            #else
            outfile << (*it)->Id() << " [" << weights.size() << "] (";
            #endif
            for(std::size_t i = 0; i < weights.size() - 1; ++i)
                outfile << weights[i] << ",";
            outfile << weights.back() << ")\n";
        }
        outfile << "End ElementalData\n\n";

        // read extraction operator from element in the form of triplet CSR
        // this requires that the Id of the cell must be unique
        std::vector<int> rowPtr;
        std::vector<int> colInd;
        std::vector<double> values;
        std::map<int, std::vector<int> > rowPtrMap;
        std::map<int, std::vector<int> > colIndMap;
        std::map<int, std::vector<double> > valuesMap;
        for(typename cell_container_t::iterator it = mpCellManager->begin(); it != mpCellManager->end(); ++it)
        {
            rowPtr.clear();
            colInd.clear();
            values.clear();
            (*it)->GetExtractionOperator(rowPtr, colInd, values);
            rowPtrMap[(*it)->Id()] = rowPtr;
            colIndMap[(*it)->Id()] = colInd;
            valuesMap[(*it)->Id()] = values;
        }

        // write extraction operator
        outfile << "Begin ElementalData EXTRACTION_OPERATOR_CSR_ROWPTR\n";
        for(std::map<int, std::vector<int> >::iterator it = rowPtrMap.begin(); it != rowPtrMap.end(); ++it)
        {
            #ifdef MDPA_CELL_RENUMBERING
            outfile << OldToNewElemId[(*it).first] << " [" << (*it).second.size() << "] (";
            #else
            outfile << (*it).first << " [" << (*it).second.size() << "] (";
            #endif
            for(std::size_t i = 0; i < (*it).second.size() - 1; ++i)
                outfile << (*it).second[i] << ",";
            outfile << (*it).second.back() << ")\n";
        }
        outfile << "End ElementalData\n\n";

        outfile << "Begin ElementalData EXTRACTION_OPERATOR_CSR_COLIND\n";
        for(std::map<int, std::vector<int> >::iterator it = colIndMap.begin(); it != colIndMap.end(); ++it)
        {
            #ifdef MDPA_CELL_RENUMBERING
            outfile << OldToNewElemId[(*it).first] << " [" << (*it).second.size() << "] (";
            #else
            outfile << (*it).first << " [" << (*it).second.size() << "] (";
            #endif
            for(std::size_t i = 0; i < (*it).second.size() - 1; ++i)
                outfile << (*it).second[i] << ",";
            outfile << (*it).second.back() << ")\n";
        }
        outfile << "End ElementalData\n\n";

        outfile << "Begin ElementalData EXTRACTION_OPERATOR_CSR_VALUES\n";
        for(std::map<int, std::vector<double> >::iterator it = valuesMap.begin(); it != valuesMap.end(); ++it)
        {
            #ifdef MDPA_CELL_RENUMBERING
            outfile << OldToNewElemId[(*it).first] << " [" << (*it).second.size() << "] (";
            #else
            outfile << (*it).first << " [" << (*it).second.size() << "] (";
            #endif
            for(std::size_t i = 0; i < (*it).second.size() - 1; ++i)
                outfile << (*it).second[i] << ",";
            outfile << (*it).second.back() << ")\n";
        }
        outfile << "End ElementalData\n\n";

        // write the degree
        outfile << "Begin ElementalData NURBS_DEGREE_1\n";
        for(typename cell_container_t::iterator it = mpCellManager->begin(); it != mpCellManager->end(); ++it)
            #ifdef MDPA_CELL_RENUMBERING
            outfile << OldToNewElemId[(*it)->Id()] << " " << mOrder1 << std::endl;
            #else
            outfile << (*it)->Id() << " " << mOrder1 << std::endl;
            #endif
        outfile << "End ElementalData\n\n";

        outfile << "Begin ElementalData NURBS_DEGREE_2\n";
        for(typename cell_container_t::iterator it = mpCellManager->begin(); it != mpCellManager->end(); ++it)
            #ifdef MDPA_CELL_RENUMBERING
            outfile << OldToNewElemId[(*it)->Id()] << " " << mOrder2 << std::endl;
            #else
            outfile << (*it)->Id() << " " << mOrder2 << std::endl;
            #endif
        outfile << "End ElementalData\n\n";

        if(TDim == 3)
        {
            outfile << "Begin ElementalData NURBS_DEGREE_3\n";
            for(typename cell_container_t::iterator it = mpCellManager->begin(); it != mpCellManager->end(); ++it)
                #ifdef MDPA_CELL_RENUMBERING
                outfile << OldToNewElemId[(*it)->Id()] << " " << mOrder3 << std::endl;
                #else
                outfile << (*it)->Id() << " " << mOrder3 << std::endl;
                #endif
            outfile << "End ElementalData\n\n";
        }

        // write the division
        outfile << "Begin ElementalData NUM_DIVISION_1\n";
        for(typename cell_container_t::iterator it = mpCellManager->begin(); it != mpCellManager->end(); ++it)
            #ifdef MDPA_CELL_RENUMBERING
            outfile << OldToNewElemId[(*it)->Id()] << " 1" << std::endl;
            #else
            outfile << (*it)->Id() << " 1" << std::endl;
            #endif
        outfile << "End ElementalData\n\n";

        outfile << "Begin ElementalData NUM_DIVISION_2\n";
        for(typename cell_container_t::iterator it = mpCellManager->begin(); it != mpCellManager->end(); ++it)
            #ifdef MDPA_CELL_RENUMBERING
            outfile << OldToNewElemId[(*it)->Id()] << " 1" << std::endl;
            #else
            outfile << (*it)->Id() << " 1" << std::endl;
            #endif
        outfile << "End ElementalData\n\n";

        if(TDim == 3)
        {
            outfile << "Begin ElementalData NUM_DIVISION_3\n";
            for(typename cell_container_t::iterator it = mpCellManager->begin(); it != mpCellManager->end(); ++it)
                #ifdef MDPA_CELL_RENUMBERING
                outfile << OldToNewElemId[(*it)->Id()] << " 1" << std::endl;
                #else
                outfile << (*it)->Id() << " 1" << std::endl;
                #endif
            outfile << "End ElementalData\n\n";
        }

        outfile.close();
        std::cout << "Export MDPA to " << fn << " completed" << std::endl;
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::ExportMDPA2(std::string fn) const
    {
        std::ofstream outfile(fn.c_str());

        // write header
        outfile << "//KRATOS isogeometric application data file for hierarchical B-Splines\n";
        outfile << "//(c) 2018 Hoang Giang Bui, Ruhr-University Bochum\n";

        IsogeometricMathUtils::timestamp(outfile);

        // write model_part data section
        outfile << "Begin ModelPartData\n";
        outfile << "End ModelPartData\n\n";

        // write properties
        outfile << "Begin Properties 1\n";
        outfile << "End Properties\n\n";

        // write nodes
        outfile << "Begin Nodes\n";
        int BfId = 0;
        #ifdef MDPA_NODE_RENUMBERING
        std::map<int, int> OldToNewBfId; // map from old basis function id to new basis function id
        #endif
        for(typename bf_container_t::iterator it_bf = mBasisFuncs.begin(); it_bf != mBasisFuncs.end(); ++it_bf)
        {
            #ifdef MDPA_NODE_RENUMBERING
            outfile << ++BfId << "\t" << (*it_bf)->GetControlPoint().X() << "\t" << (*it_bf)->GetControlPoint().Y() << "\t" << (*it_bf)->GetControlPoint().Z() << std::endl; // assign new basis function id
            OldToNewBfId[(*it_bf)->Id()] = BfId;
            #else
            outfile << (*it_bf)->Id() << "\t" << (*it_bf)->GetControlPoint().X() << "\t" << (*it_bf)->GetControlPoint().Y() << "\t" << (*it_bf)->GetControlPoint().Z() << std::endl; // assign new basis function id
            #endif
        }
        outfile << "End Nodes\n\n";

        // write bezier block
        outfile << "Begin BezierBlock\n";

        outfile << "    Begin IsogeometricBezierData\n";
        outfile << "    // geom_id  number_of_anchors  local_dim  workspace_dim  order_1  order_2  order_3\n";
        std::vector<int> rowPtr;
        std::vector<int> colInd;
        std::vector<double> values;
        for(typename cell_container_t::iterator it = mpCellManager->begin(); it != mpCellManager->end(); ++it)
        {
            outfile << "        "
                    << (*it)->Id() << " " << (*it)->NumberOfAnchors() << " ";

            if(TDim == 2)
                outfile << "2 2 ";
            else if(TDim == 3)
                outfile << "3 3 ";

            outfile << mOrder1 << " " << mOrder2;
            if(TDim == 3)
                outfile << " " << mOrder3 << std::endl;
            else
                outfile << " 0" << std::endl;

            const std::vector<double>& weights = (*it)->GetAnchorWeights();
            outfile << "        [" << weights.size() << "] (";
            for(std::size_t i = 0; i < weights.size() - 1; ++i)
                outfile << weights[i] << ",";
            outfile << weights.back() << ")\n";

            outfile << "        CSR\n";

            rowPtr.clear();
            colInd.clear();
            values.clear();
            (*it)->GetExtractionOperator(rowPtr, colInd, values);

            outfile << "        [" << rowPtr.size() << "] (";
            for(std::size_t i = 0; i < rowPtr.size() - 1; ++i)
                outfile << rowPtr[i] << ",";
            outfile << rowPtr.back() << ")\n";

            outfile << "        [" << colInd.size() << "] (";
            for(std::size_t i = 0; i < colInd.size() - 1; ++i)
                outfile << colInd[i] << ",";
            outfile << colInd.back() << ")\n";

            outfile << "        [" << values.size() << "] (";
            for(std::size_t i = 0; i < values.size() - 1; ++i)
                outfile << values[i] << ",";
            outfile << values.back() << ")\n";

            outfile << std::endl;
        }
        outfile << "    End IsogeometricBezierData\n\n";

        outfile << "    Begin ElementsWithGeometry";
        if(TDim == 2)
            outfile << " KinematicLinearBezier2D\n";
        else if(TDim == 3)
            outfile << " KinematicLinearBezier3D\n";
        for(typename cell_container_t::iterator it = mpCellManager->begin(); it != mpCellManager->end(); ++it)
        {
            const std::vector<std::size_t>& supported_bfs = (*it)->GetSupportedAnchors();
            outfile << "        ";
            #ifdef MDPA_CELL_RENUMBERING
            outfile << ++ElemId << " 1"; // assign new element Id
            #else
            outfile << (*it)->Id() << " 1"; // assign new element Id
            #endif
            outfile << " " << (*it)->Id(); // geometry id
            for(std::size_t i = 0; i < supported_bfs.size(); ++i)
                #ifdef MDPA_NODE_RENUMBERING
                outfile << " " << OldToNewBfId[supported_bfs[i]];
                #else
                outfile << " " << supported_bfs[i];
                #endif
            outfile << std::endl;
        }
        outfile << "    End ElementsWithGeometry\n";

        outfile << "End BezierBlock\n";

        outfile.close();
        std::cout << "Export MDPA to " << fn << " completed" << std::endl;
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::ExportPostMDPA(std::string fn, int NumDivision1, int NumDivision2, int NumDivision3)
    {
        // create a vector of sampling knots in each direction
        std::vector<double> SamplingKnots1;
        std::vector<double> SamplingKnots2;
        std::vector<double> SamplingKnots3;

        int min_xi = (*(mKnots1.begin()))->Value();
        int max_xi = (*(mKnots1.end() - 1))->Value();
        for(std::size_t i = 0; i < NumDivision1 + 1; ++i)
            SamplingKnots1.push_back(min_xi + (double) i / NumDivision1 * (max_xi - min_xi));

        int min_eta = (*(mKnots2.begin()))->Value();
        int max_eta = (*(mKnots2.end() - 1))->Value();
        for(std::size_t i = 0; i < NumDivision2 + 1; ++i)
            SamplingKnots2.push_back(min_eta + (double) i / NumDivision2 * (max_eta - min_eta));

        if(TDim == 3)
        {
            int min_zeta = (*(mKnots3.begin()))->Value();
            int max_zeta = (*(mKnots3.end() - 1))->Value();
            for(std::size_t i = 0; i < NumDivision3 + 1; ++i)
                SamplingKnots3.push_back(min_zeta + (double) i / NumDivision3 * (max_zeta - min_zeta));
        }

        // generate a list of point based on sampling points
        std::vector<int> point_list;
        std::map<int, double> X_list;
        std::map<int, double> Y_list;
        std::map<int, double> Z_list;
        std::map<int, double> xi_list;
        std::map<int, double> eta_list;
        std::map<int, double> zeta_list;
        std::map<int, int> parent_element_list;
        std::vector<std::vector<unsigned int> > connectivities_list;
        if(TDim == 2)
        {
            int cnt = 0;
            Vector W;
            for(std::size_t i = 0; i < SamplingKnots1.size(); ++i)
            {
                for(std::size_t j = 0; j < SamplingKnots2.size(); ++j)
                {
                    // check if this point lie in which cell
                    bool found = false;
                    for(cell_container_t::iterator it_cell = mpCellManager->begin(); it_cell != mpCellManager->end(); ++it_cell)
                    {
                        if((*it_cell)->IsCoverage(SamplingKnots1[i], SamplingKnots2[j]))
                        {
                            found = true;

                            // compute the local coordinates of the sampling point
                            double local_xi = (SamplingKnots1[i] - (*it_cell)->LeftValue()) / ((*it_cell)->RightValue() - (*it_cell)->LeftValue());
                            double local_eta = (SamplingKnots2[j] - (*it_cell)->DownValue()) / ((*it_cell)->UpValue() - (*it_cell)->DownValue());
//                            KRATOS_WATCH(local_xi)
//                            KRATOS_WATCH(local_eta)

                            // compute the Bernstein basis function on the local coordinates
                            Vector B((mOrder1 + 1) * (mOrder2 + 1));
                            for(unsigned int k = 0; k < mOrder1 + 1; ++k)
                            {
                                for(unsigned int l = 0; l < mOrder2 + 1; ++l)
                                {
                                    unsigned int num = k * (mOrder2 + 1) + l;
                                    double B1 = BezierUtils::bernstein2(k, mOrder1, local_xi);
                                    double B2 = BezierUtils::bernstein2(l, mOrder2, local_eta);
                                    B(num) = B1 * B2;
                                }
                            }

                            // get the extraction operator on this cell
                            Matrix C = (*it_cell)->GetExtractionOperator();

                            // compute the B-splines shape function values at the local coordinates
                            Vector N((*it_cell)->NumberOfAnchors());
                            noalias(N) = prod(C, B);

                            // compute the NURBS shape function at the local coordinates
                            (*it_cell)->GetAnchorWeights(W);
                            double Denom = inner_prod(W, N);
                            Vector R((*it_cell)->NumberOfAnchors());
                            for(std::size_t k = 0; k < (*it_cell)->NumberOfAnchors(); ++k)
                                R(k) = W(k) * N(k) / Denom;

                            // get the list of supported basis function of this cell
                            const std::vector<std::size_t>& bfs_id = (*it_cell)->GetSupportedAnchors();

                            // compute the glocal coordinates at the local coordinates
                            double X = 0.0, Y = 0.0, Z = 0.0;
                            for(std::size_t k = 0; k < (*it_cell)->NumberOfAnchors(); ++k)
                            {
                                std::size_t func_id = static_cast<std::size_t>(bfs_id[k]);
                                X += R(k) * mBasisFuncs[func_id]->GetControlPoint().X();
                                Y += R(k) * mBasisFuncs[func_id]->GetControlPoint().Y();
                                Z += R(k) * mBasisFuncs[func_id]->GetControlPoint().Z();
                            }

                            // add point to list
                            point_list.push_back(++cnt);
                            X_list[cnt] = X;
                            Y_list[cnt] = Y;
                            Z_list[cnt] = Z;

                            xi_list[cnt] = local_xi;
                            eta_list[cnt] = local_eta;
                            zeta_list[cnt] = 0.0;

                            parent_element_list[cnt] = (*it_cell)->Id();
                            break;
                        }
                    }
                    if(found == false)
                        KRATOS_THROW_ERROR(std::logic_error, "The sampling points is not detected in the hierarchical mesh", "")
                }
            }

            // generate the element connectivities
            std::vector<unsigned int> connectivities(4);
            for(std::size_t i = 0; i < NumDivision1; ++i)
            {
                for(std::size_t j = 0; j < NumDivision2; ++j)
                {
//                    connectivities[0] = i * (NumDivision2 + 1) + j + 1;
//                    connectivities[1] = connectivities[0] + 1;
//                    connectivities[3] = connectivities[0] + NumDivision2 + 1;
//                    connectivities[2] = connectivities[3] + 1;
//                    connectivities_list.push_back(connectivities);

                    unsigned int Node1 = i * (NumDivision2 + 1) + j + 1;
                    unsigned int Node2 = i * (NumDivision2 + 1) + j + 2;
                    unsigned int Node3 = (i + 1) * (NumDivision2 + 1) + j + 1;
                    unsigned int Node4 = (i + 1) * (NumDivision2 + 1) + j + 2;
                    connectivities[0] = Node1;
                    connectivities[1] = Node2;
                    connectivities[2] = Node4;
                    connectivities[3] = Node3;
                    connectivities_list.push_back(connectivities);
                }
            }
        }
        else if(TDim == 3)
        {
            int cnt = 0;
            Vector W;
            for(std::size_t i = 0; i < SamplingKnots1.size(); ++i)
            {
                for(std::size_t j = 0; j < SamplingKnots2.size(); ++j)
                {
                    for(std::size_t k = 0; k < SamplingKnots3.size(); ++k)
                    {
                        // check if this point lie in which cell
                        bool found = false;
                        for(cell_container_t::iterator it_cell = mpCellManager->begin(); it_cell != mpCellManager->end(); ++it_cell)
                        {
                            if((*it_cell)->IsCoverage(SamplingKnots1[i], SamplingKnots2[j], SamplingKnots3[k]))
                            {
                                found = true;

                                // compute the local coordinates of the sampling point
                                double local_xi = (SamplingKnots1[i] - (*it_cell)->LeftValue()) / ((*it_cell)->RightValue() - (*it_cell)->LeftValue());
                                double local_eta = (SamplingKnots2[j] - (*it_cell)->DownValue()) / ((*it_cell)->UpValue() - (*it_cell)->DownValue());
                                double local_zeta = (SamplingKnots3[k] - (*it_cell)->BelowValue()) / ((*it_cell)->AboveValue() - (*it_cell)->BelowValue());
//                                KRATOS_WATCH(local_xi)
//                                KRATOS_WATCH(local_eta)
//                                KRATOS_WATCH(local_zeta)

                                // compute the Bernstein basis function on the local coordinates
                                Vector B((mOrder1 + 1) * (mOrder2 + 1) * (mOrder3 + 1));
                                for(unsigned int i1 = 0; i1 < mOrder1 + 1; ++i1)
                                    for(unsigned int i2 = 0; i2 < mOrder2 + 1; ++i2)
                                        for(unsigned int i3 = 0; i3 < mOrder3 + 1; ++i3)
                                        {
                                            unsigned int num = (i1 * (mOrder2 + 1) + i2) * (mOrder3 + 1) + i3;
                                            double B1 = BezierUtils::bernstein2(i1, mOrder1, local_xi);
                                            double B2 = BezierUtils::bernstein2(i2, mOrder2, local_eta);
                                            double B3 = BezierUtils::bernstein2(i3, mOrder3, local_zeta);
                                            B(num) = B1 * B2 * B3;
                                        }

                                // get the extraction operator on this cell
                                Matrix C = (*it_cell)->GetExtractionOperator();

                                // compute the B-splines shape function values at the local coordinates
                                Vector N((*it_cell)->NumberOfAnchors());
                                noalias(N) = prod(C, B);

                                // compute the NURBS shape function at the local coordinates
                                (*it_cell)->GetAnchorWeights(W);
                                double Denom = inner_prod(W, N);
                                Vector R((*it_cell)->NumberOfAnchors());
                                for(std::size_t k = 0; k < (*it_cell)->NumberOfAnchors(); ++k)
                                    R(k) = W(k) * N(k) / Denom;

                                // get the list of supported basis function of this cell
                                const std::vector<std::size_t>& bfs_id = (*it_cell)->GetSupportedAnchors();

                                // compute the glocal coordinates at the local coordinates
                                double X = 0.0, Y = 0.0, Z = 0.0;
                                for(std::size_t k = 0; k < (*it_cell)->NumberOfAnchors(); ++k)
                                {
                                    std::size_t func_id = static_cast<std::size_t>(bfs_id[k]);
                                    X += R(k) * mBasisFuncs[func_id]->GetControlPoint().X();
                                    Y += R(k) * mBasisFuncs[func_id]->GetControlPoint().Y();
                                    Z += R(k) * mBasisFuncs[func_id]->GetControlPoint().Z();
                                }

                                // add point to list
                                point_list.push_back(++cnt);
                                X_list[cnt] = X;
                                Y_list[cnt] = Y;
                                Z_list[cnt] = Z;

                                xi_list[cnt] = local_xi;
                                eta_list[cnt] = local_eta;
                                zeta_list[cnt] = 0.0;

                                parent_element_list[cnt] = (*it_cell)->Id();
                                break;
                            }
                        }
                        if(found == false)
                            KRATOS_THROW_ERROR(std::logic_error, "The sampling points is not detected in the hierarchical mesh", "")
                    }
                }
            }

            // generate the element connectivities
            std::vector<unsigned int> connectivities(8);
            for(std::size_t i = 0; i < NumDivision1; ++i)
            {
                for(std::size_t j = 0; j < NumDivision2; ++j)
                {
                    for(std::size_t k = 0; k < NumDivision3; ++k)
                    {
                        unsigned int Node1 = (i * (NumDivision2 + 1) + j) * (NumDivision3 + 1) + k + 1;
                        unsigned int Node2 = (i * (NumDivision2 + 1) + j + 1) * (NumDivision3 + 1) + k + 1;
                        unsigned int Node3 = ((i + 1) * (NumDivision2 + 1) + j) * (NumDivision3 + 1) + k + 1;
                        unsigned int Node4 = ((i + 1) * (NumDivision2 + 1) + j + 1) * (NumDivision3 + 1) + k + 1;
                        unsigned int Node5 = Node1 + 1;
                        unsigned int Node6 = Node2 + 1;
                        unsigned int Node7 = Node3 + 1;
                        unsigned int Node8 = Node4 + 1;

                        connectivities[0] = Node1;
                        connectivities[1] = Node2;
                        connectivities[2] = Node4;
                        connectivities[3] = Node3;
                        connectivities[4] = Node5;
                        connectivities[5] = Node6;
                        connectivities[6] = Node8;
                        connectivities[7] = Node7;

                        connectivities_list.push_back(connectivities);
                    }
                }
            }
        }

        // export to post MDPA
        std::ofstream outfile(fn.c_str());

        // write header
        outfile << "//KRATOS isogeometric application data file for post processing of Bezier-based discretization\n";
        outfile << "//(c) 2018 Hoang Giang Bui, Ruhr-University Bochum\n";

        IsogeometricMathUtils::timestamp(outfile);

        // write model_part data section
        outfile << "Begin ModelPartData\n";
        outfile << "End ModelPartData\n\n";

        // write properties
        outfile << "Begin Properties 1\n";
        outfile << "End Properties\n\n";

        // write nodes
        outfile << "Begin Nodes\n";
        for(std::size_t i = 0; i < point_list.size(); ++i)
            outfile << point_list[i] << " " << X_list[point_list[i]]
                                     << " " << Y_list[point_list[i]]
                                     << " " << Z_list[point_list[i]]
                                     << "\n";
        outfile << "End Nodes\n\n";

        // write elements
        if(TDim == 2)
            outfile << "Begin Elements KinematicLinear2D4N\n";
        else if(TDim == 3)
            outfile << "Begin Elements KinematicLinear3D8N\n";
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid Dimension", "")
        for(std::size_t i = 0; i < connectivities_list.size(); ++i)
        {
            outfile << (i + 1) << " 1";
            for(std::size_t j = 0; j < connectivities_list[i].size(); ++j)
                outfile << " " << connectivities_list[i][j];
            outfile << std::endl;
        }
        outfile << "End Elements\n\n";

        // write nodal data
        outfile << "Begin NodalData LOCAL_COORDINATES\n";
        for(std::size_t i = 0; i < point_list.size(); ++i)
            outfile << point_list[i] << " 0 [3] (" << xi_list[point_list[i]]
                                     << "," << eta_list[point_list[i]]
                                     << "," << zeta_list[point_list[i]] << ")\n";
        outfile << "End NodalData\n\n";

        outfile << "Begin NodalData PARENT_ELEMENT_ID\n";
        for(std::size_t i = 0; i < point_list.size(); ++i)
            outfile << point_list[i] << " 0 " << parent_element_list[point_list[i]] << "\n";
        outfile << "End NodalData\n\n";

        std::cout << "Export post MDPA to " << fn << " completed" << std::endl;
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::ExportCellGeologyAsPostMDPA(std::string fn)
    {
        // generate cell geology
        std::vector<unsigned int> point_list;
        std::map<unsigned int, double> X_list;
        std::map<unsigned int, double> Y_list;
        std::map<unsigned int, double> Z_list;
        std::map<unsigned int, double> xi_list;
        std::map<unsigned int, double> eta_list;
        std::map<unsigned int, double> zeta_list;
        std::map<unsigned int, unsigned int> cell_list;
        std::map<unsigned int, std::vector<std::vector<unsigned int> > > connectivities_list;
        GenerateCellGeology(point_list, X_list, Y_list, Z_list, xi_list, eta_list, zeta_list, cell_list, connectivities_list);

        // export to post MDPA
        std::ofstream outfile(fn.c_str());

        // write header
        outfile << "//KRATOS isogeometric application data file for post processing of Bezier-based discretization\n";
        outfile << "//(c) 2018 Hoang Giang Bui, Ruhr-University Bochum\n";

        IsogeometricMathUtils::timestamp(outfile);

        // write model_part data section
        outfile << "Begin ModelPartData\n";
        outfile << "End ModelPartData\n\n";

        // write properties
        outfile << "Begin Properties 1\n";
        outfile << "End Properties\n\n";

        // write nodes
        outfile << "Begin Nodes\n";
        for(std::size_t i = 0; i < point_list.size(); ++i)
            outfile << point_list[i] << " " << X_list[point_list[i]]
                                     << " " << Y_list[point_list[i]]
                                     << " " << Z_list[point_list[i]]
                                     << "\n";
        outfile << "End Nodes\n\n";

        // write elements
        if(TDim == 2)
        {
            unsigned int ElementId = 0;

            // write triangle elements
            outfile << "Begin Elements KinematicLinear2D3N\n";
            for(std::map<unsigned int, std::vector<std::vector<unsigned int> > >::iterator it = connectivities_list.begin(); it != connectivities_list.end(); ++it) // iterate through cells
                for(std::size_t i = 0; i < it->second.size(); ++i) // iterate through elements in cell
                    if(it->second[i].size() == 3)
                    {
                        outfile << ++ElementId << " 1";
                        for(std::size_t j = 0; j < it->second[i].size(); ++j)
                            outfile << " " << it->second[i][j];
                        outfile << std::endl;
                    }
            outfile << "End Elements\n\n";

            // write quadrilateral elements
            outfile << "Begin Elements KinematicLinear2D4N\n";
            for(std::map<unsigned int, std::vector<std::vector<unsigned int> > >::iterator it = connectivities_list.begin(); it != connectivities_list.end(); ++it) // iterate through cells
                for(std::size_t i = 0; i < it->second.size(); ++i) // iterate through elements in cell
                    if(it->second[i].size() == 4)
                    {
                        outfile << ++ElementId << " 1";
                        for(std::size_t j = 0; j < it->second[i].size(); ++j)
                            outfile << " " << it->second[i][j];
                        outfile << std::endl;
                    }
            outfile << "End Elements\n\n";
        }
        else if(TDim == 3)
        {
            unsigned int ElementId = 0;

            // write tetrahedra elements
            outfile << "Begin Elements KinematicLinear3D4N\n";
            for(std::map<unsigned int, std::vector<std::vector<unsigned int> > >::iterator it = connectivities_list.begin(); it != connectivities_list.end(); ++it) // iterate through cells
                for(std::size_t i = 0; i < it->second.size(); ++i) // iterate through elements in cell
                    if(it->second[i].size() == 4)
                    {
                        outfile << ++ElementId << " 1";
                        for(std::size_t j = 0; j < it->second[i].size(); ++j)
                            outfile << " " << it->second[i][j];
                        outfile << std::endl;
                    }
            outfile << "End Elements\n\n";

            // write hexahedra elements
            outfile << "Begin Elements KinematicLinear3D8N\n";
            for(std::map<unsigned int, std::vector<std::vector<unsigned int> > >::iterator it = connectivities_list.begin(); it != connectivities_list.end(); ++it) // iterate through cells
                for(std::size_t i = 0; i < it->second.size(); ++i) // iterate through elements in cell
                    if(it->second[i].size() == 8)
                    {
                        outfile << ++ElementId << " 1";
                        for(std::size_t j = 0; j < it->second[i].size(); ++j)
                            outfile << " " << it->second[i][j];
                        outfile << std::endl;
                    }
            outfile << "End Elements\n\n";
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid Dimension", "")

        // write nodal data
        outfile << "Begin NodalData LOCAL_COORDINATES\n";
        for(std::size_t i = 0; i < point_list.size(); ++i)
            outfile << point_list[i] << " 0 [3] (" << xi_list[point_list[i]]
                                     << "," << eta_list[point_list[i]]
                                     << "," << zeta_list[point_list[i]] << ")\n";
        outfile << "End NodalData\n\n";

        outfile << "Begin NodalData PARENT_ELEMENT_ID\n";
        for(std::size_t i = 0; i < point_list.size(); ++i)
            outfile << point_list[i] << " 0 " << cell_list[point_list[i]] << "\n";
        outfile << "End NodalData\n\n";

        std::cout << "Export post MDPA to " << fn << " completed" << std::endl;
    }

    template<int TDim>
    void DeprecatedHBMesh<TDim>::GenerateCellGeology(std::vector<unsigned int>& point_list,
                                     std::map<unsigned int, double>& X_list,
                                     std::map<unsigned int, double>& Y_list,
                                     std::map<unsigned int, double>& Z_list,
                                     std::map<unsigned int, double>& xi_list,
                                     std::map<unsigned int, double>& eta_list,
                                     std::map<unsigned int, double>& zeta_list,
                                     std::map<unsigned int, unsigned int>& cell_list,
                                     std::map<unsigned int, std::vector<std::vector<unsigned int> > >& Connectivities)
    {
        // extract all the cell's vertices and lines
        AutoCollapseSpatialBinning Binning(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0e-6);
        std::map<unsigned int, cell_t> MapVertexToCell; // this container map from vertex id to the cell contain it
        std::map<unsigned int, std::vector<unsigned int> > MapCellToNodes; // this container contains the vertex Ids in each cell
        for(cell_container_t::iterator it = mpCellManager->begin(); it != mpCellManager->end(); ++it)
        {
            if(TDim == 2)
            {
                unsigned int V1 = Binning.AddNode((*it)->LeftValue(), (*it)->DownValue(), 0.0);
                unsigned int V2 = Binning.AddNode((*it)->RightValue(), (*it)->DownValue(), 0.0);
                unsigned int V3 = Binning.AddNode((*it)->LeftValue(), (*it)->UpValue(), 0.0);
                unsigned int V4 = Binning.AddNode((*it)->RightValue(), (*it)->UpValue(), 0.0);

                MapVertexToCell[V1] = (*it);
                MapVertexToCell[V2] = (*it);
                MapVertexToCell[V3] = (*it);
                MapVertexToCell[V4] = (*it);
                MapCellToNodes[(*it)->Id()].push_back(V1);
                MapCellToNodes[(*it)->Id()].push_back(V2);
                MapCellToNodes[(*it)->Id()].push_back(V4);
                MapCellToNodes[(*it)->Id()].push_back(V3);
            }
            else if(TDim == 3)
            {
                unsigned int V1 = Binning.AddNode((*it)->LeftValue(), (*it)->DownValue(), (*it)->BelowValue());
                unsigned int V2 = Binning.AddNode((*it)->RightValue(), (*it)->DownValue(), (*it)->BelowValue());
                unsigned int V3 = Binning.AddNode((*it)->LeftValue(), (*it)->UpValue(), (*it)->BelowValue());
                unsigned int V4 = Binning.AddNode((*it)->RightValue(), (*it)->UpValue(), (*it)->BelowValue());
                unsigned int V5 = Binning.AddNode((*it)->LeftValue(), (*it)->DownValue(), (*it)->AboveValue());
                unsigned int V6 = Binning.AddNode((*it)->RightValue(), (*it)->DownValue(), (*it)->AboveValue());
                unsigned int V7 = Binning.AddNode((*it)->LeftValue(), (*it)->UpValue(), (*it)->AboveValue());
                unsigned int V8 = Binning.AddNode((*it)->RightValue(), (*it)->UpValue(), (*it)->AboveValue());

                MapVertexToCell[V1] = (*it);
                MapVertexToCell[V2] = (*it);
                MapVertexToCell[V3] = (*it);
                MapVertexToCell[V4] = (*it);
                MapVertexToCell[V5] = (*it);
                MapVertexToCell[V6] = (*it);
                MapVertexToCell[V7] = (*it);
                MapVertexToCell[V8] = (*it);
                MapCellToNodes[(*it)->Id()].push_back(V1);
                MapCellToNodes[(*it)->Id()].push_back(V2);
                MapCellToNodes[(*it)->Id()].push_back(V4);
                MapCellToNodes[(*it)->Id()].push_back(V3);
                MapCellToNodes[(*it)->Id()].push_back(V5);
                MapCellToNodes[(*it)->Id()].push_back(V6);
                MapCellToNodes[(*it)->Id()].push_back(V8);
                MapCellToNodes[(*it)->Id()].push_back(V7);
            }
        }

        // compute the global coordinates of the vertex
        for(unsigned int i = 0; i < Binning.NumberOfNodes(); ++i)
        {
            unsigned int Id = i + 1;
            cell_t p_cell = MapVertexToCell[Id];

            // compute the local coordinates of the vertex
            double xi, eta, zeta, local_xi, local_eta, local_zeta;
            xi = Binning.GetX(Id);
            eta = Binning.GetY(Id);
            if(TDim == 2)
                zeta = 0.0;
            else if(TDim == 3)
                zeta = Binning.GetZ(Id);

            local_xi = (xi - p_cell->LeftValue()) / (p_cell->RightValue() - p_cell->LeftValue());
            local_eta = (eta - p_cell->DownValue()) / (p_cell->UpValue() - p_cell->DownValue());
            if(TDim == 2)
                local_zeta = 0.0;
            else if(TDim == 3)
                local_zeta = (zeta - p_cell->BelowValue()) / (p_cell->AboveValue() - p_cell->BelowValue());

            // compute the Bernstein basis function on the local coordinates
            Vector B;
            if(TDim == 2)
            {
                B.resize((mOrder1 + 1) * (mOrder2 + 1));
                for(unsigned int k = 0; k < mOrder1 + 1; ++k)
                    for(unsigned int l = 0; l < mOrder2 + 1; ++l)
                    {
                        unsigned int num = k * (mOrder2 + 1) + l;
                        double B1 = BezierUtils::bernstein2(k, mOrder1, local_xi);
                        double B2 = BezierUtils::bernstein2(l, mOrder2, local_eta);
                        B(num) = B1 * B2;
                    }
            }
            else if(TDim == 3)
            {
                B.resize((mOrder1 + 1) * (mOrder2 + 1) * (mOrder3 + 1));
                for(unsigned int k = 0; k < mOrder1 + 1; ++k)
                    for(unsigned int l = 0; l < mOrder2 + 1; ++l)
                        for(unsigned int m = 0; m < mOrder3 + 1; ++m)
                        {
                            unsigned int num = (k * (mOrder2 + 1) + l) * (mOrder3 + 1) + m;
                            double B1 = BezierUtils::bernstein2(k, mOrder1, local_xi);
                            double B2 = BezierUtils::bernstein2(l, mOrder2, local_eta);
                            double B3 = BezierUtils::bernstein2(m, mOrder3, local_zeta);
                            B(num) = B1 * B2 * B3;
                        }
            }

            // get the extraction operator on this cell
            Matrix C = p_cell->GetExtractionOperator();

            // compute the B-splines shape function values at the local coordinates
            Vector N(p_cell->NumberOfAnchors());
            noalias(N) = prod(C, B);

            // compute the NURBS shape function at the local coordinates
            Vector W;
            p_cell->GetAnchorWeights(W);
            double Denom = inner_prod(W, N);
            Vector R(p_cell->NumberOfAnchors());
            for(std::size_t k = 0; k < p_cell->NumberOfAnchors(); ++k)
                R(k) = W(k) * N(k) / Denom;

            // get the list of supported basis function of this cell
            const std::vector<std::size_t>& bfs_id = p_cell->GetSupportedAnchors();

            // compute the glocal coordinates at the local coordinates
            double X = 0.0, Y = 0.0, Z = 0.0;
            for(std::size_t k = 0; k < p_cell->NumberOfAnchors(); ++k)
            {
                std::size_t func_id = static_cast<std::size_t>(bfs_id[k]);
                X += R(k) * mBasisFuncs[func_id]->GetControlPoint().X();
                Y += R(k) * mBasisFuncs[func_id]->GetControlPoint().Y();
                Z += R(k) * mBasisFuncs[func_id]->GetControlPoint().Z();
            }

            // add point to list
            X_list[Id] = X;
            Y_list[Id] = Y;
            Z_list[Id] = Z;
            point_list.push_back(Id);
            xi_list[Id] = local_xi;
            eta_list[Id] = local_eta;
            zeta_list[Id] = local_zeta;
            cell_list[Id] = p_cell->Id();
        }

        // compute the middle nodes in cell and generate the internal mesh for each cell
        for(cell_container_t::iterator it = mpCellManager->begin(); it != mpCellManager->end(); ++it)
        {
            // we check each edge in each cell, if the edge/face contain any vertex, that vertex will be includes in the middle nodes. These cells will be filled with triangular mesh.

            // global Node Ids of cell
            std::vector<unsigned int>& CellNodes = MapCellToNodes[(*it)->Id()];

            // detect middle nodes
            std::set<unsigned int> MiddleNodes;
            if(TDim == 2)
            {
                // extract edges
                typedef std::pair<unsigned int, unsigned int> edge_t;
                std::vector<edge_t> Edges;

                Edges.push_back(edge_t(CellNodes[0], CellNodes[1]));
                Edges.push_back(edge_t(CellNodes[1], CellNodes[2]));
                Edges.push_back(edge_t(CellNodes[2], CellNodes[3]));
                Edges.push_back(edge_t(CellNodes[3], CellNodes[0]));

                // detect middle nodes on edges
                for(std::size_t i = 0; i < Edges.size(); ++i)
                {
                    unsigned int n1 = Edges[i].first;
                    unsigned int n2 = Edges[i].second;

                    // check if the line contain any vertex, then these vertices are middle nodes
                    double xi1  = std::min(Binning.GetX(n1), Binning.GetX(n2));
                    double eta1 = std::min(Binning.GetY(n1), Binning.GetY(n2));
                    double xi2  = std::max(Binning.GetX(n1), Binning.GetX(n2));
                    double eta2 = std::max(Binning.GetY(n1), Binning.GetY(n2));
                    for(unsigned int i = 0; i < Binning.NumberOfNodes(); ++i)
                    {
                        unsigned int Id = i + 1;
                        double xi = Binning.GetX(Id);
                        double eta = Binning.GetY(Id);
                        bool is_inside = false;
                        is_inside = is_inside || (xi1<xi && xi<xi2) && (eta1==eta && eta==eta2); // 1st direction
                        is_inside = is_inside || (xi1==xi && xi==xi2) && (eta1<eta && eta<eta2); // 2nd direction
                        if(is_inside)
                            MiddleNodes.insert(Id);
                    }
                }
            }
            else if(TDim == 3)
            {
                // extract faces
                typedef std::pair<unsigned int, unsigned int> edge_t;
                typedef std::pair<edge_t, edge_t> face_t;
                std::vector<face_t> Faces;

                Faces.push_back(face_t(edge_t(CellNodes[0], CellNodes[1]), edge_t(CellNodes[2], CellNodes[3])));
                Faces.push_back(face_t(edge_t(CellNodes[4], CellNodes[5]), edge_t(CellNodes[6], CellNodes[7])));
                Faces.push_back(face_t(edge_t(CellNodes[0], CellNodes[1]), edge_t(CellNodes[5], CellNodes[4])));
                Faces.push_back(face_t(edge_t(CellNodes[1], CellNodes[2]), edge_t(CellNodes[6], CellNodes[5])));
                Faces.push_back(face_t(edge_t(CellNodes[2], CellNodes[3]), edge_t(CellNodes[7], CellNodes[6])));
                Faces.push_back(face_t(edge_t(CellNodes[3], CellNodes[0]), edge_t(CellNodes[4], CellNodes[7])));

                // detect middle nodes on faces
                for(std::size_t i = 0; i < Faces.size(); ++i)
                {
                    unsigned int n1 = Faces[i].first.first;
                    unsigned int n2 = Faces[i].first.second;
                    unsigned int n3 = Faces[i].second.first;
                    unsigned int n4 = Faces[i].second.second;

                    // check if the face contain any vertex, then these vertices are middle nodes
                    double xi_min  = std::min(std::min(Binning.GetX(n1), Binning.GetX(n2)), std::min(Binning.GetX(n3), Binning.GetX(n4)));
                    double xi_max  = std::max(std::max(Binning.GetX(n1), Binning.GetX(n2)), std::max(Binning.GetX(n3), Binning.GetX(n4)));
                    double eta_min  = std::min(std::min(Binning.GetY(n1), Binning.GetY(n2)), std::min(Binning.GetY(n3), Binning.GetY(n4)));
                    double eta_max  = std::max(std::max(Binning.GetY(n1), Binning.GetY(n2)), std::max(Binning.GetY(n3), Binning.GetY(n4)));
                    double zeta_min  = std::min(std::min(Binning.GetZ(n1), Binning.GetZ(n2)), std::min(Binning.GetZ(n3), Binning.GetZ(n4)));
                    double zeta_max  = std::max(std::max(Binning.GetZ(n1), Binning.GetZ(n2)), std::max(Binning.GetZ(n3), Binning.GetZ(n4)));
                    for(unsigned int i = 0; i < Binning.NumberOfNodes(); ++i)
                    {
                        unsigned int Id = i + 1;
                        double xi = Binning.GetX(Id);
                        double eta = Binning.GetY(Id);
                        double zeta = Binning.GetZ(Id);
                        bool is_inside = false;
                        if(xi_min == xi_max)
                        {
                            is_inside = is_inside || xi == xi_min &&
                                (   ((eta_min<eta && eta<eta_max) && (zeta_min<zeta && zeta<zeta_max))
                                 || ((eta_min==eta || eta==eta_max) && (zeta_min<zeta && zeta<zeta_max))
                                 || ((eta_min<eta && eta<eta_max) && (zeta_min==zeta || zeta==zeta_max))   );
                        }
                        else if(eta_min == eta_max)
                        {
                            is_inside = is_inside || eta == eta_min &&
                                (   ((xi_min<xi && xi<xi_max) && (zeta_min<zeta && zeta<zeta_max))
                                 || ((xi_min==xi || xi==xi_max) && (zeta_min<zeta && zeta<zeta_max))
                                 || ((xi_min<xi && xi<xi_max) && (zeta_min==zeta || zeta==zeta_max))   );
                        }
                        else if(zeta_min == zeta_max)
                        {
                            is_inside = is_inside || zeta == zeta_min &&
                                (   ((xi_min<xi && xi<xi_max) && (eta_min<eta && eta<eta_max))
                                 || ((xi_min==xi || xi==xi_max) && (eta_min<eta && eta<eta_max))
                                 || ((xi_min<xi && xi<xi_max) && (eta_min==eta || eta==eta_max))   );
                        }
                        else
                            KRATOS_THROW_ERROR(std::logic_error, "Something wrong with the surface, it must align with 1 of the 3 base surfaces", "")

                        if(is_inside)
                            MiddleNodes.insert(Id);
                    }
                }
            }

            if(MiddleNodes.size() == 0)
            {
                // the cell has no middle nodes, now generate the connectivities
                std::vector<unsigned int> Conn(CellNodes.begin(), CellNodes.end());
                Connectivities[(*it)->Id()].push_back(Conn);
            }
            else
            {
                // insert the middle vertices to the cell
                CellNodes.insert(CellNodes.end(), MiddleNodes.begin(), MiddleNodes.end());

                // generate triangles in case of 2D
                if(TDim == 2)
                {
                    // prepare the list of points
                    std::vector<double> XYlist;
                    for(std::size_t i = 0; i < CellNodes.size(); ++i)
                    {
                        XYlist.push_back(Binning.GetX(CellNodes[i]));
                        XYlist.push_back(Binning.GetY(CellNodes[i]));
                    }

                    // generate the triangulation
                    TriangulationUtils util;
                    std::vector<std::vector<unsigned int> > tris;
                    util.ComputeDelaunayTriangulation(XYlist, tris);

                    // map the triangle connectivities to node Id
                    std::vector<std::vector<unsigned int> > Conns(tris.size());
                    for(std::size_t i = 0; i < tris.size(); ++i)
                        for(std::size_t j = 0; j < tris[i].size(); ++j)
                            Conns[i].push_back(CellNodes[tris[i][j]]);

                    Connectivities[(*it)->Id()] = Conns;
                }
                // generate tetrahedrons in case of 3D
                else if(TDim == 3)
                {
                    #if !defined(ISOGEOMETRIC_USE_TETGEN)
                    KRATOS_THROW_ERROR(std::runtime_error, __FUNCTION__, "requires Tetgen to generate tetrahedrons for cells")
                    #endif

                    tetgenio in, out;

                    // initialize tetgen I/O
                    in.initialize();
                    out.initialize();

                    // All indices start from 1.
                    in.firstnumber = 1;

                    // fill in points
                    in.numberofpoints = CellNodes.size();
                    in.pointlist = new REAL[in.numberofpoints * 3];
                    for(std::size_t i = 0; i < CellNodes.size(); ++i)
                    {
                        in.pointlist[3*i    ] = Binning.GetX(CellNodes[i]);
                        in.pointlist[3*i + 1] = Binning.GetY(CellNodes[i]);
                        in.pointlist[3*i + 2] = Binning.GetZ(CellNodes[i]);
                    }

                    // tetrahedralize
                    std::string settings("Q"); // silent all output
                    tetrahedralize((char*)settings.c_str(), &in, &out);

                    // extract the tetrahedrons
                    unsigned int cnt = 0;
                    std::vector<std::vector<unsigned int> > Conns(out.numberoftetrahedra);
                    for(std::size_t i = 0; i < out.numberoftetrahedra; ++i)
                        for(std::size_t j = 0; j < 4; ++j)
                            Conns[i].push_back(CellNodes[out.tetrahedronlist[cnt++] - 1]);

                    Connectivities[(*it)->Id()] = Conns;
                }
            }
        }
    }

    /**
     * template instantiation
     */
    template class DeprecatedHBMesh<0>;
    template class DeprecatedHBMesh<1>;
    template class DeprecatedHBMesh<2>;
    template class DeprecatedHBMesh<3>;

} // end namespace Kratos

#undef MDPA_NODE_RENUMBERING
#undef MDPA_CELL_RENUMBERING
#undef ENABLE_PROFILING
#undef DEBUG_REFINE

