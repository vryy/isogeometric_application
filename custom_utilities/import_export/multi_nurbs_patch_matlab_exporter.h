//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_NURBS_PATCH_MATLAB_EXPORTER_H_INCLUDED)
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_NURBS_PATCH_MATLAB_EXPORTER_H_INCLUDED

// System includes
#include <vector>
#include <fstream>
#include <iomanip>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/nurbs/bsplines_fespace.h"
#include "custom_utilities/import_export/multipatch_exporter.h"

namespace Kratos
{

class MultiNURBSPatchMatlabExporterHelper
{
public:
    template<int TDim>
    static void WriteMatlabControlPoints(std::ostream& rOStream, typename Patch<TDim>::Pointer pPatch, const std::string& var_name)
    {
        KRATOS_THROW_ERROR(std::logic_error, "WriteMatlabControlPoints is not implemented for dimension", TDim)
    }
};


/**
Export NURBS patch/multipatch to Matlab to visualize with NURBS toolbox by M. Spink
 */
template<int TDim>
class MultiNURBSPatchMatlabExporterWriter : public MultiPatchExporter<TDim>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiNURBSPatchMatlabExporterWriter);

    /// Type definition
    typedef MultiPatchExporter<TDim> BaseType;
    typedef typename BaseType::knot_container_t knot_container_t;
    typedef typename knot_container_t::knot_t knot_t;

    /// Default constructor
    MultiNURBSPatchMatlabExporterWriter() : BaseType() {}

    /// Destructor
    virtual ~MultiNURBSPatchMatlabExporterWriter() {}

    /// Export a single patch
    virtual void Export(typename Patch<TDim>::Pointer pPatch, std::ostream& rOStream) const
    {
        rOStream << std::setprecision(BaseType::Accuracy());
        this->ExportMatlab(rOStream, pPatch, pPatch->Name());
    }

    /// Export a multipatch
    virtual void Export(typename MultiPatch<TDim>::Pointer pMultiPatch, std::ostream& rOStream) const
    {
        rOStream << std::setprecision(BaseType::Accuracy());

        for (typename MultiPatch<TDim>::PatchContainerType::ptr_const_iterator it = pMultiPatch->Patches().ptr_begin(); it != pMultiPatch->Patches().ptr_end(); ++it)
        {
            std::stringstream patch_name;
            patch_name << "patch" << (*it)->Id();
            this->ExportMatlab(rOStream, *it, patch_name.str());
        }
    }

private:

    void ExportMatlab(std::ostream& rOStream, typename Patch<TDim>::Pointer pPatch, const std::string& patch_name) const
    {
        if (pPatch->pFESpace()->Type() != BSplinesFESpace<TDim>::StaticType())
            KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "does not support non-NURBS patch")

        typename BSplinesFESpace<TDim>::Pointer pFESpace = boost::dynamic_pointer_cast<BSplinesFESpace<TDim> >(pPatch->pFESpace());
        if (pFESpace == NULL)
            KRATOS_THROW_ERROR(std::runtime_error, "The cast to BSplinesFESpace is failed.", "")

        if (TDim == 1)
        {
            rOStream << "knots = [";
            for (std::size_t i = 0; i < pFESpace->KnotVector(0).size(); ++i)
                rOStream << " " << pFESpace->KnotVector(0)[i];
            rOStream << "];\n";
        }
        else
        {
            rOStream << "knots = {};\n";
            for (std::size_t dim = 0; dim < TDim; ++dim)
            {
                rOStream << "knots{" << dim+1 << "} = [";
                for (std::size_t i = 0; i < pFESpace->KnotVector(dim).size(); ++i)
                    rOStream << " " << pFESpace->KnotVector(dim)[i];
                rOStream << "];\n";
            }
        }

        rOStream << "coefs = zeros(4";
        for (std::size_t dim = 0; dim < TDim; ++dim)
            rOStream << "," << pFESpace->Number(dim);
        rOStream << ");\n";

        MultiNURBSPatchMatlabExporterHelper::WriteMatlabControlPoints<TDim>(rOStream, pPatch, std::string("coefs"));

        rOStream <<  patch_name << " = nrbmak(coefs,knots);\n";
        // rOStream << "xlabel('x');\n";
        // rOStream << "ylabel('y');\n";
        // rOStream << "zlabel('z');\n";

        rOStream << patch_name << "_number = [";
        std::vector<std::size_t> func_indices = pPatch->pFESpace()->FunctionIndices();
        for (std::size_t i = 0; i < func_indices.size(); ++i)
            rOStream << " " << CONVERT_INDEX_IGA_TO_KRATOS(func_indices[i]);
        rOStream << "];\n";

        rOStream << std::endl;
    }

}; // end class MultiNURBSPatchMatlabExporterWriter

template<>
void MultiNURBSPatchMatlabExporterHelper::WriteMatlabControlPoints<1>(std::ostream& rOStream, Patch<1>::Pointer pPatch, const std::string& var_name)
{
    typedef Patch<1>::ControlPointType ControlPointType;
    typename StructuredControlGrid<1, ControlPointType>::ConstPointer pControlPointGrid
        = boost::dynamic_pointer_cast<const StructuredControlGrid<1, ControlPointType> >(pPatch->pControlPointGridFunction()->pControlGrid());
    if (pControlPointGrid == NULL)
        KRATOS_THROW_ERROR(std::runtime_error, "The cast to StructuredControlGrid is failed.", "")

    for (std::size_t nu = 0; nu < pControlPointGrid->Size(0); ++nu)
    {
        rOStream << var_name << "(:," << nu+1 << ") = [";
        for (std::size_t dim = 0; dim < 3; ++dim)
            rOStream << " " << pControlPointGrid->GetValue(nu)[dim];
        rOStream << " " << pControlPointGrid->GetValue(nu)[3] << "];\n";
    }
}

template<>
void MultiNURBSPatchMatlabExporterHelper::WriteMatlabControlPoints<2>(std::ostream& rOStream, Patch<2>::Pointer pPatch, const std::string& var_name)
{
    typedef Patch<2>::ControlPointType ControlPointType;
    typename StructuredControlGrid<2, ControlPointType>::ConstPointer pControlPointGrid
        = boost::dynamic_pointer_cast<const StructuredControlGrid<2, ControlPointType> >(pPatch->pControlPointGridFunction()->pControlGrid());
    if (pControlPointGrid == NULL)
        KRATOS_THROW_ERROR(std::runtime_error, "The cast to StructuredControlGrid is failed.", "")

    for (std::size_t nv = 0; nv < pControlPointGrid->Size(1); ++nv)
    {
        for (std::size_t nu = 0; nu < pControlPointGrid->Size(0); ++nu)
        {
            rOStream << var_name << "(:," << nu+1 << "," << nv+1 << ") = [";
            for (std::size_t dim = 0; dim < 3; ++dim)
                rOStream << " " << pControlPointGrid->GetValue(nu, nv)[dim];
            rOStream << " " << pControlPointGrid->GetValue(nu, nv)[3] << "];\n";
        }
    }
}

template<>
void MultiNURBSPatchMatlabExporterHelper::WriteMatlabControlPoints<3>(std::ostream& rOStream, Patch<3>::Pointer pPatch, const std::string& var_name)
{
    typedef Patch<3>::ControlPointType ControlPointType;
    typename StructuredControlGrid<3, ControlPointType>::ConstPointer pControlPointGrid
        = boost::dynamic_pointer_cast<const StructuredControlGrid<3, ControlPointType> >(pPatch->pControlPointGridFunction()->pControlGrid());
    if (pControlPointGrid == NULL)
        KRATOS_THROW_ERROR(std::runtime_error, "The cast to StructuredControlGrid is failed.", "")

    for (std::size_t nw = 0; nw < pControlPointGrid->Size(2); ++nw)
    {
        for (std::size_t nv = 0; nv < pControlPointGrid->Size(1); ++nv)
        {
            for (std::size_t nu = 0; nu < pControlPointGrid->Size(0); ++nu)
            {
                rOStream << var_name << "(:," << nu+1 << "," << nv+1 << "," << nw+1 << ") = [";
                for (std::size_t dim = 0; dim < 3; ++dim)
                    rOStream << " " << pControlPointGrid->GetValue(nu, nv, nw)[dim];
                rOStream << " " << pControlPointGrid->GetValue(nu, nv, nw)[3] << "];\n";
            }
        }
    }
}

class MultiNURBSPatchMatlabExporter
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MultiNURBSPatchMatlabExporter);

    template<int TDim>
    static void Export(typename Patch<TDim>::Pointer pPatch, const std::string& filename)
    {
        std::ofstream outfile;
        outfile.open(filename, std::ios::out);

        MultiNURBSPatchMatlabExporterWriter<TDim> dummy;
        dummy.Export(pPatch, outfile);

        outfile.close();
        std::cout << pPatch->Type() << " " << pPatch->Id() << " is exported to " << filename << " successfully" << std::endl;
    }

    template<int TDim>
    static void Export(typename MultiPatch<TDim>::Pointer pMultiPatch, const std::string& filename)
    {
        std::ofstream outfile;
        outfile.open(filename, std::ios::out);

        MultiNURBSPatchMatlabExporterWriter<TDim> dummy;
        dummy.Export(pMultiPatch, outfile);

        outfile.close();
        std::cout << "Multipatch is exported to " << filename << " successfully" << std::endl;
    }

    /// Information
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MultiNURBSPatchMatlabExporter";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
    }
};


/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const MultiNURBSPatchMatlabExporter& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}



} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTI_NURBS_PATCH_MATLAB_EXPORTER_H_INCLUDED defined

