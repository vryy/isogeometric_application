//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 25 Nov 2025 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_WRAPPER_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_WRAPPER_H_INCLUDED

#include "custom_utilities/multipatch.h"

namespace Kratos
{

/**
 * Wrapper class of MultiPatch. This can be used to wrap and serialize all types of multipatch.
 * Upon loading, the concrete type of the original multipatch will be returned.
 */
class MultiPatchWrapper
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultiPatchWrapper);

    /// Empty constructor, for deserialization
    MultiPatchWrapper()
    {}

    /// Constructor with multipatch, for serialization
    MultiPatchWrapper(BaseMultiPatch::Pointer pMultiPatch)
    : mpMultiPatch(pMultiPatch)
    {}

    /// Create an empty instance of the multipatch wrapper
    static MultiPatchWrapper::Pointer Create()
    {
        return MultiPatchWrapper::Pointer(new MultiPatchWrapper());
    }

    /// Create a new instance of the multipatch wrapper
    static MultiPatchWrapper::Pointer Create(BaseMultiPatch::Pointer pMultiPatch)
    {
        return MultiPatchWrapper::Pointer(new MultiPatchWrapper(pMultiPatch));
    }

    /// Get the underlying multipatch pointer
    BaseMultiPatch::Pointer Get() const
    {
        return mpMultiPatch;
    }

    /// Create an empty multipatch pointer based on template data type identifier
    static BaseMultiPatch::Pointer CreateEmptyMultiPatch(int dim, const std::string& lc_str, const std::string& c_str, const std::string& d_str)
    {
        typedef typename PatchSelector<1>::RealPatch::DataType DoubleType;
        typedef typename PatchSelector<1>::ComplexPatch::DataType ComplexType;

        constexpr const char* double_str = DataTypeToString<DoubleType>::Get();
        constexpr const char* complex_str = DataTypeToString<ComplexType>::Get();

        typedef PatchSelector<1>::RealMultiPatch MultiPatch1DType;
        typedef PatchSelector<2>::RealMultiPatch MultiPatch2DType;
        typedef PatchSelector<3>::RealMultiPatch MultiPatch3DType;

        typedef PatchSelector<1>::ComplexMultiPatch ComplexMultiPatch1DType;
        typedef PatchSelector<2>::ComplexMultiPatch ComplexMultiPatch2DType;
        typedef PatchSelector<3>::ComplexMultiPatch ComplexMultiPatch3DType;

        typedef PatchSelector<1>::GComplexMultiPatch GComplexMultiPatch1DType;
        typedef PatchSelector<2>::GComplexMultiPatch GComplexMultiPatch2DType;
        typedef PatchSelector<3>::GComplexMultiPatch GComplexMultiPatch3DType;

        BaseMultiPatch::Pointer pMultiPatch = nullptr;
        if      ((dim == 1) && (lc_str == double_str) && (c_str == double_str) && (d_str == double_str))    pMultiPatch = MultiPatch1DType::Create();
        else if ((dim == 2) && (lc_str == double_str) && (c_str == double_str) && (d_str == double_str))    pMultiPatch = MultiPatch2DType::Create();
        else if ((dim == 3) && (lc_str == double_str) && (c_str == double_str) && (d_str == double_str))    pMultiPatch = MultiPatch3DType::Create();
        else if ((dim == 1) && (lc_str == double_str) && (c_str == double_str) && (d_str == complex_str))   pMultiPatch = ComplexMultiPatch1DType::Create();
        else if ((dim == 2) && (lc_str == double_str) && (c_str == double_str) && (d_str == complex_str))   pMultiPatch = ComplexMultiPatch2DType::Create();
        else if ((dim == 3) && (lc_str == double_str) && (c_str == double_str) && (d_str == complex_str))   pMultiPatch = ComplexMultiPatch3DType::Create();
        else if ((dim == 1) && (lc_str == double_str) && (c_str == complex_str) && (d_str == complex_str))  pMultiPatch = GComplexMultiPatch1DType::Create();
        else if ((dim == 2) && (lc_str == double_str) && (c_str == complex_str) && (d_str == complex_str))  pMultiPatch = GComplexMultiPatch2DType::Create();
        else if ((dim == 3) && (lc_str == double_str) && (c_str == complex_str) && (d_str == complex_str))  pMultiPatch = GComplexMultiPatch3DType::Create();

        return pMultiPatch;
    }

private:

    BaseMultiPatch::Pointer mpMultiPatch = nullptr;

    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        int dim = 0;
        std::string lc_str = "", c_str = "", d_str = "";
        if (mpMultiPatch != nullptr)
        {
            dim = mpMultiPatch->Dimension();
            lc_str = mpMultiPatch->LocalCoordinateTypeStr();
            c_str = mpMultiPatch->CoordinateTypeStr();
            d_str = mpMultiPatch->DataTypeStr();
        }

        rSerializer.save("dim", dim);
        rSerializer.save("lc_str", lc_str);
        rSerializer.save("c_str", c_str);
        rSerializer.save("d_str", d_str);

        if (mpMultiPatch != nullptr)
            rSerializer.save("mpatch", *mpMultiPatch);
    }

    void load(Serializer& rSerializer)
    {
        int dim = 0;
        std::string lc_str = "", c_str = "", d_str = "";
        rSerializer.load("dim", dim);
        rSerializer.load("lc_str", lc_str);
        rSerializer.load("c_str", c_str);
        rSerializer.load("d_str", d_str);

        mpMultiPatch = CreateEmptyMultiPatch(dim, lc_str, c_str, d_str);
        if (mpMultiPatch != nullptr) rSerializer.load("mpatch", *mpMultiPatch);
    }
    ///@}
};

} // end namespace Kratos

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_MULTIPATCH_WRAPPER_H_INCLUDED
