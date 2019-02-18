//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Feb 2019 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ISOGEOMETRIC_APPLICATION_IGA_DEFINE_PYTHON_H_INCLUDED )
#define  KRATOS_ISOGEOMETRIC_APPLICATION_IGA_DEFINE_PYTHON_H_INCLUDED

namespace Kratos
{

namespace Python
{

/**
 * Helper class to get the reference to the patch pointer
 */
template<class TPatchType>
TPatchType& GetReference(typename TPatchType::Pointer& dummy)
{
    return *dummy;
}

/**
 * Helper function
 */
template<typename TType>
std::size_t Isogeometric_GetId(TType& rDummy)
{
    return rDummy.Id();
}

/**
 * Helper function
 */
template<typename TType>
void Isogeometric_SetId(TType& rDummy, std::size_t Id)
{
    rDummy.SetId(Id);
}

/**
 * Helper function
 */
template<typename TType>
void Isogeometric_DoNotSetId(TType& rDummy, std::size_t Id)
{
    // DO NOTHING
}

/**
 * Helper function
 */
template<typename TType>
std::size_t Isogeometric_GetEquationId(TType& rDummy)
{
    return rDummy.EquationId();
}

/**
 * Helper function
 */
template<typename TType>
void Isogeometric_SetEquationId(TType& rDummy, std::size_t EquationId)
{
    rDummy.SetEquationId(EquationId);
}

/**
 * Helper function
 */
template<typename TFESpaceType>
boost::python::list FESpace_ExtractBoundaryBfsByFlag(TFESpaceType& rDummy, std::size_t boundary_id)
{
    typedef typename TFESpaceType::bf_t bf_t;

    std::vector<bf_t> bf_list = rDummy.ExtractBoundaryBfsByFlag(boundary_id);

    boost::python::list Output;
    for (std::size_t i = 0; i < bf_list.size(); ++i)
        Output.append(bf_list[i]);

    return Output;
}

/**
 * Helper function
 */
template<typename TFESpaceType>
typename TFESpaceType::bf_t FESpace_GetItem(TFESpaceType& rDummy, std::size_t i)
{
    return rDummy[i];
}

} // namespace Python

} // namespace Kratos.

#endif // KRATOS_ISOGEOMETRIC_APPLICATION_IGA_DEFINE_PYTHON_H_INCLUDED defined

