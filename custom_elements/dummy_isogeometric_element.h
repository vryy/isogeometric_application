//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Jan 21 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_DUMMY_ISOGEOMETRIC_ELEMENT_H_INCLUDED )
#define  KRATOS_DUMMY_ISOGEOMETRIC_ELEMENT_H_INCLUDED

// External includes

// Project includes
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "custom_geometries/isogeometric_geometry.h"

namespace Kratos
{
/**
 * This condition does nothing. But it gives the isogeometric geometry for other condition that contains it.
 */
class DummyIsogeometricElement : public Element
{
public:
    // Counted pointer of DummyIsogeometricElement
    KRATOS_CLASS_POINTER_DEFINITION(DummyIsogeometricElement);

    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef IsogeometricGeometry<GeometryType::PointType> IsogeometricGeometryType;

    /**
     * Default constructor.
     */
    DummyIsogeometricElement();
    DummyIsogeometricElement( IndexType NewId, GeometryType::Pointer pGeometry);
    DummyIsogeometricElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /**
     * Destructor.
     */
    virtual ~DummyIsogeometricElement();

    /**
     * Operations.
     */

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const final;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const final;

    IntegrationMethod GetIntegrationMethod() const final;

    /**
     * Calculates the local system contributions for this contact element
     */
    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const ProcessInfo& rCurrentProcessInfo) final;

    void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                 const ProcessInfo& rCurrentProcessInfo) final;

    void EquationIdVector( EquationIdVectorType& rResult,
                           const ProcessInfo& rCurrentProcessInfo) const final;

    void GetDofList( DofsVectorType& ConditionalDofList,
                     const ProcessInfo& CurrentProcessInfo) const final;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) final;

    /// Turn back information as a string.
    std::string Info() const final
    {
        return "DummyIsogeometricElement";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << "DummyIsogeometricElement #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        mpIsogeometricGeometry->PrintData(rOStream);
        rOStream << " (IntegrationMethod: " << static_cast<int>(mThisIntegrationMethod) << ")" << std::endl;
    }

private:

    IsogeometricGeometryType::Pointer mpIsogeometricGeometry;

    IntegrationMethod mThisIntegrationMethod;

    friend class Serializer;

    void save ( Serializer& rSerializer ) const final
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Element )
    }

    void load ( Serializer& rSerializer ) final
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Element )
    }

    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       const ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag,
                       bool CalculateResidualVectorFlag);

}; // Class DummyIsogeometricElement

}  // namespace Kratos.

#endif // KRATOS_DUMMY_ISOGEOMETRIC_ELEMENT_H_INCLUDED defined
