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
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
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

        virtual Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const;

        virtual Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
                                PropertiesType::Pointer pProperties) const;

        IntegrationMethod GetIntegrationMethod() const;

        /**
         * Calculates the local system contributions for this contact element
         */
        virtual void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                   VectorType& rRightHandSideVector,
                                   ProcessInfo& rCurrentProcessInfo);

        virtual void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                     ProcessInfo& rCurrentProcessInfo);

        virtual void EquationIdVector( EquationIdVectorType& rResult,
                               ProcessInfo& rCurrentProcessInfo);

        virtual void GetDofList( DofsVectorType& ConditionalDofList,
                         ProcessInfo& CurrentProcessInfo);

        virtual void Initialize(const ProcessInfo& rCurrentProcessInfo);

        /// Turn back information as a string.
        virtual std::string Info() const
        {
            return "DummyIsogeometricElement";
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "DummyIsogeometricElement #" << Id();
        }

        /// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const
        {
            mpIsogeometricGeometry->PrintData(rOStream);
            rOStream << " (IntegrationMethod: " << mThisIntegrationMethod << ")" << std::endl;
        }

    private:

        IsogeometricGeometryType::Pointer mpIsogeometricGeometry;

        IntegrationMethod mThisIntegrationMethod;

        friend class Serializer;

        virtual void save ( Serializer& rSerializer ) const
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Element )
        }

        virtual void load ( Serializer& rSerializer )
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Element )
        }

        void CalculateAll( MatrixType& rLeftHandSideMatrix,
                           VectorType& rRightHandSideVector,
                           ProcessInfo& rCurrentProcessInfo,
                           bool CalculateStiffnessMatrixFlag,
                           bool CalculateResidualVectorFlag);

}; // Class DummyIsogeometricElement

}  // namespace Kratos.


#endif // KRATOS_DUMMY_ISOGEOMETRIC_ELEMENT_H_INCLUDED defined

