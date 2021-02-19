//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 7 Dec 17 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_DUMMY_ISOGEOMETRIC_CONDITION_H_INCLUDED )
#define  KRATOS_DUMMY_ISOGEOMETRIC_CONDITION_H_INCLUDED


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
class DummyIsogeometricCondition : public Condition
{
    public:
        // Counted pointer of DummyIsogeometricCondition
        KRATOS_CLASS_POINTER_DEFINITION(DummyIsogeometricCondition);

        typedef GeometryData::IntegrationMethod IntegrationMethod;

        typedef IsogeometricGeometry<GeometryType::PointType> IsogeometricGeometryType;

        /**
         * Default constructor.
         */
        DummyIsogeometricCondition();
        DummyIsogeometricCondition( IndexType NewId, GeometryType::Pointer pGeometry);
        DummyIsogeometricCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        /**
         * Destructor.
         */
        virtual ~DummyIsogeometricCondition();

        /**
         * Operations.
         */

        Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const final;

        Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
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
            return "DummyIsogeometricCondition";
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const final
        {
            rOStream << "DummyIsogeometricCondition #" << Id();
        }

        /// Print object's data.
        void PrintData(std::ostream& rOStream) const final
        {
            mpIsogeometricGeometry->PrintData(rOStream);
            rOStream << " (IntegrationMethod: " << mThisIntegrationMethod << ")" << std::endl;
        }

    private:

        IsogeometricGeometryType::Pointer mpIsogeometricGeometry;

        IntegrationMethod mThisIntegrationMethod;

        friend class Serializer;

        void save ( Serializer& rSerializer ) const final
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Condition )
        }

        void load ( Serializer& rSerializer ) final
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Condition )
        }

        void CalculateAll( MatrixType& rLeftHandSideMatrix,
                           VectorType& rRightHandSideVector,
                           const ProcessInfo& rCurrentProcessInfo,
                           bool CalculateStiffnessMatrixFlag,
                           bool CalculateResidualVectorFlag);

}; // Class DummyIsogeometricCondition

}  // namespace Kratos.


#endif // KRATOS_DUMMY_ISOGEOMETRIC_CONDITION_H_INCLUDED defined

