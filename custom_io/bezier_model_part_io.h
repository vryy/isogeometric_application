//
//   Project Name:        KratosIsogeometricApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 19 Oct 2015 $
//
//

#if !defined(KRATOS_BEZIER_MODEL_PART_IO_H_INCLUDED )
#define  KRATOS_BEZIER_MODEL_PART_IO_H_INCLUDED

// System includes
#include <string>
#include <fstream>
#include <set>

// External includes

// Project includes
#include "includes/model_part_io.h"

namespace Kratos
{

class BezierInfo : public IndexedObject
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BezierInfo);
    std::size_t n;
    int local_space_dim;
    int global_space_dim;
    int p1, p2, p3;
    Vector weights;
    int mat_type; // 0: full, 1: MCSR, 2: CSR
    Matrix C; // extraction operator

    BezierInfo() : IndexedObject(0) {}

    BezierInfo(std::size_t NewId) : IndexedObject(NewId) {}

    void Print(std::ostream& rOStream) const
    {
        rOStream << "id: " << Id() << std::endl;
        rOStream << "n: " << n << std::endl;
        rOStream << "local_space_dim: " << local_space_dim << std::endl;
        rOStream << "global_space_dim: " << global_space_dim << std::endl;
        rOStream << "p1: " << p1 << std::endl;
        rOStream << "p2: " << p2 << std::endl;
        rOStream << "p3: " << p3 << std::endl;
        rOStream << "weights: " << weights << std::endl;
        rOStream << "C: " << C << std::endl;
    }

private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override
    {
    }

    void load(Serializer& rSerializer) override
    {
    }
};

/// output stream function
inline std::ostream & operator <<(std::ostream& rOStream, const BezierInfo& rThis)
{
    rThis.Print(rOStream);
    return rOStream;
}

class KRATOS_API(ISOGEOMETRIC_APPLICATION) BezierModelPartIO : public ModelPartIO<>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BezierModelPartIO);

    typedef ModelPartIO<> BaseType;

    typedef typename BaseType::NodeType NodeType;

    typedef typename BaseType::MeshType MeshType;

    typedef typename BaseType::NodesContainerType NodesContainerType;

    typedef typename BaseType::PropertiesContainerType PropertiesContainerType;

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    typedef typename BaseType::ConditionsContainerType ConditionsContainerType;

    typedef typename BaseType::ConnectivitiesContainerType ConnectivitiesContainerType;

    typedef typename BaseType::SizeType SizeType;

    typedef PointerVectorSet<BezierInfo, IndexedObject> BezierInfoContainerType;

    /// Default Constructor with  filenames.
    BezierModelPartIO(std::string const& Filename,
#ifdef SD_APP_FORWARD_COMPATIBILITY
                      const Flags Options = BaseIO::READ
#else
                      const Flags Options = BaseIO::READ | BaseIO::NOT_IGNORE_VARIABLES_ERROR
#endif
                     );

    /// Destructor.
    ~BezierModelPartIO() override;

    /// Read the data and initialize the model part
    void ReadModelPart(ModelPart & rThisModelPart) override;

private:

    BezierInfoContainerType::Pointer mpBezierInfoContainer;

    void ReadBezierBlock(ModelPart & rThisModelPart);

    void ReadIsogeometricBezierDataBlock(BezierInfoContainerType& rThisBezierInfo);

    void ReadElementsWithGeometryBlock(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, BezierInfoContainerType& rGeometryInfo, ElementsContainerType& rThisElements);

    void ReadConditionsWithGeometryBlock(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, BezierInfoContainerType& rGeometryInfo, ConditionsContainerType& rThisConditions);

};

}

#endif
