/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 19 Oct 2015 $
//   Revision:            $Revision: 1.0 $
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

class KRATOS_API(ISOGEOMETRIC_APPLICATION) BezierModelPartIO : public ModelPartIO
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BezierModelPartIO);

    typedef ModelPartIO BaseType;

    typedef BaseType::NodeType NodeType;

    typedef BaseType::MeshType MeshType;

    typedef BaseType::NodesContainerType NodesContainerType;

    typedef BaseType::PropertiesContainerType PropertiesContainerType;

    typedef BaseType::ElementsContainerType ElementsContainerType;

    typedef BaseType::ConditionsContainerType ConditionsContainerType;

    typedef BaseType::ConnectivitiesContainerType ConnectivitiesContainerType;

    typedef BaseType::SizeType SizeType;

    typedef PointerVectorSet<BezierInfo, IndexedObject> BezierInfoContainerType;

    /// Default Constructor with  filenames.
    BezierModelPartIO(std::string const& Filename,
#ifdef SD_APP_FORWARD_COMPATIBILITY
                      const Flags Options = IO::READ
#else
                      const Flags Options = IO::READ | IO::NOT_IGNORE_VARIABLES_ERROR
#endif
                     );

    /// Destructor.
    virtual ~BezierModelPartIO();

    /// Read the data and initialize the model part
    virtual void ReadModelPart(ModelPart & rThisModelPart);

private:

    BezierInfoContainerType::Pointer mpBezierInfoContainer;

    void ReadBezierBlock(ModelPart & rThisModelPart);

    void ReadIsogeometricBezierDataBlock(BezierInfoContainerType& rThisBezierInfo);

    void ReadElementsWithGeometryBlock(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, BezierInfoContainerType& rGeometryInfo, ElementsContainerType& rThisElements);

    void ReadConditionsWithGeometryBlock(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, BezierInfoContainerType& rGeometryInfo, ConditionsContainerType& rThisConditions);

};

}

#endif
