//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Sep 2014 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_AUTO_COLLAPSE_SPATIAL_BINNING_H_INCLUDED )
#define  KRATOS_AUTO_COLLAPSE_SPATIAL_BINNING_H_INCLUDED

// System includes
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <map>

// External includes 

// Project includes
#include "includes/define.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

class SpatialPoint
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SpatialPoint);
    SpatialPoint(std::size_t id, double x, double y, double z) : mId(id), mX(x), mY(y), mZ(z) {}
    ~SpatialPoint() {}
    void SetId(std::size_t id) {mId = id;}
    std::size_t GetId() const {return mId;}
    double GetX() const {return mX;}
    double GetY() const {return mY;}
    double GetZ() const {return mZ;}
private:
    double mX, mY, mZ;
    std::size_t mId;
}; // class SpatialPoint


/*** Detail class definition.
 * This class defines a key in space. Used for spatial binning algorithm.
 */
class SpatialKey
{
    public:
        KRATOS_CLASS_POINTER_DEFINITION(SpatialKey);
        
        SpatialKey(int ix, int iy, int iz) : x(ix), y(iy), z(iz) {}
        ~SpatialKey() {}
        bool operator<(const SpatialKey& rOther) const
        {
            if(x == rOther.x)
            {
                if(y == rOther.y)
                {
                    return z < rOther.z;
                }
                else
                    return y < rOther.y;
            }
            else
                return x < rOther.x;
        }
        int kx() const {return x;}
        int ky() const {return y;}
        int kz() const {return z;}
    private:
        int x, y, z;
};


/// Short class definition.
/*** Detail class definition.
 * This utility class supports for spatial binning with auto collapsing functionality. This class used SpatialPoint as point data and will collapse the conincident point.
 */
class AutoCollapseSpatialBinning
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;    
    
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(AutoCollapseSpatialBinning);
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AutoCollapseSpatialBinning(double X0, double Y0, double Z0, double Dx, double Dy, double Dz, double tol)
    : mX0(X0), mY0(Y0), mZ0(Z0), mDx(Dx), mDy(Dy), mDz(Dz), mTol(tol)
    {
        mLastNode = 0;
    }

    /// Destructor.
    virtual ~AutoCollapseSpatialBinning()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    
    /**
     * This function add node to the spatial binning and return the id of this node in the bin
     */
    unsigned int AddNode(double X, double Y, double Z)
    {
        // find the cell containing point
        int ix = (int) floor((X - mX0) / mDx);
        int iy = (int) floor((Y - mY0) / mDy);
        int iz = (int) floor((Z - mZ0) / mDz);
        
        // check if the spatial key exist
        SpatialKey key(ix, iy, iz);
    //    std::cout << "ix = " << ix << ", iy = " << iy << ", iz = " << iz << std::endl;
        std::map<SpatialKey, std::vector<int> >::iterator it = mBin.find(key);
        if(it != mBin.end())
        {
            // check if node already exist in the cell
    //        std::cout << "check if node already exist in the cell" << std::endl;
            for(int i = 0; i < it->second.size(); ++i)
            {
                double xi = mPointList[it->second[i] - 1]->GetX();
                double yi = mPointList[it->second[i] - 1]->GetY();
                double zi = mPointList[it->second[i] - 1]->GetZ();
                double d = sqrt(pow(X - xi, 2) + pow(Y - yi, 2) + pow(Z - zi, 2));
                if(d < mTol)
                {
    //                std::cout << "node exist in cell, return its id" << std::endl;
                    return mPointList[it->second[i] - 1]->GetId();
                }
            }
            // node does not exist in cell, insert the node into spatial bin
            mPointList.push_back(boost::shared_ptr<SpatialPoint>(new SpatialPoint(++mLastNode, X, Y, Z)));
            mBin[key].push_back(mLastNode);
    //        std::cout << "node does not exist in cell, insert the node into spatial bin, mLastNode = " << mLastNode << ", mBin[key].size() = " << mBin[key].size() << std::endl;
            return mLastNode;
        }
        else
        {
            // insert the node into spatial bin
            mPointList.push_back(boost::shared_ptr<SpatialPoint>(new SpatialPoint(++mLastNode, X, Y, Z)));
            mBin[key].push_back(mLastNode);
    //        std::cout << "insert the node into spatial bin, mLastNode = " << mLastNode << ", mBin[key].size() = " << mBin[key].size() << std::endl;
            return mLastNode;
        }
    }

    std::size_t NumberOfNodes() const {return mPointList.size();}
    double GetX(unsigned int id) const {return mPointList[id - 1]->GetX();}
    double GetY(unsigned int id) const {return mPointList[id - 1]->GetY();}
    double GetZ(unsigned int id) const {return mPointList[id - 1]->GetZ();}

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "AutoCollapseSpatialBinning";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "AutoCollapseSpatialBinning";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {}

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{
    double mX0, mY0, mZ0, mDx, mDy, mDz;
    unsigned int mLastNode;
    double mTol;

    std::map<SpatialKey, std::vector<int> > mBin;
    std::vector<boost::shared_ptr<SpatialPoint> > mPointList;
    
    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    
    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{
    
    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    AutoCollapseSpatialBinning& operator=(AutoCollapseSpatialBinning const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    AutoCollapseSpatialBinning(AutoCollapseSpatialBinning const& rOther)
    {
    }

    ///@}

}; // Class AutoCollapseSpatialBinning

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
//inline std::istream& operator >>(std::istream& rIStream, AutoCollapseSpatialBinning& rThis)
//{
//    return rIStream;
//}

///// output stream function
//inline std::ostream& operator <<(std::ostream& rOStream,
//        const AutoCollapseSpatialBinning& rThis)
//{
//    rThis.PrintInfo(rOStream);
//    rOStream << std::endl;
//    rThis.PrintData(rOStream);

//    return rOStream;
//}
///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_AUTO_COLLAPSE_SPATIAL_BINNING_H_INCLUDED defined

