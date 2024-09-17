#include <vector>
#include<cmath>

//
// produces a 3-D blockmesh
// does not store coords etc..
// keeps a record of local partition
// and calculates on the fly
//
template<typename num, typename iter>
class blockMesh{
  private:
    num dx[3];          //The resolution of the blockmesh
    iter nx[3], nxe[3]; //The size of the blockmesh
    iter nStart, nEnd;  //Local start and end nodes
    iter neStart, neEnd;//Local start and end element
  public:
    blockMesh(){

    }

}
