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
    iter nx[3], nxe[3]; //The size of the blockmesh
    num dx[3];     //The resolution of the blockmesh
  public:
    


}
