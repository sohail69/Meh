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
    iter nx, ny, nz;    //Number of nodes in each direction
    iter nxe, nye, nze; //Number of elements in each direction
    num dx, dy, dz;     //the uniform sizes in each direction
  public:

}
