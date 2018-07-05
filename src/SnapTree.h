// SnapTree.h

#include "Image.h"
#include "signal.hpp"

// define how many 2D images to be loaded each time
#define MAX_IMAGES_STREAM 16

// 16-bit data bits shift method:
// 0: >> 4
// 1: >> 2
//
#define METHOD_BITS_SHIFT 0

// meta
// Folder: Y/Y_X
// File: Y_X_Z.tif

// cube
class Cube
{
public:
    Cube();
    ~Cube();

public:
    int offset_D; // z

    string fileName; // 000000_000000_000000.tif
    string filePath; // ourdir/RESXXxXXxXX/000000/000000_000000/000000_000000_000000.tif
    uint32 depth; // 256
};

// folder
class YXFolder
{
public:
    YXFolder();
    ~YXFolder();

public:
    int offset_V; // y
    int offset_H; // x

    uint16 lengthFileName; // 25 len("000000_000000_000000.tif") + 1
    uint16 lengthDirName; // 21 len("000000/000000_000000") + 1

    string dirName; // 000000/000000_000000
    string xDirPath; // ourdir/RESXXxXXxXX/000000
    string yDirPath; // ourdir/RESXXxXXxXX/000000/000000_000000

    uint32 height, width; // 256x256

    map<int,Cube> cubes;
};

// layer
class Layer
{
public:
    Layer();
    ~Layer();

public:
    uint16 rows, cols, ncubes; // floor(dim_V/height)+1, floor(dim_H/width)+1
    float vs_x, vs_y, vs_z; // voxel sizes
    uint32 dim_V, dim_H, dim_D; // dimensions y, x, z

    string layerName; // outdir/RESXXxXXxXX

    map<string, YXFolder> yxfolders; // <dirName, YXFolder>
};

// tree
class Tree
{
public:
    Tree();
    ~Tree();

public:
    int init(string outdir, long dimx, long dimy, long dimz, long dimc, int bytes, int resolutions);

public:
    float org_V, org_H, org_D; // offsets (0, 0, 0)
    axis reference_V, reference_H, reference_D; // vertical, horizonal, depth
    float mdata_version; // 2

    uint32 color, bytesPerVoxel; //
    uint32 cubex, cubey, cubez;

    map<int, Layer> layers;
};

//
class SnapTree
{
public:
    SnapTree(string inputdir, string outputdir, int scales, int genMetaInfo=0, int sx=0, int sy=0, int sz=0, int startz=-1, int endz=-1);
    ~SnapTree();

public:
    int init();
    uint8* load(long zs, long ze);
    int reformat();
    int index();
    void resume(int startz, int endz);

public:
    string srcdir, dstdir;
    int resolutions;
    uint32 width, height, depth; // 3D image stacks
    uint32 zstart, zend; // split for parallel reformatting

    set<string> input2DTIFFs;
    uint32 block_width, block_height, block_depth;
    uint16 datatype, outDatatype;
    uint32 color;
    uint8 *ubuffer;
    int nbits;
    string resumeConfigFile;

    bool inputImageDimensions;
    bool split;

    Tree meta;
    int genMetaInfoOnly;
};
