// SnapTree.cpp

#include "SnapTree.h"

// cube
Cube::Cube()
{

}

Cube::~Cube()
{

}

YXFolder::YXFolder()
{
    lengthFileName = 25;
    lengthDirName = 21;
}

YXFolder::~YXFolder()
{

}

Layer::Layer()
{
    vs_x = 1;
    vs_y = 1;
    vs_z = 1;
}

Layer::~Layer()
{

}

Tree::Tree()
{
    org_V = 0;
    org_H = 0;
    org_D = 0;

    reference_V = vertical; // y
    reference_H = horizontal; // x
    reference_D = depth; // z

    mdata_version = 2;

    cubex = 256;
    cubey = 256;
    cubez = 256;
}

Tree::~Tree()
{

}

int Tree::init(string outdir, long dimx, long dimy, long dimz, long dimc, int bytes, int resolutions)
{
    //
    color = dimc;
    bytesPerVoxel = bytes; // 1 for 8-bit; 2 for 16-bit

    //
    for(int res_i=0; res_i<resolutions; res_i++)
    {
        //
        Layer layer;

        layer.dim_V = (uint32)(dimy/(uint32)pow(2,res_i));
        layer.dim_H = (uint32)(dimx/(uint32)pow(2,res_i));
        layer.dim_D = (uint32)(dimz/(uint32)pow(2,res_i));

        layer.vs_x = pow(2,res_i);
        layer.vs_y = pow(2,res_i);
        layer.vs_z = pow(2,res_i);

        layer.rows = (layer.dim_V%cubey==0)?floor(layer.dim_V/cubey):floor(layer.dim_V/cubey)+1;
        layer.cols = (layer.dim_H%cubex==0)?floor(layer.dim_H/cubex):floor(layer.dim_H/cubex)+1;
        layer.ncubes = (layer.dim_D%cubez==0)?floor(layer.dim_D/cubez):floor(layer.dim_D/cubez)+1;

        stringstream filepath;
        filepath <<outdir<<"/RES"<<layer.dim_V<<"x"<<layer.dim_H<<"x"<<layer.dim_D;

        layer.layerName = filepath.str();

        //
        for(int row = 0; row < layer.rows; row++)
        {
            //
            stringstream x_pos;
            x_pos.width(6);
            x_pos.fill('0');
            x_pos << row*cubey*(int)pow(2.0, res_i) * 10;

            stringstream V_DIR_path;
            V_DIR_path << layer.layerName << "/" << x_pos.str();

            uint32 length;

            //
            if(row==layer.rows-1)
            {
                if(layer.dim_V%cubey==0)
                {
                    length = cubey;
                }
                else
                {
                    length = layer.dim_V%cubey;
                }
            }
            else
            {
                length = cubey;
            }

            for(int col = 0; col < layer.cols; col++)
            {
                //
                YXFolder yxfolder;

                yxfolder.offset_V = row*cubey;
                yxfolder.offset_H = col*cubex;

                yxfolder.height = length;

                //
                stringstream y_pos;
                y_pos.width(6);
                y_pos.fill('0');
                y_pos << yxfolder.offset_H*(int)pow(2.0,res_i) * 10;

                stringstream H_DIR_path, dir_name;
                H_DIR_path << V_DIR_path.str() << "/" << x_pos.str() << "_" << y_pos.str();
                dir_name << x_pos.str() << "/" << x_pos.str() << "_" << y_pos.str();

                //
                yxfolder.xDirPath = V_DIR_path.str();
                yxfolder.yDirPath = H_DIR_path.str();
                yxfolder.dirName = dir_name.str();

                if(col==layer.cols-1)
                {
                    if(layer.dim_H%cubex==0)
                    {
                        yxfolder.width = cubex;
                    }
                    else
                    {
                        yxfolder.width = layer.dim_H%cubex;
                    }
                }
                else
                {
                    yxfolder.width = cubex;
                }

                //
                for(int z = 0; z<layer.ncubes; z++)
                {
                    //
                    Cube cube;

                    cube.offset_D = z*cubez;

                    stringstream z_pos;
                    z_pos.width(6);
                    z_pos.fill('0');
                    z_pos << cube.offset_D*(int)pow(2.0,res_i) * 10;

                    cube.fileName = x_pos.str() + "_" + y_pos.str() + "_" + z_pos.str() + ".tif";
                    cube.filePath = yxfolder.yDirPath + "/" + cube.fileName;

                    if(z==layer.ncubes-1)
                    {
                        if(layer.dim_D%cubez==0)
                        {
                            cube.depth = cubez;
                        }
                        else
                        {
                            cube.depth = layer.dim_D%cubez;
                        }
                    }
                    else
                    {
                        cube.depth = cubez;
                    }

                    //
                    yxfolder.cubes.insert(make_pair(cube.offset_D, cube));
                }

                //
                layer.yxfolders.insert(make_pair(yxfolder.dirName, yxfolder));

            } // col
        }// row

        //
        layers.insert(make_pair(res_i, layer));

    } // res_i

    //
    return 0;
}

//
SnapTree::SnapTree(string inputdir, string outputdir, int scales, int genMetaInfo, int sx, int sy, int sz, int startz, int endz)
{
    // default parameters settings
    block_width = 256;
    block_height = 256;
    block_depth = 32; //

    outDatatype = 1; // 8-bit

    nbits = 4; // bit shift for 16-bit

    ubuffer = NULL;

    // inputs
    resolutions = scales;

    genMetaInfoOnly = genMetaInfo;

    inputImageDimensions = false;
    if(sx>0 && sy>0 && sz>0)
    {
        width = sx;
        height = sy;
        depth = sz;

        inputImageDimensions = true;

        color = 1;
        datatype = 2; // default

        if(block_depth > depth)
        {
            block_depth = depth;
        }
    }

    meta.cubex = block_width;
    meta.cubey = block_height;
    meta.cubez = block_depth;

    split = false;
    if(startz>=0 && endz>0)
    {
        zstart = startz;
        zend = endz;

        split = true;
    }

    if(resolutions<1)
    {
        cout<<"Invalide resolutions setting \n";
        return;
    }

    //
    srcdir.assign(inputdir);
    dstdir.assign(outputdir);

    //
    omp_set_num_threads(omp_get_max_threads());

    // generate meta info only
    if(genMetaInfo && inputImageDimensions)
    {
        index();
        exit(0);
    }

    //
    if(init())
    {
        cout<<"fail in init() \n";
        exit(-1);
    }

    //
    if(reformat())
    {
        cout<<"fail in reformat() \n";
        exit(-1);
    }
}

SnapTree::~SnapTree()
{
    input2DTIFFs.clear();
}

int SnapTree::init()
{
    // assuming the input folder has all images
    // SnapTree will reformat images between zstart and zend

    //
    DIR *indir = opendir(srcdir.c_str());
    if(indir == NULL)
    {
        cout<< srcdir <<": No such file or directory"  << endl;
        closedir(indir);
        return -1;
    }
    else
    {
        // get image list
        struct dirent *dirinfo = readdir(indir);
        while(dirinfo)
        {
            if(!strcmp(dirinfo->d_name,".") || !strcmp(dirinfo->d_name,".."))
            {
                dirinfo = readdir(indir);
                continue;
            }
            input2DTIFFs.insert(srcdir + "/" + dirinfo->d_name); // absolute path

            dirinfo = readdir(indir);
        }
        closedir(indir);
    }

    if(input2DTIFFs.size()<1)
    {
        cout<<"No TIFF file found \n";
        return -1;
    }

    //
    if(inputImageDimensions == false)
    {
        string firstfilepath = *input2DTIFFs.begin();
        loadTiffMetaInfo(const_cast<char*>(firstfilepath.c_str()), width, height, depth, color, datatype);
        depth = input2DTIFFs.size();

        if(block_depth > depth)
        {
            block_depth = depth;
            meta.cubez = block_depth;
        }

        cout<<"Image Info obtained from "<<firstfilepath<<endl;
    }

    if(split==false)
    {
        zstart = 0;
        zend = depth;
    }

    cout<<"Image Size "<<width<<"x"<<height<<"x"<<depth<<"x"<<color<<" with "<<datatype<<endl;

    //
    if(index())
    {
        cout<<"fail in generating mdata"<<endl;
        return -1;
    }

    // for the case with small z
    float w = width;
    float h = height;
    float d = depth;
    long n = 1;
    for(size_t i=0; i<resolutions; i++)
    {
        w *= 0.5;
        h *= 0.5;
        d *= 0.5;

        if(w>=1 && h>=1 && d>=1)
        {
            n++;
        }
        else
        {
            break;
        }
    }

    if(n<resolutions)
        resolutions = n;

    cout<<"resolutions "<<resolutions<<endl;

    // Make Hierarchical Dirs and blank TIFF images between (zstart, zend)
    map<int, Layer>::iterator res_iter = meta.layers.begin();
    while(res_iter != meta.layers.end())
    {
        //
        uint32 res_i = res_iter->first;

        uint32 startZ = (uint32)(zstart/(uint32)pow(2,res_i));
        uint32 endZ = (uint32)(zend/(uint32)pow(2,res_i));

        //
        Layer layer = (res_iter++)->second;

        int zs = getN(startZ, meta.cubez);
        int ze = getN(endZ, meta.cubez);

        cout<<"zs "<<zs<<" ze "<<ze<<" "<<layer.ncubes<<endl;

        //
        map<string, YXFolder>::iterator iter = layer.yxfolders.begin();
        while(iter != layer.yxfolders.end())
        {
            //
            YXFolder yxfolder = (iter++)->second;

            if(makeDir(yxfolder.xDirPath.c_str()))
            {
                cout<<"fail in mkdir "<<yxfolder.xDirPath<<endl;
                return -1;
            }

            if(makeDir(yxfolder.yDirPath.c_str()))
            {
                cout<<"fail in mkdir "<<yxfolder.yDirPath<<endl;
                return -1;
            }

            // process related cubes within [zstart, zend]
            for(int z=zs; z<ze; z++)
            {
                int offset = z*meta.cubez;

                auto it = yxfolder.cubes.find(offset);
                if(it != yxfolder.cubes.end())
                {
                    Cube cube = it->second;

                    // cube has not been created on the storage
                    struct stat info;
                    if( stat( cube.filePath.c_str(), &info ) != 0 )
                    {
                        if(initTiff3DFile((char *)(cube.filePath.c_str()),yxfolder.width,yxfolder.height,cube.depth,meta.color,1) != 0)
                        {
                            cout<<"fail in initTiff3DFile\n";
                            return -1;
                        }
                    }
                }
            } // cube
        } // yxfolder
    } // layer

    //
    return 0;
}

uint8 *SnapTree::load(long zs, long ze)
{
    //
    long sbv_V, sbv_H, sbv_D;

    sbv_V = height;
    sbv_H = width;
    sbv_D = ze - zs;

    //
    uint8 *subvol = NULL;

    try
    {
        subvol = new uint8 [sbv_V * sbv_H * sbv_D * datatype];
    }
    catch(...)
    {
        cout<<"failed to alloc memory for subvol \n";
        return NULL;
    }

    // fstream TIFFs from disk to memory
    vector<stringstream*> dataInMemory;
    vector<uint8*> imgList;

    //
    int k;
    for(k=0; k<sbv_D; k++)
    {
        //building image path
        string slicepath = *next(input2DTIFFs.begin(), zs + k);

        cout<<"load ... "<<slicepath<<endl;

        //
        ifstream inFile;
        inFile.open(slicepath.c_str());
        if (!inFile) {
            cerr << "Unable to open file "<<slicepath<<endl;
            return NULL;
        }

        //
        dataInMemory.push_back(new stringstream);
        *dataInMemory[k] << inFile.rdbuf();

        //
        inFile.close();

        //
        uint8 *slice = subvol + (k*sbv_V*sbv_H*datatype);
        imgList.push_back(slice);
    }

    // multithreaded read TIFFs from memory
    #pragma omp parallel
    {
        #pragma omp for
        for(k=0; k<sbv_D; k++)
        {
            unsigned int sx, sy;
            readTiff(dataInMemory[k],imgList[k],sx,sy,0,0,0,sbv_V-1,0,sbv_H-1);
        }
    }

    //
    return subvol;
}

int SnapTree::reformat()
{
    // reformat [zstart, zend] TIFF images

    //
    int zstep = MAX_IMAGES_STREAM;

    //
    for(int z=zstart; z<=zend; z+=zstep)
    {
        // load a mount of images [zs, ze]
        int zs = z;
        int ze = (z+zstep <= zend) ? (z+zstep) : zend;
        long zdepth = ze - zs;

        cout<<"reformat a sub volume ["<<zs<<", "<<ze<<"]"<<endl;

        //
        auto start_load = std::chrono::high_resolution_clock::now();

        ubuffer = load(zs, ze);

        auto end_load = std::chrono::high_resolution_clock::now();
        cout<<"load a sub volume takes "<<std::chrono::duration_cast<std::chrono::milliseconds>(end_load - start_load).count()<<" ms."<<endl;

        // bit-shift for 16-bit input data
        if(datatype>1 && nbits)
        {
            long totalvoxels = (height * width * zdepth) * color;
            if ( datatype == 2 )
            {
                #pragma omp parallel
                {
                    uint16 *ptr = (uint16 *) ubuffer;
                    long i;
                    #pragma omp for
                    for(i=0; i<totalvoxels; i++ )
                    {
                        ptr[i] = ptr[i] >> nbits;
                    }
                }
            }
        }

        //
        auto start_reformat = std::chrono::high_resolution_clock::now();

        //
        int pre_V, pre_H, pre_D;
        int cur_V, cur_H, cur_D;

        //
        map<int, Layer>::iterator res_iter = meta.layers.begin();
        while(res_iter != meta.layers.end())
        {
            //
            uint32 res_i = res_iter->first;

            if(res_i>0)
            {
                pre_V = cur_V;
                pre_H = cur_H;
                pre_D = cur_D;
            }

            //
            Layer layer = (res_iter++)->second;

            //
            int zs_i = zs / int(pow(2, res_i));
            int zdepth_i = zdepth / int(pow(2, res_i));

            cur_V = layer.dim_V;
            cur_H = layer.dim_H;
            cur_D = zdepth_i;

            // downsample
            if(res_i>0)
            {
                halveSample(ubuffer,pre_V,pre_H,pre_D,HALVE_BY_MAX,datatype);
            }

            cout<<"processing scale "<<res_i<<" z "<<zs_i<<" + "<<zdepth_i<<endl;

            //
            long imgwidth = layer.dim_H;

            //
            map<string, YXFolder>::iterator iter = layer.yxfolders.begin();
            while(iter != layer.yxfolders.end())
            {
                //
                YXFolder yxfolder = (iter++)->second;

                // locate the related cube: assuming [zs, ze] within a cube (does not cover many cubes)
                int zcube = zs_i / meta.cubez;

                int offset = zcube*meta.cubez;

                auto it = yxfolder.cubes.find(offset);
                if(it != yxfolder.cubes.end())
                {
                    long sx = yxfolder.width;
                    long sy = yxfolder.height;

                    long ystart = yxfolder.offset_V;
                    long xstart = yxfolder.offset_H;

                    Cube cube = it->second;

                    //cout<<"save sub volume to "<<cube.filePath<<endl;

                    //
                    void *fhandle = 0;
                    if(openTiff3DFile((char *)(cube.filePath.c_str()),(char *)("a"),fhandle,true))
                    {
                        cout<<"fail in openTiff3DFile "<<cube.filePath<<endl;
                        return -1;
                    }

                    long szChunk = yxfolder.width*yxfolder.height;
                    unsigned char *p = NULL;

                    try
                    {
                        p = new unsigned char [szChunk];
                        // memset(p, 0, szChunk);
                    }
                    catch(...)
                    {
                        cout<<"fail to alloc memory \n";
                        return -1;
                    }

                    // copy a 2D section and append it to the corresponding chunk image
                    for(long zbuf=0; zbuf<zdepth_i; zbuf++)
                    {
                        //
                        long offsetInput = zbuf*layer.dim_H*layer.dim_V;

                        //
                        if(datatype == 2)
                        {
                            // 16-bit input
                            uint16 *raw_ch16 = (uint16 *) ubuffer + offsetInput;

                            if(outDatatype == 1)
                            {
                                // 8-bit output

                                // copy
                                #pragma omp parallel for collapse(2)
                                for(long y=0; y<sy; y++)
                                {
                                    for(long x=0; x<sx; x++)
                                    {
                                        p[y*sx+x] = raw_ch16[(y+ystart)*imgwidth + (x+xstart)];
                                    }
                                }

                                // save
                                appendSlice2Tiff3DFile(fhandle,cube.offset_D+zbuf,(unsigned char *)p,sx,sy,meta.color,8,cube.depth);
                            }
                            else
                            {
                                // 16-bit output
                                uint16 *out_ch16 = (uint16 *) p;

                                // copy
                                #pragma omp parallel for collapse(2)
                                for(long y=0; y<sy; y++)
                                {
                                    for(long x=0; x<sx; x++)
                                    {
                                        out_ch16[y*sx+x] = raw_ch16[(y+ystart)*imgwidth + (x+xstart)];
                                    }
                                }

                                // save
                                appendSlice2Tiff3DFile(fhandle,cube.offset_D+zbuf,(unsigned char *)out_ch16,sx,sy,meta.color,16,cube.depth);
                            }

                        }
                        else if(datatype == 1)
                        {
                            // 8-bit input
                            uint8 *raw_ch8 = (uint8 *) ubuffer + offsetInput;

                            if(outDatatype == 1)
                            {
                                // 8-bit output

                                //cout<<"copying ... ("<<xstart<<", "<<ystart<<", "<<zbuf<<") -> "<<zs_i-cube.offset_D+zbuf<<" of "<<cube.depth<<endl;

                                // copy
                                #pragma omp parallel for collapse(2)
                                for(long y=0; y<sy; y++)
                                {
                                    for(long x=0; x<sx; x++)
                                    {
                                        p[y*sx+x] = raw_ch8[(y+ystart)*imgwidth + (x+xstart)];
                                    }
                                }

                                // save
                                appendSlice2Tiff3DFile(fhandle,zs_i-cube.offset_D+zbuf,(unsigned char *)p,sx,sy,meta.color,8,cube.depth);
                            }
                            else
                            {
                                // 16-bit output
                            }
                        }
                        else
                        {
                            // other datatypes
                        }
                    }

                    //
                    del1dp(p);

                    // close(fhandle) i.e. currently opened file
                    TIFFClose((TIFF *) fhandle);

                } // cube
            } // yxfolder
        } // layer

        // release allocated memory
        del1dp(ubuffer);

        //
        auto end_reformat = std::chrono::high_resolution_clock::now();
        cout<<"reformat and write chunk images takes "<<std::chrono::duration_cast<std::chrono::milliseconds>(end_reformat - start_reformat).count()<<" ms."<<endl;

    } // z

    //
    return 0;
}

int SnapTree::index()
{
    // saving mdata.bin for fast indexing image blocks instead of re-scan files every time

    // voxel size 1 micron by default
    // original offsets 0 mm by default

    //
    bool mDataDebug = false;

    //
    DIR *outdir = opendir(dstdir.c_str());
    if(outdir == NULL)
    {
        // mkdir outdir
        if(makeDir(dstdir.c_str()))
        {
            cout<<"fail in mkdir "<<dstdir<<endl;
            return -1;
        }
    }
    else
    {
        closedir(outdir);
    }

    auto start = std::chrono::high_resolution_clock::now();

    //
    meta.init(dstdir, width, height, depth, color, outDatatype, resolutions);

    //
    map<int, Layer>::iterator res_iter = meta.layers.begin();
    while(res_iter != meta.layers.end())
    {
        //
        Layer layer = (res_iter++)->second;

        //cout<<"res_i "<<res_iter->first<<endl;

        if(makeDir(layer.layerName.c_str()))
        {
            cout<<"fail in mkdir "<<layer.layerName<<endl;
            return -1;
        }

        //
        string filename = layer.layerName + "/mdata.bin";

        struct stat info;

        // mdata.bin does not exist
        if( stat( filename.c_str(), &info ) != 0 )
        {
            // save mdata.bin
            FILE *file;

            file = fopen(filename.c_str(), "w");

            fwrite(&(meta.mdata_version), sizeof(float), 1, file);
            fwrite(&(meta.reference_V), sizeof(axis), 1, file);
            fwrite(&(meta.reference_H), sizeof(axis), 1, file);
            fwrite(&(meta.reference_D), sizeof(axis), 1, file);
            fwrite(&(layer.vs_x), sizeof(float), 1, file);
            fwrite(&(layer.vs_y), sizeof(float), 1, file);
            fwrite(&(layer.vs_z), sizeof(float), 1, file);
            fwrite(&(layer.vs_x), sizeof(float), 1, file);
            fwrite(&(layer.vs_y), sizeof(float), 1, file);
            fwrite(&(layer.vs_z), sizeof(float), 1, file);
            fwrite(&(meta.org_V), sizeof(float), 1, file);
            fwrite(&(meta.org_H), sizeof(float), 1, file);
            fwrite(&(meta.org_D), sizeof(float), 1, file);
            fwrite(&(layer.dim_V), sizeof(uint32), 1, file);
            fwrite(&(layer.dim_H), sizeof(uint32), 1, file);
            fwrite(&(layer.dim_D), sizeof(uint32), 1, file);
            fwrite(&(layer.rows), sizeof(uint16), 1, file);
            fwrite(&(layer.cols), sizeof(uint16), 1, file);

            //
            if(mDataDebug)
            {
                cout<<"filename "<<filename<<endl;

                cout<<"meta.mdata_version "<<meta.mdata_version<<endl;
                cout<<"meta.reference_V "<<meta.reference_V<<endl;
                cout<<"meta.reference_H "<<meta.reference_H<<endl;
                cout<<"meta.reference_D "<<meta.reference_D<<endl;
                cout<<"layer.vs_x "<<layer.vs_x<<endl;
                cout<<"layer.vs_y "<<layer.vs_y<<endl;
                cout<<"layer.vs_z "<<layer.vs_z<<endl;
                cout<<"layer.vs_x "<<layer.vs_x<<endl;
                cout<<"layer.vs_y "<<layer.vs_y<<endl;
                cout<<"layer.vs_z "<<layer.vs_z<<endl;
                cout<<"meta.org_V "<<meta.org_V<<endl;
                cout<<"meta.org_H "<<meta.org_H<<endl;
                cout<<"meta.org_D "<<meta.org_D<<endl;
                cout<<"layer.dim_V "<<layer.dim_V<<endl;
                cout<<"layer.dim_H "<<layer.dim_H<<endl;
                cout<<"layer.dim_D "<<layer.dim_D<<endl;
                cout<<"layer.rows "<<layer.rows<<endl;
                cout<<"layer.cols "<<layer.cols<<endl;
            }

            //
            map<string, YXFolder>::iterator iter = layer.yxfolders.begin();
            while(iter != layer.yxfolders.end())
            {
                //
                YXFolder yxfolder = (iter++)->second;

                //
                fwrite(&(yxfolder.height), sizeof(uint32), 1, file);
                fwrite(&(yxfolder.width), sizeof(uint32), 1, file);
                fwrite(&(layer.dim_D), sizeof(uint32), 1, file); // depth of all blocks
                fwrite(&layer.ncubes, sizeof(uint32), 1, file);
                fwrite(&(meta.color), sizeof(uint32), 1, file);
                fwrite(&(yxfolder.offset_V), sizeof(int), 1, file);
                fwrite(&(yxfolder.offset_H), sizeof(int), 1, file);
                fwrite(&(yxfolder.lengthDirName), sizeof(uint16), 1, file);
                fwrite(const_cast<char *>(yxfolder.dirName.c_str()), yxfolder.lengthDirName, 1, file);

                //
                if(mDataDebug)
                {
                    cout<<"... "<<endl;
                    cout<<"HEIGHT "<<yxfolder.height<<endl;
                    cout<<"WIDTH "<<yxfolder.width<<endl;
                    cout<<"DEPTH "<<layer.dim_D<<endl;
                    cout<<"N_BLOCKS "<<layer.ncubes<<endl;
                    cout<<"N_CHANS "<<meta.color<<endl;
                    cout<<"ABS_V "<<yxfolder.offset_V<<endl;
                    cout<<"ABS_H "<<yxfolder.offset_H<<endl;
                    cout<<"str_size "<<yxfolder.lengthDirName<<endl;
                    cout<<"DIR_NAME "<<yxfolder.dirName<<endl;
                }

                //
                map<int, Cube>::iterator it = yxfolder.cubes.begin();
                while(it != yxfolder.cubes.end())
                {
                    //
                    Cube cube = (it++)->second;

                    //
                    fwrite(&(yxfolder.lengthFileName), sizeof(uint16), 1, file);
                    fwrite(const_cast<char *>(cube.fileName.c_str()), yxfolder.lengthFileName, 1, file);
                    fwrite(&(cube.depth), sizeof(uint32), 1, file);
                    fwrite(&(cube.offset_D), sizeof(int), 1, file);

                    //
                    if(mDataDebug)
                    {
                        cout<<"... ..."<<endl;
                        cout<<"str_size "<<yxfolder.lengthFileName<<endl;
                        cout<<"FILENAMES["<<it->first<<"] "<<cube.fileName<<endl;
                        cout<<"BLOCK_SIZE+i "<<cube.depth<<endl;
                        cout<<"BLOCK_ABS_D+i "<<cube.offset_D<<endl;
                    }
                }
                fwrite(&(meta.bytesPerVoxel), sizeof(uint32), 1, file);

                if(mDataDebug)
                {
                    cout<<"N_BYTESxCHAN "<<meta.bytesPerVoxel<<endl;
                }
            }
            fclose(file);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    cout<<"generated the meta info (mdata.bin) ... "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()<<" ms."<<endl;

    //
    return 0;
}

