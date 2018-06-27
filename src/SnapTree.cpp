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

        layer.dim_V = (uint32)(dimy/pow(2,res_i));
        layer.dim_H = (uint32)(dimx/pow(2,res_i));
        layer.dim_D = (uint32)(dimz/pow(2,res_i));

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

    meta.cubex = block_width;
    meta.cubey = block_height;
    meta.cubez = block_depth; //

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
    }

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

        uint32 startZ = (uint32)(zstart/pow(2,res_i));
        uint32 endZ = (uint32)(zend/pow(2,res_i));

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

                    if(initTiff3DFile((char *)(cube.filePath.c_str()),yxfolder.width,yxfolder.height,cube.depth,meta.color,1) != 0)
                    {
                        cout<<"fail in initTiff3DFile\n";
                        return -1;
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
    // meta
//    vector<LAYER> layers;

//    //
//    int stack_block[TMITREE_MAX_HEIGHT];
//    int slice_start[TMITREE_MAX_HEIGHT];
//    int slice_end[TMITREE_MAX_HEIGHT];

//    int nzsize[TMITREE_MAX_HEIGHT];

//    for(int res_i=0; res_i< resolutions; res_i++)
//    {
//        stack_block[res_i] = 0;
//        slice_start[res_i] = 0;
//        slice_end[res_i] = slice_start[res_i] + stacks_D[res_i][0][0][0] - 1;

//        //cout<<"slice_end["<<res_i<<"] "<<slice_end[res_i]<<endl;

//        nzsize[res_i] = 0;
//    }

//    //
//    for(long z=0, z_parts=1; z<depth; z+=z_max_res, z_parts++)
//    {
//        if(!genMetaInfoOnly && !genZeroDataOnly)
//        {
//            auto start = std::chrono::high_resolution_clock::now();

//            ubuffer = load(z,(z+z_max_res <= depth) ? (z+z_max_res) : depth);

//            auto end = std::chrono::high_resolution_clock::now();

//            cout<<"load a sub volume takes "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()<<" ms."<<endl;

//            //
//            if(datatype>1 && nbits)
//            {
//                long totalvoxels = (height * width * ((z_ratio>0) ? z_max_res : (depth%z_max_res)))*color;
//                if ( datatype == 2 )
//                {
//                    #pragma omp parallel
//                    {
//                        uint16 *ptr = (uint16 *) ubuffer;
//                        #pragma omp for
//                        for(long i=0; i<totalvoxels; i++ )
//                        {
//                            // ptr[i] = ptr[i] >> nbits << nbits; // 16-bit
//                            ptr[i] = ptr[i] >> nbits;
//                        }
//                    }
//                }
//            }

//        } // genMetaInfoOnly

//        //saving the sub volume
//        auto start = std::chrono::high_resolution_clock::now();
//        for(int i=0; i< resolutions; i++)
//        {
//            //cout<<"resolution "<<i<<endl;

//            // meta info for index()
//            LAYER layer;
//            layer.rows = n_stacks_V[i];
//            layer.cols = n_stacks_H[i];

//            layer.vs_x = pow(2, i);
//            layer.vs_y = pow(2, i);
//            layer.vs_z = pow(2, halve_pow2[i]);

//            layer.dim_V = (uint32)(height/layer.vs_y);
//            layer.dim_H = (uint32)(width/layer.vs_x);
//            layer.dim_D = (uint32)(depth/layer.vs_z);

//            //
//            long nCopies = 0;
//            string srcFile;

//            // check if current block is changed
//            if ( (z / pow(2,halve_pow2[i])) > slice_end[i] ) {
//                stack_block[i]++;
//                slice_start[i] = slice_end[i] + 1;
//                slice_end[i] += stacks_D[i][0][0][stack_block[i]];
//            }

//            // find abs_pos_z at resolution i
//            std::stringstream abs_pos_z;
//            abs_pos_z.width(6);
//            abs_pos_z.fill('0');
//            abs_pos_z << (int)((pow(2,halve_pow2[i])*slice_start[i]) * 10);

//            // compute the number of slice of previous groups at resolution i
//            // note that z_parts in the number and not an index (starts from 1)
//            long n_slices_pred  = (z_parts - 1) * z_max_res / pow(2,halve_pow2[i]);

//            // buffer size along D is different when the remainder of the subdivision by z_max_res is considered
//            long z_size = (z_ratio>0) ? z_max_res : (depth%z_max_res);

//            //cout<<"z_parts "<<z_parts<<" z_ratio "<<z_ratio<<" z_max_res "<<z_max_res<<" depth "<<depth<<endl;

//            //halvesampling resolution if current resolution is not the deepest one
//            if(i!=0)
//            {
//                if(!genMetaInfoOnly && !genZeroDataOnly)
//                {
//                    if ( halve_pow2[i] == (halve_pow2[i-1]+1) )
//                    {
//                        //cout<<"3D downsampling \n";

//                        // 3D
//                        halveSample(ubuffer,(int)height/(pow(2,i-1)),(int)width/(pow(2,i-1)),(int)z_size/(pow(2,halve_pow2[i-1])),HALVE_BY_MAX,datatype);

//                        // debug
//                        // writeTiff3DFile("test.tif", ubuffer, (int)width/(pow(2,i)), (int)height/(pow(2,i)), (int)z_size/(pow(2,halve_pow2[i])), 1, datatype);
//                    }
//                    else if ( halve_pow2[i] == halve_pow2[i-1] )
//                    {
//                        //cout<<"2D downsampling \n";

//                        // 2D
//                        halveSample2D(ubuffer,(int)height/(pow(2,i-1)),(int)width/(pow(2,i-1)),(int)z_size/(pow(2,halve_pow2[i-1])),HALVE_BY_MAX,datatype);
//                    }
//                    else
//                    {
//                        cout<<"halve sampling level "<<halve_pow2[i]<<" not supported at resolution "<<i<<endl;
//                        return -1;
//                    }
//                }
//            }

//            //saving at current resolution if it has been selected and iff buffer is at least 1 voxel (Z) deep
//            if((z_size/(pow(2,halve_pow2[i]))) > 0)
//            {
//                //storing in 'base_path' the absolute path of the directory that will contain all stacks
//                std::stringstream base_path;
//                base_path << filePaths[i] << "/";

//                //cout<<"base_path "<<base_path.str()<<endl;

//                // meta info for index()
////                bool addMeta = false;

////                if(z == 0)
////                {
////                    addMeta = true;
////                }
////                else if( (int)(slice_start[i] / stacks_D[i][0][0][stack_block[i]]) > nzsize[i])
////                {
////                    nzsize[i]++;
////                    addMeta = true;
////                }

//                //looping on new stacks
//                for(int stack_row = 0, start_height = 0, end_height = 0; stack_row < n_stacks_V[i]; stack_row++)
//                {
//                    //incrementing end_height
//                    end_height = start_height + stacks_V[i][stack_row][0][0]-1;

//                    //computing V_DIR_path and creating the directory the first time it is needed
//                    std::stringstream multires_merging_x_pos;
//                    multires_merging_x_pos.width(6);
//                    multires_merging_x_pos.fill('0');
//                    multires_merging_x_pos << start_height*(int)pow(2.0,i) * 10;

//                    std::stringstream V_DIR_path;
//                    V_DIR_path << base_path.str() << multires_merging_x_pos.str();

//                    //cout<<"V_DIR_path "<<V_DIR_path.str()<<endl;

//                    if(z==0)
//                    {
//                        if(makeDir(V_DIR_path.str().c_str()))
//                        {
//                            cout<<" unable to create V_DIR "<<V_DIR_path.str()<<endl;
//                            return -1;
//                        }
//                    }

//                    int sz[4];
//                    // int datatype_out = 2; // changed to 16-bit 6/1/2018 yy
//                    int datatype_out = 1;

//                    //
//                    for(int stack_column = 0, start_width=0, end_width=0; stack_column < n_stacks_H[i]; stack_column++)
//                    {
//                        //
//                        end_width  = start_width  + stacks_H[i][stack_row][stack_column][0]-1;

//                        //computing H_DIR_path and creating the directory the first time it is needed
//                        std::stringstream multires_merging_y_pos;
//                        multires_merging_y_pos.width(6);
//                        multires_merging_y_pos.fill('0');
//                        multires_merging_y_pos << start_width*(int)pow(2.0,i) * 10;

//                        std::stringstream H_DIR_path;
//                        H_DIR_path << V_DIR_path.str() << "/" << multires_merging_x_pos.str() << "_" << multires_merging_y_pos.str();

//                        //cout<<"H_DIR_path "<<H_DIR_path.str()<<endl;

//                        //
//                        if(z==0)
//                        {
//                            if(makeDir(H_DIR_path.str().c_str()))
//                            {
//                                cout<<" unable to create H_DIR "<<H_DIR_path.str()<<endl;
//                                return -1;
//                            }
//                            else
//                            {
//                                // the directory has been created for the first time
//                                // initialize block files
//                                sz[0] = stacks_H[i][stack_row][stack_column][0];
//                                sz[1] = stacks_V[i][stack_row][stack_column][0];
//                                sz[3] = 1;

//                                //
//                                int slice_start_temp = 0;
//                                for ( int j=0; j < n_stacks_D[i]; j++ )
//                                {
//                                    bool copying = false;

//                                    if(sz[2] == stacks_D[i][stack_row][stack_column][j])
//                                        copying = true;

//                                    sz[2] = stacks_D[i][stack_row][stack_column][j];

//                                    //cout<<" ... "<<j<<" sz[2] "<<sz[2]<<" slice_end "<<slice_end[i]<<endl;

//                                    std::stringstream abs_pos_z_temp;
//                                    abs_pos_z_temp.width(6);
//                                    abs_pos_z_temp.fill('0');
//                                    abs_pos_z_temp << (int)((pow(2,halve_pow2[i])*slice_start_temp) * 10);

//                                    std::stringstream img_path_temp;
//                                    img_path_temp << H_DIR_path.str() << "/" << multires_merging_x_pos.str() << "_" << multires_merging_y_pos.str() << "_" << abs_pos_z_temp.str()<<".tif";

//                                    //cout<<"when z=0: z "<<z<<" ("<<sz[0]<<", "<<sz[1]<<", "<<sz[2]<<") "<<abs_pos_z_temp.str()<<endl;

//                                    if(!genMetaInfoOnly)
//                                    {
//                                        // auto start_init = std::chrono::high_resolution_clock::now();
//                                        if(nCopies==0)
//                                        {
//                                            if(initTiff3DFile((char *)img_path_temp.str().c_str(),sz[0],sz[1],sz[2],sz[3],datatype_out) != 0)
//                                            {
//                                                cout<<"fail in initTiff3DFile\n";
//                                                return -1;
//                                            }
//                                            srcFile = img_path_temp.str();
//                                        }
//                                        else if(copying)
//                                        {
//                                            copyFile(srcFile.c_str(), img_path_temp.str().c_str());
//                                        }
//                                        else
//                                        {
//                                            if(initTiff3DFile((char *)img_path_temp.str().c_str(),sz[0],sz[1],sz[2],sz[3],datatype_out) != 0)
//                                            {
//                                                cout<<"fail in initTiff3DFile\n";
//                                                return -1;
//                                            }
//                                        }
//                                        nCopies++;
//                                    }// genMetaInfoOnly

//                                    //
//                                    //auto end_init = std::chrono::high_resolution_clock::now();
//                                    //cout<<"writing chunk images takes "<<std::chrono::duration_cast<std::chrono::milliseconds>(end_init - start_init).count()<<" ms."<<endl;

//                                    //
//                                    slice_start_temp += (int)sz[2];
//                                } // j
//                            } // after create H_DIR
//                        } // z

//                        //saving HERE

//                        //
//                        std::stringstream partial_img_path;
//                        partial_img_path << H_DIR_path.str() << "/" << multires_merging_x_pos.str() << "_" << multires_merging_y_pos.str() << "_";

//                        int slice_ind = (int)(n_slices_pred - slice_start[i]);

//                        std::stringstream img_path;
//                        img_path << partial_img_path.str() << abs_pos_z.str() << ".tif";

//                        //cout<<"img_path "<<img_path.str()<<endl;

//                        //
//                        void *fhandle = 0;
//                        int  n_pages_block = stacks_D[i][0][0][stack_block[i]]; // number of pages of current block
//                        bool block_changed = false; // true if block is changed executing the next for cycle

//                        if(!genMetaInfoOnly && !genZeroDataOnly)
//                        {
//                            if(openTiff3DFile((char *)img_path.str().c_str(),(char *)("a"),fhandle,true))
//                            {
//                                cout<<"fail in openTiff3DFile"<<endl;
//                                return -1;
//                            }
//                        }// genMetaInfoOnly

//                        //
//                        sz[0] = end_width - start_width + 1;
//                        sz[1] = end_height - start_height + 1;
//                        sz[2] = n_pages_block;
//                        sz[3] = 1;
//                        long szChunk = sz[0]*sz[1]*sz[3]*datatype_out;
//                        unsigned char *p = NULL;

//                        if(!genMetaInfoOnly && !genZeroDataOnly)
//                        {
//                            try
//                            {
//                                p = new unsigned char [szChunk];
//                                memset(p, 0, szChunk);
//                            }
//                            catch(...)
//                            {
//                                cout<<"fail to alloc memory \n";
//                            }
//                        }// genMetaInfoOnly

//                        // meta info for index()
////                        BLOCK block;

////                        if(addMeta)
////                        {
////                            block.width = sz[0];
////                            block.height = sz[1];
////                            block.depths.push_back(sz[2]);
////                            block.color = sz[3];
////                            block.bytesPerVoxel = datatype_out;

////                            block.dirName = multires_merging_x_pos.str() + "/" + multires_merging_x_pos.str() + "_" + multires_merging_y_pos.str();
////                            block.offset_H = start_width;
////                            block.offset_V = start_height;
////                            block.fileNames.push_back(multires_merging_x_pos.str() + "_" + multires_merging_y_pos.str() + "_" + abs_pos_z.str() + ".tif");
////                            block.offsets_D.push_back(slice_start[i]);
////                        }

//                        bool blocksaved = false;

//                        //cout<<"z "<<z<<endl;

//                        // WARNING: assumes that block size along z is not less that z_size/(powInt(2,i))
//                        for(int buffer_z=0; buffer_z<(int)(z_size/(pow(2,halve_pow2[i]))); buffer_z++, slice_ind++)
//                        {
//                            //cout<<"buffer_z "<<buffer_z<<" slice_ind "<<slice_ind<<" z "<<z<<" z_size/(pow(2,halve_pow2[i]) "<<z_size/(pow(2,halve_pow2[i]))<<" z_size "<<z_size<<endl;
//                            //cout<<"(z / pow(2,halve_pow2[i]) + buffer_z) "<<z / pow(2,halve_pow2[i]) + buffer_z<<" slice_end["<<i<<"] "<<slice_end[i]<<endl;

//                            // z is an absolute index in volume while slice index should be computed on a relative basis
//                            if ( (int)(z / pow(2,halve_pow2[i]) + buffer_z) > slice_end[i] && !block_changed)
//                            {
//                                //cout<<"block changed "<<slice_end[i]<<endl;

//                                // start a new block along z
//                                std::stringstream abs_pos_z_next;
//                                abs_pos_z_next.width(6);
//                                abs_pos_z_next.fill('0');
//                                abs_pos_z_next << (pow(2,halve_pow2[i])*(slice_end[i]+1)) * 10;
//                                img_path.str("");
//                                img_path << partial_img_path.str() << abs_pos_z_next.str() << ".tif";

//                                //cout<<"... img_path "<<img_path.str()<<endl;

//                                slice_ind = 0;

//                                if(!genMetaInfoOnly && !genZeroDataOnly)
//                                {
//                                    // close(fhandle) i.e. file corresponding to current block
//                                    TIFFClose((TIFF *) fhandle);
//                                    if(openTiff3DFile((char *)img_path.str().c_str(),(char *)("a"),fhandle,true))
//                                    {
//                                        cout<<"fail in openTiff3DFile"<<endl;
//                                        return -1;
//                                    }
//                                }// genMetaInfoOnly
//                                n_pages_block = stacks_D[i][0][0][stack_block[i]+1];
//                                block_changed = true;

//                                sz[2] = n_pages_block;
//                                szChunk = sz[0]*sz[1]*sz[3]*datatype_out;

//                                //
//                                if(!genMetaInfoOnly && !genZeroDataOnly)
//                                {
//                                    if(!p)
//                                    {
//                                        try
//                                        {
//                                            p = new unsigned char [szChunk];
//                                            memset(p, 0, szChunk);
//                                        }
//                                        catch(...)
//                                        {
//                                            cout<<"fail to alloc memory \n";
//                                        }
//                                    }
//                                    else
//                                    {
//                                        memset(p, 0, szChunk);
//                                    }
//                                }// genMetaInfoOnly
//                            }

//                            //
//                            if(!genMetaInfoOnly && !genZeroDataOnly)
//                            {
//                                //
//                                long raw_img_width = width/(pow(2,i));

//                                //
//                                if(datatype == 2)
//                                {
//                                    // 16-bit input
//                                    long offset = buffer_z*(long)(height/pow(2,i))*(long)(width/pow(2,i));
//                                    uint16 *raw_ch16 = (uint16 *) ubuffer + offset;

//                                    //cout<<"pointer p: "<<static_cast<void*>(p)<<endl;
//                                    //cout<<"pointer raw data: "<<static_cast<void*>(raw_ch16)<<endl;

//                                    if(datatype_out == 1)
//                                    {
//                                        // 8-bit output

//                                        //
//                                        #pragma omp parallel for collapse(2)
//                                        for(long y=0; y<sz[1]; y++)
//                                        {
//                                            for(long x=0; x<sz[0]; x++)
//                                            {
//                                                p[y*sz[0]+x] = raw_ch16[(y+start_height)*(raw_img_width) + (x+start_width)];
//                                            }
//                                        }

//                                        // temporary save all the way (version 1.01 5/25/2018)
////                                        int temp_n_chans = color;
////                                        if(temp_n_chans==2)
////                                            temp_n_chans++;

//                                        appendSlice2Tiff3DFile(fhandle,slice_ind,(unsigned char *)p,sz[0],sz[1],color,8,sz[2]);
//                                        blocksaved = true;

//                                        //
////                                        int numNonZeros = 0;
////                                        int saveVoxelThresh = 1;

////                                        #pragma omp parallel for reduction(+:numNonZeros)
////                                        for(int x=0; x<szChunk; x++)
////                                        {
////                                            if(p[x]>0)
////                                                numNonZeros++;
////                                        }

////                                        if(numNonZeros>saveVoxelThresh)
////                                        {
////                                            int temp_n_chans = color;
////                                            if(temp_n_chans==2)
////                                                temp_n_chans++;

////                                            appendSlice2Tiff3DFile(fhandle,slice_ind,(unsigned char *)p,sz[0],sz[1],temp_n_chans,8,sz[2]);
////                                            blocksaved = true;
////                                        }
//                                    }
//                                    else
//                                    {
//                                        // 16-bit output

//                                        uint16 *out_ch16 = (uint16 *) p;

//                                        //
//                                        #pragma omp parallel for collapse(2)
//                                        for(long y=0; y<sz[1]; y++)
//                                        {
//                                            for(long x=0; x<sz[0]; x++)
//                                            {
//                                                out_ch16[y*sz[0]+x] = raw_ch16[(y+start_height)*(raw_img_width) + (x+start_width)];
//                                            }
//                                        }

//                                        // temporary save all the way (version 1.01 5/25/2018)
//                                        int temp_n_chans = color;
//                                        if(temp_n_chans==2)
//                                            temp_n_chans++;

//                                        appendSlice2Tiff3DFile(fhandle,slice_ind,(unsigned char *)out_ch16,sz[0],sz[1],temp_n_chans,16,sz[2]);
//                                        blocksaved = true;

//                                    }

//                                }
//                                else if(datatype == 1)
//                                {
//                                    // 8-bit input
//                                    long offset = buffer_z*(long)(height/pow(2,i))*(long)(width/pow(2,i));
//                                    uint8 *raw_ch8 = (uint8 *) ubuffer + offset;

//                                    if(datatype_out == 1)
//                                    {
//                                        // 8-bit output

//                                        //
//                                        #pragma omp parallel for collapse(2)
//                                        for(long y=0; y<sz[1]; y++)
//                                        {
//                                            for(long x=0; x<sz[0]; x++)
//                                            {
//                                                p[y*sz[0]+x] = raw_ch8[(y+start_height)*(raw_img_width) + (x+start_width)];
//                                            }
//                                        }

//                                        // temporary save all the way (version 1.01 5/25/2018)
//                                        int temp_n_chans = color;
//                                        if(temp_n_chans==2)
//                                            temp_n_chans++;

//                                        appendSlice2Tiff3DFile(fhandle,slice_ind,(unsigned char *)p,sz[0],sz[1],temp_n_chans,8,sz[2]);
//                                        blocksaved = true;

//                                        //
////                                        int numNonZeros = 0;
////                                        int saveVoxelThresh = 1;

////                                        #pragma omp parallel for reduction(+:numNonZeros)
////                                        for(int x=0; x<szChunk; x++)
////                                        {
////                                            if(p[x]>0)
////                                                numNonZeros++;
////                                        }

////                                        //cout<<"... raw_img_width "<<raw_img_width<<" offset "<<offset<<" height/pow(2,i) "<<height/pow(2,i)<<" width/pow(2,i) "<<width/pow(2,i)<<endl;

////                                        if(numNonZeros>saveVoxelThresh)
////                                        {
////                                            int temp_n_chans = color;
////                                            if(temp_n_chans==2)
////                                                temp_n_chans++;

////                                            //cout<<"... save slice_ind: "<<slice_ind<<endl;
////                                            appendSlice2Tiff3DFile(fhandle,slice_ind,(unsigned char *)p,sz[0],sz[1],temp_n_chans,8,sz[2]);
////                                            blocksaved = true;
////                                        }
//                                    }
//                                    else
//                                    {
//                                        // 16-bit output

//                                    }
//                                }
//                                else
//                                {
//                                    // other datatypes
//                                }

//                            }// genMetaInfoOnly
//                        }

//                        if(!genMetaInfoOnly && !genZeroDataOnly)
//                        {
//                            //
//                            del1dp(p);

//                            // close(fhandle) i.e. currently opened file
//                            TIFFClose((TIFF *) fhandle);
//                        }

//                        //
//                        start_width  += stacks_H[i][stack_row][stack_column][0];

//                        //
////                        if(addMeta)
////                        {
////                            block.nonZeroBlocks.push_back(blocksaved);
////                            block.nBlocksPerDir = block.fileNames.size();
////                            layer.blocks.push_back(block);
////                            layer.n_scale = i;
////                        }

//                    }
//                    start_height += stacks_V[i][stack_row][0][0];
//                }
//            }

//            //
////            if(!layer.blocks.empty())
////            {
////                layers.push_back(layer);
////            }
//        }

//        //releasing allocated memory
//        del1dp(ubuffer);

//        //
//        auto end = std::chrono::high_resolution_clock::now();
//        cout<<"writing sub volume's chunk images takes "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()<<" ms."<<endl;
//    }

    //
    return 0;
}

int SnapTree::index()
{
    // saving mdata.bin for fast indexing image blocks instead of re-scan files every time

    // voxel size 1 micron by default
    // original offsets 0 mm by default

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
    meta.init(dstdir, width, height, depth, color, datatype, resolutions);

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

        // save
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

//        cout<<"filename "<<filename<<endl;

//        cout<<"meta.mdata_version "<<meta.mdata_version<<endl;
//        cout<<"meta.reference_V "<<meta.reference_V<<endl;
//        cout<<"meta.reference_H "<<meta.reference_H<<endl;
//        cout<<"meta.reference_D "<<meta.reference_D<<endl;
//        cout<<"layer.vs_x "<<layer.vs_x<<endl;
//        cout<<"layer.vs_y "<<layer.vs_y<<endl;
//        cout<<"layer.vs_z "<<layer.vs_z<<endl;
//        cout<<"layer.vs_x "<<layer.vs_x<<endl;
//        cout<<"layer.vs_y "<<layer.vs_y<<endl;
//        cout<<"layer.vs_z "<<layer.vs_z<<endl;
//        cout<<"meta.org_V "<<meta.org_V<<endl;
//        cout<<"meta.org_H "<<meta.org_H<<endl;
//        cout<<"meta.org_D "<<meta.org_D<<endl;
//        cout<<"layer.dim_V "<<layer.dim_V<<endl;
//        cout<<"layer.dim_H "<<layer.dim_H<<endl;
//        cout<<"layer.dim_D "<<layer.dim_D<<endl;
//        cout<<"layer.rows "<<layer.rows<<endl;
//        cout<<"layer.cols "<<layer.cols<<endl;

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

//            cout<<"... "<<endl;
//            cout<<"HEIGHT "<<yxfolder.height<<endl;
//            cout<<"WIDTH "<<yxfolder.width<<endl;
//            cout<<"DEPTH "<<layer.dim_D<<endl;
//            cout<<"N_BLOCKS "<<layer.ncubes<<endl;
//            cout<<"N_CHANS "<<meta.color<<endl;
//            cout<<"ABS_V "<<yxfolder.offset_V<<endl;
//            cout<<"ABS_H "<<yxfolder.offset_H<<endl;
//            cout<<"str_size "<<yxfolder.lengthDirName<<endl;
//            cout<<"DIR_NAME "<<yxfolder.dirName<<endl;


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

//                cout<<"... ..."<<endl;
//                cout<<"str_size "<<yxfolder.lengthFileName<<endl;
//                cout<<"FILENAMES["<<it->first<<"] "<<cube.fileName<<endl;
//                cout<<"BLOCK_SIZE+i "<<cube.depth<<endl;
//                cout<<"BLOCK_ABS_D+i "<<cube.offset_D<<endl;

            }
            fwrite(&(meta.bytesPerVoxel), sizeof(uint32), 1, file);

            //cout<<"N_BYTESxCHAN "<<meta.bytesPerVoxel<<endl;
        }
        fclose(file);
    }

    auto end = std::chrono::high_resolution_clock::now();
    cout<<"generated the meta info (mdata.bin) ... "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()<<" ms."<<endl;

    //
    return 0;
}

