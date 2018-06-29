// SnapTree
// Reformatting large-scale dataset into a hierarchical file organization
// gnayuy, 6/22/2018
//
//
// Here is simple implementation to demonstrate the idea of SnapTree
//


//
# include <stdio.h>
# include <stdlib.h>
#include <iostream>

#include "cxxopts.hpp"
#include "SnapTree.h"

//
int main (int argc, const char *argv[])
{
    // Usage:
    // SnapTree -i <input_DIR> -o <output_DIR> -n <Number_of_Resolutions_of_TMITREE>

    // assuming input 2D TIFF (LZW compressed) images and convert to 3D TIFF blocks
    // 3D block with 256x256x256

    if(argc<2)
    {
        cout<<"snaptree version 1.0\n";
        cout<<"snaptree -h\n";
        return 0;
    }

    //
    int scales = 3;
    string inputDir, outputDir;
    int mdata = 0;
    int sx=0, sy=0, sz=0;
    int startz=-1, endz=-1;

    //
    try
    {
        cxxopts::Options options(argv[0], "snaptree -i <input_DIR> -o <output_DIR> -n <Number_of_Resolutions>");
        options
                .positional_help("[optional args]")
                .show_positional_help();

        options.add_options()
                ("h,help", "SnapTree Version 1.0")
                ("i,input", "Input DIR", cxxopts::value<std::string>(inputDir))
                ("o,output", "Output DIR", cxxopts::value<std::string>(outputDir))
                ("n,resolutions", "N Resolutions", cxxopts::value<int>(scales))
                ("m,meta", "generate meta info only", cxxopts::value<int>(mdata))
                ("x,dimx", "x-dimension", cxxopts::value<int>(sx))
                ("y,dimy", "y-dimension", cxxopts::value<int>(sy))
                ("z,dimz", "z-dimension", cxxopts::value<int>(sz))
                ("s,startz", "start z-slice", cxxopts::value<int>(startz))
                ("e,endz", "end z-slice", cxxopts::value<int>(endz))
                ;

        auto cmds = options.parse(argc, argv);

        if (cmds.count("help"))
        {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (cmds.count("input"))
        {
            std::cout << " -- Input 2D TIFF Images DIR: " << cmds["input"].as<std::string>() << std::endl;
        }

        if (cmds.count("output"))
        {
            std::cout << " -- Output 3D TIFF Images DIR: " << cmds["output"].as<std::string>() << std::endl;
        }

        if (cmds.count("resolutions"))
        {
            std::cout << " -- Convert data into " << cmds["resolutions"].as<int>() << " scales" << std::endl;
        }

        if (cmds.count("meta"))
        {
            std::cout << " -- Generate mdata.bin " << cmds["meta"].as<int>() << std::endl;
        }

        if (cmds.count("dimx"))
        {
            std::cout << " -- Input x-dimension " << cmds["dimx"].as<int>() << std::endl;
        }

        if (cmds.count("dimy"))
        {
            std::cout << " -- Input y-dimension " << cmds["dimy"].as<int>() << std::endl;
        }

        if (cmds.count("dimz"))
        {
            std::cout << " -- Input z-dimension " << cmds["dimz"].as<int>() << std::endl;
        }

        if (cmds.count("startz"))
        {
            std::cout << " -- the starting z-slice " << cmds["startz"].as<int>() << std::endl;
        }

        if (cmds.count("endz"))
        {
            std::cout << " -- the end z-slice " << cmds["endz"].as<int>() << std::endl;
        }

    }
    catch(const cxxopts::OptionException& e)
    {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }

    // BigTree
    SnapTree snaptree(inputDir, outputDir, scales, mdata, sx, sy, sz, startz, endz);

    //
    return 0;
}
