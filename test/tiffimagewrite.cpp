// g++ -o tiffimagewrite tiffimagewrite.cpp -ltiff

//
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <climits>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <string>
#include <tuple>
#include <stack>
#include <fstream>
#include <iostream>
#include <stdio.h>

// #include <omp.h>

#include "tiffio.h"

using namespace std;


//
int main(int argc, char *argv[])
{
    if(argc<2)
    {
        cout<<"tiffimagewrite output.tif"<<endl;
        return -1;
    }
    
    //
    unsigned char *p = NULL;
    
    uint32 width = 256;
    uint32 length = 256;
    uint32 depth = 256;
    
    long imgsz = width*length*depth;
    
    uint32  rowsperstrip = (uint32) -1;
    uint16    photometric = PHOTOMETRIC_MINISBLACK;
    uint16    config = PLANARCONFIG_CONTIG;
    
    TIFFDataType dtype = TIFF_BYTE;
    
    try {
        p = new unsigned char [imgsz];
        memset(p, 0, imgsz);
    } catch (...) {
        cout<<"fail to allocate memory for p"<<endl;
        return -1;
    }
    
    //
    TIFF *out = TIFFOpen(argv[1], "w");
    
    
    if(depth>1)
    {
        // 3D
        for(uint32 slice=0; slice<depth; slice++)
        {
            TIFFSetDirectory(out,slice);
            
            TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);
            TIFFSetField(out, TIFFTAG_IMAGELENGTH, length);
            TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);
            TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1);
            TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, length);
            TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
            TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
            TIFFSetField(out, TIFFTAG_PLANARCONFIG,PLANARCONFIG_CONTIG);
            //TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
            TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
            
            // We are writing single page of the multipage file
            TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
            TIFFSetField(out, TIFFTAG_PAGENUMBER, slice, depth);
            
            // the file has been already opened: rowsPerStrip it is not too large for this image width
            TIFFWriteEncodedStrip(out, 0, p, width*length);
            
            TIFFWriteDirectory(out);
        }
        
    }
    else
    {
        // 2D
        //
        TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);
        TIFFSetField(out, TIFFTAG_IMAGELENGTH, length);
        TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
        TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1);
        TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);
        // TIFFSetField(out, TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);
        TIFFSetField(out, TIFFTAG_PLANARCONFIG, config);
        TIFFSetField(out, TIFFTAG_PHOTOMETRIC, photometric);
        
        TIFFSetField(out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
        
        TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
        
        TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, length);
        
        //    for (uint32 row = 0; row < length; row++)
        //    {
        //        unsigned char * buf = p + row*width;
        //        if (TIFFWriteScanline(out, p, row, 0) < 0)
        //        {
        //            cout<<"Failed in writing row "<<row<<endl;
        //            return -1;
        //        }
        //    }
        
        //TIFFWriteEncodedStrip(out, 0, p, width*length);
        
        if(TIFFWriteTile(out, p, width, length, 0, 0) < 0)
        {
            cout<<"Failed in writing tiled tiff "<<endl;
        }
    }
    
    //
    TIFFClose(out);
    
    
    //
    if(p)
    {
        delete []p;
    }
    
    //
    return 0;
}
