// g++ -std=c++11 -fopenmp -o downsample2images downsample2images.cpp -ltiff
// gnayuy, 6/22/2018

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

#include <omp.h>

#include "tiffio.h"

using namespace std;

template <class InputType, class OutputType>
int maxdownsamplep(InputType *pimg1, InputType *pimg2, OutputType *p, long sx, long sy)
{
    long ssx = sx/2;
    long ssy = sy/2;
    
    // step 1. init with pimg1[2*y*sx + 2*x]
    
    long x,y;
#pragma omp parallel for private(x)
    for(y=0; y<ssy; y++)
    {
        for(x=0; x<ssx; x++)
        {
            p[y*ssx + x] = (OutputType) pimg1[2*y*sx + 2*x];
        }
    }
    
    // step 2. max ( p[y*ssx + x], pimg1[2*y*sx + 2*x+1] )
    
#pragma omp parallel for private(x)
    for(y=0; y<ssy; y++)
    {
        for(x=0; x<ssx; x++)
        {
            if(pimg1[2*y*sx + 2*x+1] > p[y*ssx + x])
                p[y*ssx + x] = (OutputType) pimg1[2*y*sx + 2*x+1];
        }
    }
    
    // step 3. max ( p[y*ssx + x], pimg1[(2*y+1)*sx + 2*x] )
    
#pragma omp parallel for private(x)
    for(y=0; y<ssy; y++)
    {
        for(x=0; x<ssx; x++)
        {
            if(pimg1[(2*y+1)*sx + 2*x] > p[y*ssx + x])
                p[y*ssx + x] = (OutputType) pimg1[(2*y+1)*sx + 2*x];
        }
    }
    
    // step 4. max ( p[y*ssx + x], pimg1[(2*y+1)*sx + 2*x+1] )
    
#pragma omp parallel for private(x)
    for(y=0; y<ssy; y++)
    {
        for(x=0; x<ssx; x++)
        {
            if(pimg1[(2*y+1)*sx + 2*x+1] > p[y*ssx + x])
                p[y*ssx + x] = (OutputType) pimg1[(2*y+1)*sx + 2*x+1];
        }
    }
    
    // step 5. max ( p[y*ssx + x], pimg2[2*y*sx + 2*x] )
    
#pragma omp parallel for private(x)
    for(y=0; y<ssy; y++)
    {
        for(x=0; x<ssx; x++)
        {
            if(pimg2[2*y*sx + 2*x] > p[y*ssx + x])
                p[y*ssx + x] = (OutputType) pimg2[2*y*sx + 2*x];
        }
    }
    
    // step 6. max ( p[y*ssx + x], pimg2[2*y*sx + 2*x+1] )
    
#pragma omp parallel for private(x)
    for(y=0; y<ssy; y++)
    {
        for(x=0; x<ssx; x++)
        {
            if(pimg2[2*y*sx + 2*x+1] > p[y*ssx + x])
                p[y*ssx + x] = (OutputType) pimg2[2*y*sx + 2*x+1];
        }
    }
    
    // step 7. max ( p[y*ssx + x], pimg2[(2*y+1)*sx + 2*x] )
    
#pragma omp parallel for private(x)
    for(y=0; y<ssy; y++)
    {
        for(x=0; x<ssx; x++)
        {
            if(pimg2[(2*y+1)*sx + 2*x] > p[y*ssx + x])
                p[y*ssx + x] = (OutputType) pimg2[(2*y+1)*sx + 2*x];
        }
    }
    
    // step 8. max ( p[y*ssx + x], pimg2[(2*y+1)*sx + 2*x+1] )
    
#pragma omp parallel for private(x)
    for(y=0; y<ssy; y++)
    {
        for(x=0; x<ssx; x++)
        {
            if(pimg2[(2*y+1)*sx + 2*x+1] > p[y*ssx + x])
                p[y*ssx + x] = (OutputType) pimg2[(2*y+1)*sx + 2*x+1];
        }
    }
    
    return 0;
}

//
/// testing openmp sections
//
// atom function
template <class InputType, class OutputType>
void maxdownsample_func(InputType *pimg1, InputType *pimg2, OutputType *poutput, long sx, long ssx, long endx, long endy)
{
    for(long y=0; y<endy; y++)
    {
        long ofy = y*ssx;

        long ofy1 = 2*y*sx;
        long ofy2 = (2*y+1)*sx;

        for(long x=0; x<endx; x++)
        {
            //
            int A = pimg1[ofy1 + 2*x];
            int B = pimg1[ofy1 + 2*x+1];
            if ( B > A ) A = B;

            B = pimg1[ofy2 + 2*x];
            if ( B > A ) A = B;

            B = pimg1[ofy2 + 2*x+1];
            if ( B > A ) A = B;

            B = pimg2[ofy1 + 2*x];
            if ( B > A ) A = B;

            B = pimg2[ofy1 + 2*x+1];
            if ( B > A ) A = B;

            B = pimg2[ofy2 + 2*x];
            if ( B > A ) A = B;

            B = pimg2[ofy2 + 2*x+1];
            if ( B > A ) A = B;

            // computing max
            poutput[ofy + x] = (OutputType)(A);
        }
    }
}

// parallel func
template <class InputType, class OutputType>
int maxdownsample_parallel(InputType *pimg1, InputType *pimg2, OutputType *poutput, long sx, long sy)
{
    long ssx = sx/2;
    long ssy = sy/2;

    long endx = ssx;
    long endy = ssy/8;

    #pragma omp parallel sections
    {
        // thread 1
        #pragma omp section
        maxdownsample_func<InputType,OutputType>(pimg1, pimg2, poutput, sx, ssx, endx, endy);

        // thread 2
        #pragma omp section
        maxdownsample_func<InputType,OutputType>(pimg1+endy*2*sx, pimg2+endy*2*sx, poutput+endy*ssx, sx, ssx, endx, endy);

        // thread 3
        #pragma omp section
        maxdownsample_func<InputType,OutputType>(pimg1+endy*4*sx, pimg2+endy*4*sx, poutput+endy*2*ssx, sx, ssx, endx, endy);

        // thread 4
        #pragma omp section
        maxdownsample_func<InputType,OutputType>(pimg1+endy*6*sx, pimg2+endy*6*sx, poutput+endy*3*ssx, sx, ssx, endx, endy);

        // thread 5
        #pragma omp section
        maxdownsample_func<InputType,OutputType>(pimg1+endy*8*sx, pimg2+endy*8*sx, poutput+endy*4*ssx, sx, ssx, endx, endy);

        // thread 6
        #pragma omp section
        maxdownsample_func<InputType,OutputType>(pimg1+endy*10*sx, pimg2+endy*10*sx, poutput+endy*5*ssx, sx, ssx, endx, endy);

        // thread 7
        #pragma omp section
        maxdownsample_func<InputType,OutputType>(pimg1+endy*12*sx, pimg2+endy*12*sx, poutput+endy*6*ssx, sx, ssx, endx, endy);

        // thread 8
        #pragma omp section
        maxdownsample_func<InputType,OutputType>(pimg1+endy*14*sx, pimg2+endy*14*sx, poutput+endy*7*ssx, sx, ssx, endx, endy);
    }

    //
    return 0;
}

// sequential max downsampling
template <class InputType, class OutputType>
int maxdownsamples(InputType *pimg1, InputType *pimg2, OutputType *poutput, long sx, long sy)
{
    long ssx = sx/2;
    long ssy = sy/2;
    
    long x,y;
    for(y=0; y<ssy; y++)
    {
        long ofy = y*ssx;
        
        long ofy1 = 2*y*sx;
        long ofy2 = (2*y+1)*sx;
        
        for(x=0; x<ssx; x++)
        {
            //
            int A = pimg1[ofy1 + 2*x];
            int B = pimg1[ofy1 + 2*x+1];
            if ( B > A ) A = B;
            
            B = pimg1[ofy2 + 2*x];
            if ( B > A ) A = B;
            
            B = pimg1[ofy2 + 2*x+1];
            if ( B > A ) A = B;
            
            B = pimg2[ofy1 + 2*x];
            if ( B > A ) A = B;
            
            B = pimg2[ofy1 + 2*x+1];
            if ( B > A ) A = B;
            
            B = pimg2[ofy2 + 2*x];
            if ( B > A ) A = B;
            
            B = pimg2[ofy2 + 2*x+1];
            if ( B > A ) A = B;
            
            // computing max
            poutput[ofy + x] = (OutputType)(A);
        }
    }
    
    return 0;
}

//
char *tiffread(char* filename, unsigned char *&p, uint32 &sz0, uint32  &sz1, uint32  &sz2, uint16 &datatype, uint16 &comp)
{
    //
    TIFF *input = TIFFOpen(filename,"r");
    if (!input)
    {
        return ((char *) "Cannot open the file.");
    }
    
    if (!TIFFGetField(input, TIFFTAG_IMAGEWIDTH, &sz0))
    {
        TIFFClose(input);
        return ((char *) "Image width of undefined.");
    }
    
    if (!TIFFGetField(input, TIFFTAG_IMAGELENGTH, &sz1))
    {
        TIFFClose(input);
        return ((char *) "Image length of undefined.");
    }
    
    if (!TIFFGetField(input, TIFFTAG_BITSPERSAMPLE, &datatype))
    {
        TIFFClose(input);
        return ((char *) "Undefined bits per sample.");
    }
    
    uint16 Cpage;
    if (!TIFFGetField(input, TIFFTAG_PAGENUMBER, &Cpage, &sz2) || sz2==0)
    {
        sz2 = 1;
        while ( TIFFReadDirectory(input) )
        {
            sz2++;
        }
    }
    datatype /= 8;

    long imgsz = (long)sz0*(long)sz1*(long)sz2*(long)datatype;
    
    //
    try
    {
        p = new unsigned char [imgsz];
    }
    catch(...)
    {
        return ((char*) "fail to alloc memory for loading a tiff image.");
    }
    
    //
    uint32 rps;
    int StripsPerImage,LastStripSize;
    
    //
    if (!TIFFGetField(input, TIFFTAG_ROWSPERSTRIP, &rps))
    {
        TIFFClose(input);
        return ((char *) "Undefined rowsperstrip.");
    }
    
    if (!TIFFGetField(input, TIFFTAG_COMPRESSION, &comp))
    {
        TIFFClose(input);
        return ((char *) "Undefined compression.");
    }
    
    StripsPerImage =  (sz1 + rps - 1) / rps;
    LastStripSize = sz1 % rps;
    if (LastStripSize==0)
        LastStripSize=rps;
    
    if (!TIFFSetDirectory(input, 0)) // init
    {
        TIFFClose(input);
        return ((char *) "fail to setdir.");
    }
    
    unsigned char *buf = p;

    do{

        for (int i=0; i < StripsPerImage-1; i++)
        {
            if (comp==1)
            {
                TIFFReadRawStrip(input, i, buf,  rps * sz0 * datatype);
                buf = buf + rps * sz0 * datatype;
            }
            else
            {
                TIFFReadEncodedStrip(input, i, buf, rps * sz0 * datatype);
                buf = buf + rps * sz0 * datatype;
            }
        }

        if (comp==1)
        {
            TIFFReadRawStrip(input, StripsPerImage-1, buf, LastStripSize * sz0 * datatype);
        }
        else
        {
            TIFFReadEncodedStrip(input, StripsPerImage-1, buf, LastStripSize * sz0 * datatype);
        }
        buf = buf + LastStripSize * sz0 * datatype;
        
    }
    while (TIFFReadDirectory(input)); // while (TIFFReadDirectory(input));
    
    //
    TIFFClose(input);
    
    //
    return ((char *) 0);
}

//
char *tiffwrite(char* filename, unsigned char *p, uint32 sz0, uint32  sz1, uint32  sz2, uint16 datatype, uint16 comp)
{
    
    TIFF *output = TIFFOpen(filename,"w");
    
    //
    if(sz2>1)
    {
        // 3D TIFF
        for(long slice=0; slice<sz2; slice++)
        {
            TIFFSetDirectory(output,slice);
            
            TIFFSetField(output, TIFFTAG_IMAGEWIDTH, sz0);
            TIFFSetField(output, TIFFTAG_IMAGELENGTH, sz1);
            TIFFSetField(output, TIFFTAG_BITSPERSAMPLE, (uint16)(datatype*8));
            TIFFSetField(output, TIFFTAG_SAMPLESPERPIXEL, 1);
            TIFFSetField(output, TIFFTAG_ROWSPERSTRIP, sz1);
            TIFFSetField(output, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
            TIFFSetField(output, TIFFTAG_COMPRESSION, comp);
            TIFFSetField(output, TIFFTAG_PLANARCONFIG,PLANARCONFIG_CONTIG);
            //TIFFSetField(output, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
            TIFFSetField(output, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
            
            // We are writing single page of the multipage file
            TIFFSetField(output, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
            TIFFSetField(output, TIFFTAG_PAGENUMBER, (uint16)slice, (uint16)sz2);
            
            // the file has been already opened: rowsPerStrip it is not too large for this image width
            TIFFWriteEncodedStrip(output, 0, p, sz0 * sz1 * datatype);
            
            TIFFWriteDirectory(output);
        }
    }
    else
    {
        // 2D TIFF
        TIFFSetField(output, TIFFTAG_IMAGEWIDTH, sz0);
        TIFFSetField(output, TIFFTAG_IMAGELENGTH, sz1);
        TIFFSetField(output, TIFFTAG_BITSPERSAMPLE, (uint16)(datatype*8));
        TIFFSetField(output, TIFFTAG_SAMPLESPERPIXEL, 1);
        TIFFSetField(output, TIFFTAG_ROWSPERSTRIP, sz1);
        TIFFSetField(output, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
        TIFFSetField(output, TIFFTAG_COMPRESSION, comp);
        TIFFSetField(output, TIFFTAG_PLANARCONFIG,PLANARCONFIG_CONTIG);
        TIFFSetField(output, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
        
        // the file has been already opened: rowsPerStrip it is not too large for this image width
        TIFFWriteEncodedStrip(output, 0, p, sz0 * sz1 * datatype);
        
        //
        TIFFWriteDirectory(output);
    }
    
    //
    TIFFClose(output);
    
    //
    return ((char *) 0);
}

//
int main(int argc, char *argv[])
{
    if(argc<4)
    {
        cout<<"./imagemaxds image1.tif image2.tif output.tif"<<endl;
        return -1;
    }
    
    //
    unsigned char *pImg1, *pImg2, *p;
    uint32 sx1, sx2, sy1, sy2, sz1, sz2, sx, sy;
    long imgsz, outsz;
    uint16 datatype1, datatype2;
    uint16 comp;
    
    //
    auto start = std::chrono::high_resolution_clock::now();
    
    //
    char *error_check1 = tiffread(argv[1], pImg1, sx1, sy1, sz1, datatype1, comp);
    if(error_check1)
    {
        cout<<error_check1<<endl;
        return -1;
    }
    
    char *error_check2 = tiffread(argv[2], pImg2, sx2, sy2, sz2, datatype2, comp);
    if(error_check2)
    {
        cout<<error_check2<<endl;
        return -1;
    }
    
    //
    if(datatype1 != datatype2)
    {
        cout<<"datatypes are not matched"<<endl;
        return -1;
    }
    
    //
    if(sz1 != sz2 || sy1 != sy2 || sx1 != sx2)
    {
        cout<<"image dimensions are not matched"<<endl;
        return -1;
    }
    
    //
    sx = sx1 / 2;
    sy = sy1 / 2;
    
    try
    {
        outsz = sx*sy;
        p = new unsigned char [outsz];
        // memset(p, 0, outsz);
    }
    catch(...)
    {
        cout<<"fail to alloc memory for out image\n";
        return -1;
    }
    
    //
    if(datatype1 == 2)
    {
        unsigned short *p1 = (unsigned short *) (pImg1);
        unsigned short *p2 = (unsigned short *) (pImg2);
        
        // bit-shift
        int nbits = 2;

        #pragma omp parallel
        {
            imgsz = sx1*sy1;
            long i;
            #pragma omp for
            for (i=0; i<imgsz; i++)
            {
                p1[i] = p1[i] >> nbits;
                p2[i] = p2[i] >> nbits;

                if(p1[i]>255)
                    p1[i] = 255;

                if(p2[i]>255)
                    p2[i] = 255;
            }
        }
        
        // downsample
        maxdownsample_parallel<unsigned short, unsigned char>(p1, p2, p, sx1, sy1);
    }
    else if(datatype1 == 1)
    {
        // downsample
        unsigned char *p1 = (unsigned char *) (pImg1);
        unsigned char *p2 = (unsigned char *) (pImg2);
        
        maxdownsamplep<unsigned char, unsigned char>(p1, p2, p, sx1, sy1);
    }
    else
    {
        cout<<"unsupported data type"<<endl;
        return -1;
    }
    
    // save result
    char *error_check = tiffwrite(argv[3], p, sx, sy, 1, 1, comp);
    if(error_check)
    {
        cout<<error_check<<endl;
        return -1;
    }
    
    //
    auto end = std::chrono::high_resolution_clock::now();
    cout<<"downsample 2 images takes "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()<<" ms."<<endl;
    
    //
    if(pImg1)
    {
        delete []pImg1;
    }
    
    if(pImg2)
    {
        delete []pImg2;
    }
    
    if(p)
    {
        delete []p;
    }
    
    //
    return 0;
}

