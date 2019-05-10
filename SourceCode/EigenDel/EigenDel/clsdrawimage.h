#ifndef CLSDRAWIMAGE_H
#define CLSDRAWIMAGE_H

#include <iostream>
#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include <stdio.h>

using namespace std;
using namespace cv;

struct St_ImgPos
{
    int iCurRefStart;
    int iCurRefEnd;

    int iPreRefStart;
    int iPreRefEnd;

    int iCurQueryStart;
    int iCurQueryEnd;

    int iPreQueryEnd;
    int iPreQueryStart;

    int iVerticalOffset;
    int iPreviousYEnd;
    int iPreviousYStart;

    int iAmplify;
    bool bPrimery; // This is always the first part -->

    int iBaseLine; // For reference --> for example: subject sequence

    bool bGoodDualAlign;
    bool bBadDualAlign;

    St_ImgPos()
    {
        Clear();
        iAmplify = -1;
        iBaseLine = -1;
        bGoodDualAlign = false;
        bBadDualAlign = false;
    }

    void Clear()
    {
        iCurRefStart = -1;
        iCurRefEnd = -1;

        iPreRefStart = -1;
        iPreRefEnd = -1;

        iCurQueryStart = -1;
        iCurQueryEnd = -1;

        iPreQueryEnd = -1;
        iPreQueryStart = -1;

        iVerticalOffset = -1;
        iPreviousYEnd = -1;
        iPreviousYStart = -1;

        bPrimery = false;
    }
};

class ClsDrawImage
{
public:
    ClsDrawImage();

public:
    void DrawDiagnal(Mat& img, St_ImgPos& stImgPos);
    void ShowImage(Mat& img);
    void SaveImage(Mat& img, string strFilePath);
};

#endif // CLSDRAWIMAGE_H
