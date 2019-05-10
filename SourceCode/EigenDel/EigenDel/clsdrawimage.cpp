#include "clsdrawimage.h"

ClsDrawImage::ClsDrawImage()
{

}

void ClsDrawImage::DrawDiagnal(Mat& img, St_ImgPos& stImgPos)
{
    //This is the normal case
    if(stImgPos.bPrimery) // the first time  --> we should think nothing -> that's why it is primary
    {
        int iPoint1X = (stImgPos.iCurRefStart - stImgPos.iBaseLine) * stImgPos.iAmplify;
        int iPoint1Y = iPoint1X + stImgPos.iVerticalOffset;  // like a square -->

        int iPoint2X = (stImgPos.iCurRefEnd - stImgPos.iBaseLine) * stImgPos.iAmplify;
        int iPoint2Y = iPoint2X + stImgPos.iVerticalOffset;

        stImgPos.iPreviousYStart = iPoint1Y < iPoint2Y ? iPoint1Y : iPoint2Y;
        stImgPos.iPreviousYEnd = iPoint1Y < iPoint2Y ? iPoint2Y : iPoint1Y;
        stImgPos.iPreRefStart = stImgPos.iCurRefStart;
        stImgPos.iPreRefEnd = stImgPos.iCurRefEnd;
        stImgPos.iPreQueryStart = stImgPos.iCurQueryStart;
        stImgPos.iPreQueryEnd = stImgPos.iCurQueryEnd;

        cv::Scalar cvTmpScalar(0,0,200);
        if(stImgPos.bGoodDualAlign)
            cvTmpScalar = cv::Scalar(200, 200, 0);
        else if(stImgPos.bBadDualAlign)
            cvTmpScalar = cv::Scalar(0, 0, 0);

        cv::line(img, cv::Point(iPoint1X, iPoint1Y), cv::Point(iPoint2X, iPoint2Y),
                 cvTmpScalar, 2, CV_AA);

        stImgPos.bPrimery = false;
        return;
    }

    //This case is for up to low --> not primary
    if(!stImgPos.bPrimery &&
       stImgPos.iPreRefEnd > 0 &&
       (stImgPos.iCurRefStart > stImgPos.iPreRefEnd))
    {
        int iPoint1X = (stImgPos.iCurRefStart - stImgPos.iBaseLine) * stImgPos.iAmplify;
        int iPoint1Y = stImgPos.iPreviousYEnd +
                       (stImgPos.iCurQueryStart - stImgPos.iPreQueryEnd) * stImgPos.iAmplify;

        int iPoint2X = (stImgPos.iCurRefEnd - stImgPos.iBaseLine) * stImgPos.iAmplify;
        int iPoint2Y = iPoint1Y + (stImgPos.iCurRefEnd - stImgPos.iCurRefStart) * stImgPos.iAmplify;

        cv::Scalar cvTmpScalar(0,200,0);
        if(stImgPos.bBadDualAlign)
            cvTmpScalar = cv::Scalar(0, 0, 0);

        cv::line(img, cv::Point(iPoint1X, iPoint1Y), cv::Point(iPoint2X, iPoint2Y),
                 cvTmpScalar, 2, CV_AA);
        return;
    }

    //This case is low to up
    if(!stImgPos.bPrimery &&
       stImgPos.iPreRefStart > 0 &&
       (stImgPos.iCurRefEnd < stImgPos.iPreRefStart))
    {
        int iPoint2X = (stImgPos.iCurRefEnd - stImgPos.iBaseLine) * stImgPos.iAmplify;
        int iPoint2Y = stImgPos.iPreviousYStart -
                       (stImgPos.iPreQueryStart - stImgPos.iCurQueryEnd) * stImgPos.iAmplify;

        int iPoint1X = (stImgPos.iCurRefStart - stImgPos.iBaseLine) * stImgPos.iAmplify;
        int iPoint1Y = iPoint2Y - (stImgPos.iCurRefEnd - stImgPos.iCurRefStart) * stImgPos.iAmplify;

        cv::Scalar cvTmpScalar(200,0,0);
        if(stImgPos.bBadDualAlign)
            cvTmpScalar = cv::Scalar(0, 0, 0);

        cv::line(img, cv::Point(iPoint1X, iPoint1Y), cv::Point(iPoint2X, iPoint2Y),
                 cvTmpScalar, 2, CV_AA);
        return;
    }
} //Let's put this logic into the main code of find SV -- > Go tomorrow

void ClsDrawImage::ShowImage(Mat& img)
{
    cv::namedWindow("drawing", CV_WINDOW_AUTOSIZE|CV_WINDOW_FREERATIO);
    cv::imshow("drawing", img);
    cv::waitKey(0);
}

void ClsDrawImage::SaveImage(Mat& img, string strFilePath)
{
    cv::imwrite(strFilePath.c_str(), img);

}
