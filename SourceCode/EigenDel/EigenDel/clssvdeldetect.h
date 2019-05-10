#ifndef CLSSVDELDETECT_H
#define CLSSVDELDETECT_H
#include "clsparsebam.h"


struct St_GroupBound
{
    vector<St_BorderSCR> vSCR;
    bool bLeftSVBoundGood;
    bool bRightSVBoundGood;
    int iClipPos;

    St_GroupBound()
    {
        Clear();
    }

    void Clear()
    {
        vSCR.clear();
        bLeftSVBoundGood = false;
        bRightSVBoundGood = false;
        iClipPos = 0;
    }
};

enum En_DepthType{dtZeroPure=0, dtZeroWeek, dtZeroMix, dtRepWeek, dtRepStrong, dtMax};

struct St_DepthSeg
{
    int iZeroPureDepth;
    int iZeroPureDepthLen;

    float fZeroWeekDepth;
    int iZeroWeekDepthLen;

    float fZeroMixDepth;
    int iZeroMixDepthLen;

    float fRepWeekDepth;
    int iRepWeekDepthLen;

    float fRepStrongDepth;
    int iRepStrongDepthLen;

    float fAvgDepth;
    int iAvgDepthLen;

    float fStatDepth;
    int iStatDepthLen;
    float fStatLenRatio;

    En_DepthType enDepthType;

    int iSegLen;

    St_DepthSeg()
    {
        iZeroPureDepth = 0; // this is fix
        iZeroPureDepthLen = 0;

        fZeroWeekDepth = -1;
        iZeroWeekDepthLen = 0;

        fZeroMixDepth = -1;
        iZeroMixDepthLen = 0;

        fRepWeekDepth = -1;
        iRepWeekDepthLen = 0;

        fRepStrongDepth = -1;
        iRepStrongDepthLen = 0;

        fAvgDepth = -1;
        iAvgDepthLen = 0;

        fStatDepth = -1;
        iStatDepthLen = 0;

        fStatLenRatio = 0;

        enDepthType = dtMax;

        iSegLen = -1;
    }

    int GetRegDepthLen()
    {
        return iSegLen - iZeroPureDepthLen - iZeroWeekDepthLen - iZeroMixDepthLen; //"iZeroWeekDepthLen" already includes "iZeroPureDepthLen"
    }

    string GetZeroPure()
    {
        return to_string(iZeroPureDepth) + "/" + to_string(iZeroPureDepthLen);
    }

    string GetZeroWeek()
    {
        return FloatToStr(fZeroWeekDepth) + "/" + to_string(iZeroWeekDepthLen);
    }

    string GetZeroMix()
    {
        return FloatToStr(fZeroMixDepth) + "/" + to_string(iZeroMixDepthLen);
    }

    string GetRepWeek()
    {
        return FloatToStr(fRepWeekDepth) + "/" + to_string(iRepWeekDepthLen);
    }

    string GetRepStrong()
    {
        return FloatToStr(fRepStrongDepth) + "/" + to_string(iRepStrongDepthLen);
    }

    string GetAvg()
    {
        return FloatToStr(fAvgDepth) + "/" + to_string(iAvgDepthLen);
    }

    string GetStat()
    {
        return FloatToStr(fStatDepth) + "/" + to_string(iStatDepthLen);
    }

    float GetZeroPureDepthLenRatio()
    {
        return (float)this->iZeroPureDepthLen / this->iSegLen;
    }

    float GetZeroWeekDepthLenRatio()
    {
        return (float)this->iZeroWeekDepthLen / this->iSegLen;
    }

    float GetZeroMixDepthLenRatio()
    {
        return (float)this->iZeroMixDepthLen / this->iSegLen;
    }

    float GetRepWeekDepthLenRatio()
    {
        return (float)this->iRepWeekDepthLen / this->iSegLen;
    }

    float GetRepStrongDepthLenRatio()
    {
        return (float)this->iRepStrongDepthLen / this->iSegLen;
    }
};

//Features from each group:
struct St_GroupFeatures
{
    int iStart;
    int iEnd;

    //For partital special depth range
    St_DepthSeg stDCRangeCore;
    //float fPartialDepth;
    //int iDepthRangLen;

    //purpose: get the regular coverage:
    //so the partial depth should significantly lower or hight than those the reg depth
    St_DepthSeg stDCRangeLeft;
    St_DepthSeg stDCRangeRight;

//    float fAvgRegDepthLeft;
//    int iRegDepthLeftLen;
//    float fAvgRegDepthRight;
//    int iRegDepthRightLen;

    St_DepthSeg stNoSvRangeLeft;
    St_DepthSeg stNoSvRangeRight;
    //float fNoSvDepthLeft; // Left from the left boundary
    //float fNoSvDepthRight; // Right from the right boundary

    int iStartClipReadsNum;
    int iEndClipReadsNum;
    int iDRReadsNum;
    int iSpeciReadsNum;

     St_GroupFeatures()
     {
         Clear();
     }     

     void Clear()
     {
         iStart = -1;
         iEnd = -1;

//         fPartialDepth = -1;
//         iDepthRangLen = -1;

//         fAvgRegDepthLeft = -1;
//         iRegDepthLeftLen = -1;
//         fAvgRegDepthRight = -1;
//         iRegDepthRightLen = -1;

//         fNoSvDepthLeft = -1;
//         fNoSvDepthRight = -1;

         iStartClipReadsNum = 0;
         iEndClipReadsNum = 0;
         iDRReadsNum = 0;
         iSpeciReadsNum = 0;
     }

     static bool sort_depth_small_large_func(float f1, float f2)
     {
         if(f1 < f2)
             return true;
         else
             return false;
     }

//     float GetMaxTargetDepth()
//     {
//         vector<float> vDepth{fPartialDepth, fAvgRegDepthLeft, fAvgRegDepthRight, fNoSvDepthLeft, fNoSvDepthRight};
//         sort(vDepth.begin(), vDepth.end(), sort_depth_small_large_func);
//         return *(vDepth.end() - 1);
//     }

//     float GetMinTargetDepth()
//     {
//         vector<float> vDepth{fPartialDepth, fAvgRegDepthLeft, fAvgRegDepthRight, fNoSvDepthLeft, fNoSvDepthRight};
//         sort(vDepth.begin(), vDepth.end(), sort_depth_small_large_func);
//         return *vDepth.begin();
//     }
};

// use the similar way of delly to make clusrering
struct St_BaseGroup
{
    int iStart; // valid range for SV    
    int iEnd; // valid range for SV

    int iLeftBoundary; //For real left boundary
    int iRightBoundary; //For real right boundary

    int iRealSVStart;
    int iRealSVEnd;

    int iStartAlignLen;
    int iLeftBoundaryAlignLen;

    bool bStartClipAdjusted;
    bool bEndClipAdjusted;

    St_GroupFeatures stGF; // group features

    St_BaseGroup()
    {
        Clear();
    }

    void Clear()
    {
        iStart = 0;
        iEnd = 0;
        iLeftBoundary = 0;
        iRightBoundary = 0;
        iRealSVStart = -1;
        iRealSVEnd = -1;
        iStartAlignLen = -1;
        iLeftBoundaryAlignLen = -1;
        stGF.Clear();
        bStartClipAdjusted = false;
        bEndClipAdjusted = false;
    }

    int GetLen()
    {
        return iRightBoundary - iLeftBoundary + 1;
    }

    int GetAccuLen()
    {
        int iLen = iEnd - (iStart + iStartAlignLen) + 1 ;
        if(bStartClipAdjusted)
            iLen = iEnd - iStart + 1 ;
        return iLen;
    }

    int GetStart()
    {
        if(bStartClipAdjusted)
            return iStart;
        else
            return iStart + iStartAlignLen;
    }
};

// Discordant reads group
enum En_GroupType{gtNegative=0, gtPositive, gtMax};

struct St_DRGroup: public St_BaseGroup
{
    //-->
    vector<St_DiscordantReads> vDR;  //Left-most first, Right-most last
    vector<St_BorderSCR> vDRRangedSCR; // SCR dropped into DR range
    vector<St_RegReads> vRegReads; // For one mapped and the mate unmapped
    //<--    
    St_SV stSv;
    En_GroupType enType;

    St_DRGroup()
    {
        Clear();
    }

    void Clear()
    {
        vDR.clear();
        vDRRangedSCR.clear();
        vRegReads.clear();
        St_BaseGroup::Clear();
        stSv.Clear();
        enType = gtMax;
    }
};

struct St_SCRGroup: public St_BaseGroup
{
    // -->
    vector<St_BorderSCR> vStartSCR;  // Start Border of Clipped Reads
    vector<St_BorderSCR> vEndSCR; //End Border of Clipped Reads
    // <--

    St_SCRGroup()
    {
        Clear();
    }

    void Clear()
    {
        vStartSCR.clear();
        vEndSCR.clear();
        St_BaseGroup::Clear();
    }
};

struct St_SCRBound
{
    vector<St_BorderSCR> vSCR;
    int iCount;
    int iClipPos;
    int iDiff;
    int iClosestMat;

    St_SCRBound()
    {
        iCount = 0;
        iClipPos = -1;
        iDiff = 100000;
        iClosestMat = 100000;
    }

    void Clear()
    {
        vSCR.clear();
        iCount = 0;
        iClipPos = -1;
        iDiff = 100000;
        iClosestMat = 100000;
    }
};

struct St_FinalResult
{
    int iPositiveNum;
    int iNegativeNum;
    int iStdNum;
    string strChrom;

    St_FinalResult()
    {
        Clear();
    }

    void Clear()
    {
        iPositiveNum = 0;
        iNegativeNum = 0;
        iStdNum = 0;
        strChrom = "";
    }
};

class ClsSvDelDetect
{
public:
    ClsSvDelDetect();
    ~ClsSvDelDetect();

public:
    void ClusterDiscrodantReadsDelly(vector<St_DRGroup>& vDRGroup, vector<St_DiscordantReads>& vDiscdRreads,
                                     vector<St_BorderSCR>& vBorderSCR);

    void ClusterSCReads(vector<St_SCRGroup>& vSCRGroup, vector<St_BorderSCR>& vBorderSCR);
    void GetStdSvDelRelatedGroup(vector<St_DRGroup>& vDRGroup, vector<St_SCRGroup>& vSCRGroup, vector<St_SV>& vSvDEL);

    void DepthFilter(vector<St_DRGroup>& vDRGroup, vector<St_SCRGroup>& vSCRGroup,
                     string& strBamFile, string strChrom);

    void GetSpeciReads(vector<St_DRGroup>& vDRGroup, string strBamFile);
    void CollectFeatures(vector<St_DRGroup>& vDRGroup, string strBamFile, string strChrom);

public:
    void SetMeanInsertSize(float fV1);
    void SetStdDevInsertSize(float fV1);
    void SetAvgDepth(float fV1);
    void SetChromName(string strChromName);
    void SetOfstream(ofstream& ofs);

    //void GetDepthOfRange(St_DRGroup& stDRGroup, string strBamFile);
    void GetDepthOfRange(string strBamFile,
                         int iStartPos, int iEndPos, int iRangeLen,
                         float& fTargetRangDepth, int& iTargetRangLen, string strChrom);

    void GetDepthOfRangeEX(string strBamFile,
                           int iStartPos, int iEndPos, int iRangeLen,
                           St_DepthSeg& stDepthSeg, string strChrom);

private:
    void UpgradeSCR(vector<St_BorderSCR>& vBorderSCR, vector<St_BorderSCR>& vNewSCR);
    void GetSCRGroup(vector<St_GroupBound>& vGroupBound, vector<St_SCRGroup>& vSCRGroup);
    bool GetGroupFromBoundary(vector<St_BorderSCR>& vSCR, St_SCRGroup& stSCRGroup);

    bool DepthFilterGroup(St_BaseGroup& stGroup, string& strBamFile,
                          float fThreshold, string strChrom);

    void FilterBadDR(vector<St_DiscordantReads>& vDiscdRreads);
    void UpdateStartEndByDR(St_DRGroup& stDRGroup);
    void UpdateBoundBySCR(St_DRGroup& stDRGroup);

    void GetDepthFeatures(St_DRGroup& stDRGroup, string strBamFile, string strChrom);
    void GetClipFeatures(St_DRGroup& stDRGroup);

    void UpdateStartEndBySameClipOrMate(St_DRGroup& stDRGroup,
                                        vector<St_SCRBound>& vLeftBoundSCR, vector<St_SCRBound>& vRightBoundSCR,
                                        bool& bLeftUpdated, bool& bRightUpdated);

private:
    float GetMinValidIS();
    float GetMaxValidIS();

private:
    float m_fMeanInsertSize;
    float m_fStdDevInsertSize;
    float m_fAvgDepth;
    string m_strChromName;
    ofstream* m_pOfs;
};

#endif // CLSSVDELDETECT_H
