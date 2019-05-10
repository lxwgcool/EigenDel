#include "clslearning.h"
#include "clssvdeldetect.h"
#include <algorithm>
#include <cmath>

ClsLearning::ClsLearning()
{
    m_fAvgDepth = -1;
    m_strChromName = "";
    m_strSampleName = "";
    m_pOfs = NULL;
}

float ClsLearning::GetAvgDepth(){}

void ClsLearning::SetAvgDepth(float fV1)
{
    this->m_fAvgDepth = fV1;
}

void ClsLearning::SetChromName(string strChromName)
{
    this->m_strChromName = strChromName;
}

void ClsLearning::SetSampleName(string strSampleName)
{
    this->m_strSampleName = strSampleName;
}

void ClsLearning::SetOfstream(ofstream& ofs)
{
    this->m_pOfs = &ofs;
}

bool sort_positive_group_by_sv_small2large(St_DRGroup st1, St_DRGroup st2)
{
    if(st1.stSv.GetLen() <= st2.stSv.GetLen())
        return true;
    else
        return false;
}

void ClsLearning::PrintFeatures(vector<St_DRGroup>& vDRGroup, vector<St_SCRGroup>& vSCRGroup, vector<St_SV>& vSvDEL)
{
    //(1) Create DRGroup For each SvDEL
    //What we need to collect:
    //Step 1: Collect Positive and Negative Position
    vector<St_DRGroup> vPositiveGroup;
    vector<St_DRGroup> vNegativeGroup;
    GetPositiveAndNegativeGroup(vDRGroup, vSvDEL, vPositiveGroup, vNegativeGroup); //vSCRGroup
    sort(vPositiveGroup.begin(), vPositiveGroup.end(), sort_positive_group_by_sv_small2large);
    //Step 2: Save Features into the file       

    //Save the data for positive samples
    ofstream ofs;
    ofs.open(("./Output/PositiveSamples_Chrom_" + m_strChromName + ".txt").c_str());
    ofstream ofsDepthInfo;
    ofsDepthInfo.open(("./Output/PositiveSamplesDepthInfo_Chrom_" + m_strChromName + ".txt").c_str());
    //ofsDepthInfo1.open("./PositiveSamplesDepthInfo.txt");
    //ofsDepthInfo2.open("./PositiveSamplesDepthInfo.txt");

    ofstream ofsAbandonDepthInfo;
    ofsAbandonDepthInfo.open(("./Output/AbandonPositiveSamplesDepthInfo_Chrom_" + m_strChromName + ".txt").c_str());

    (*m_pOfs) << "Positive Samples" << endl;
    int iAbandonNum = 0;
    int iGoodNum = 0;

    int arryDT[dtMax+1] = {0, 0, 0, 0, 0, 0};
    string arryDTName[dtMax+1] = {"ZeroPure", "ZeroWeek", "ZeroMix", "RepWeek", "RepStrong", "dtMax"};

    int iIndex = 0;
    for(vector<St_DRGroup>::iterator itr = vPositiveGroup.begin(); itr != vPositiveGroup.end(); itr++, iIndex++)
    {
        if(itr->stGF.stDCRangeCore.iStatDepthLen < 10)
        {
            iAbandonNum++;
            PrintCandiDepthDetail(ofsAbandonDepthInfo, itr, iIndex);
            continue;
        }

        arryDT[itr->stGF.stDCRangeCore.enDepthType]++;

        //iStart, iEnd, fPartialDepth, iDepthRangLen, fAvgRegDepthLeft, fAvgRegDepthRight, iStartClipReadsNum, iEndClipReadsNum, iDRReadsNum, iSpeciReadsNum
        //iStart and iEnd should turns to the length: abs(iEnd - iStart)
        /*
         * This is for positive samples:
         * So the real info should be like below:
         * 1: SV_start, SV_End, SV_Length
         * 2: Len (iEnd - iStart), fPartialDepth, iDepthRangLen, fAvgRegDepthLeft, fAvgRegDepthRight, iStartClipReadsNum, iEndClipReadsNum, iDRReadsNum, iSpeciReadsNum
         */
        int iStrictLen = itr->stGF.iEnd - (itr->stGF.iStart + itr->iStartAlignLen);
        if(itr->bStartClipAdjusted)
            iStrictLen = itr->stGF.iEnd - itr->stGF.iStart;;

        ofs << to_string(itr->stSv.iPos) << "\t" << to_string(itr->stSv.iEnd) << "\t" << to_string(itr->stSv.iEnd - itr->stSv.iPos) << " --- "
            << to_string(iStrictLen) << "\t"
            << FloatToStr(itr->stGF.stDCRangeCore.fStatDepth) << "\t"
            << to_string(itr->stGF.stDCRangeCore.iStatDepthLen) << "\t"

            << FloatToStr(itr->stGF.stDCRangeLeft.fStatDepth) << "\t"
            << to_string(itr->stGF.stDCRangeLeft.iStatDepthLen) << "\t"

            << FloatToStr(itr->stGF.stDCRangeRight.fStatDepth) << "\t"
            << to_string(itr->stGF.stDCRangeRight.iStatDepthLen) << "\t"

            << FloatToStr(itr->stGF.stNoSvRangeLeft.fAvgDepth) << "\t"
            << FloatToStr(itr->stGF.stNoSvRangeRight.fAvgDepth) << "\t"

            << to_string(itr->stGF.iStartClipReadsNum) << "\t"
            << to_string(itr->stGF.iEndClipReadsNum) << "\t"
            << to_string(itr->stGF.iDRReadsNum) << "\t"
            << to_string(itr->stGF.iSpeciReadsNum) << endl;

        (*m_pOfs) << to_string(iIndex) << " -- \t"
             << FloatToStr(itr->stGF.stNoSvRangeLeft.fAvgDepth) << ", "
             << FloatToStr(itr->stGF.stDCRangeLeft.fStatDepth) << ", "
             << FloatToStr(itr->stGF.stDCRangeCore.fStatDepth) << ", "
             << FloatToStr(itr->stGF.stDCRangeRight.fStatDepth) << ", "
             << FloatToStr(itr->stGF.stNoSvRangeRight.fAvgDepth) << endl;

        PrintCandiDepthDetail(ofsDepthInfo, itr, iIndex);
        iGoodNum++;
//        //For No Sv Left
//        ofsDepthInfo << itr->stGF.stNoSvRangeLeft.GetAvg() << " --- "
//                     << to_string(itr->GetAccuLen()) << " --- "
//                     << arryDTName[itr->stGF.stDCRangeCore.enDepthType] << " --- "
//                     << "(" << arryDTName[itr->stGF.stDCRangeLeft.enDepthType] << ", "
//                            << arryDTName[itr->stGF.stDCRangeRight.enDepthType] << ")" << endl;

//        //For DC Sv Left
//        ofsDepthInfo << itr->stGF.stDCRangeLeft.GetZeroPure() << ",  "
//                     << itr->stGF.stDCRangeLeft.GetZeroWeek() << ",  "
//                     << itr->stGF.stDCRangeLeft.GetZeroMix() << ",  "
//                     << itr->stGF.stDCRangeLeft.GetRepWeek() << ",  "
//                     << itr->stGF.stDCRangeLeft.GetRepStrong() << ",  "
//                     << itr->stGF.stDCRangeLeft.GetStat() << " ---  " << endl;

//        //For DC sv Core
//        ofsDepthInfo << itr->stGF.stDCRangeCore.GetZeroPure() << ",  "
//                     << itr->stGF.stDCRangeCore.GetZeroWeek() << ",  "
//                     << itr->stGF.stDCRangeCore.GetZeroMix() << ",  "
//                     << itr->stGF.stDCRangeCore.GetRepWeek() << ",  "
//                     << itr->stGF.stDCRangeCore.GetRepStrong() << ",  "
//                     << itr->stGF.stDCRangeCore.GetStat() << " ---  " << endl;

//        //For DC sv Right
//        ofsDepthInfo << itr->stGF.stDCRangeRight.GetZeroPure() << ",  "
//                     << itr->stGF.stDCRangeRight.GetZeroWeek() << ",  "
//                     << itr->stGF.stDCRangeRight.GetZeroMix() << ",  "
//                     << itr->stGF.stDCRangeRight.GetRepWeek() << ",  "
//                     << itr->stGF.stDCRangeRight.GetRepStrong() << ",  "
//                     << itr->stGF.stDCRangeRight.GetStat() << " ---  " << endl;;

//        //For No sv Right
//        ofsDepthInfo << itr->stGF.stNoSvRangeRight.GetAvg() << endl << endl << endl;

//                     << FloatToStr(itr->stGF.stNoSvRangeLeft.fAvgDepth) << ", "
//                     << FloatToStr(itr->stGF.stDCRangeLeft.fStatDepth) << ", "
//                     << FloatToStr(itr->stGF.stDCRangeCore.fStatDepth) << ", "
//                     << FloatToStr(itr->stGF.stDCRangeRight.fStatDepth) << ", "
//                     << FloatToStr(itr->stGF.stNoSvRangeRight.fAvgDepth) << endl;
    }
    (*m_pOfs) << "iAbandonNum: " << to_string(iAbandonNum) << endl;
    (*m_pOfs) << "PositiveNum: " << to_string(iGoodNum) << endl;
    (*m_pOfs) << endl;

    for(int i=0; i<dtMax+1; i++)
    {
        ofsDepthInfo << arryDTName[i] << ": " << to_string(arryDT[i]) << endl;
    }

    ofs.close();
    ofsDepthInfo.close();
    ofsAbandonDepthInfo.close();

    //Step 3: save features for negative samples
    for(int i=0; i<dtMax+1; i++)
        arryDT[i] = 0;

    ofs.open(("./Output/NegativeSamples_Chrom_" + m_strChromName + ".txt").c_str());
    ofsDepthInfo.open(("./Output/NegativeSamplesDepthInfo_Chrom_" + m_strChromName + ".txt").c_str());

    ofsAbandonDepthInfo.open(("./Output/AbandonNegativeSamplesDepthInfo_Chrom_" + m_strChromName + ".txt").c_str());

    (*m_pOfs) << "Negative Samples" << endl;
    iAbandonNum = 0;
    iGoodNum = 0;
    for(vector<St_DRGroup>::iterator itr = vNegativeGroup.begin(); itr != vNegativeGroup.end(); itr++, iIndex++)
    {
        if(itr->stGF.stDCRangeCore.iStatDepthLen < 10)
        {
            iAbandonNum++;
            PrintCandiDepthDetail(ofsAbandonDepthInfo, itr, iIndex);
            continue;
        }

        arryDT[itr->stGF.stDCRangeCore.enDepthType]++;

        //iStart, iEnd, fPartialDepth, iDepthRangLen, fAvgRegDepthLeft, fAvgRegDepthRight, iStartClipReadsNum, iEndClipReadsNum, iDRReadsNum, iSpeciReadsNum
        //iStart and iEnd should turns to the length: abs(iEnd - iStart)
        /*
         * This is for positive samples:
         * So the real info should be like below:
         * 1: SV_start, SV_End, SV_Length
         * 2: Len (iEnd - iStart), fPartialDepth, iDepthRangLen, fAvgRegDepthLeft, fAvgRegDepthRight, iStartClipReadsNum, iEndClipReadsNum, iDRReadsNum, iSpeciReadsNum
         */
        ofs //<< to_string(itr->stSv.iPos) << "\t" SvClutering.py<< to_string(itr->stSv.iEnd) << "\t" << to_string(itr->stSv.iEnd - itr->stSv.iPos) << " --- "
            << to_string(itr->stGF.iEnd - itr->stGF.iStart) << "\t"
            << FloatToStr(itr->stGF.stDCRangeCore.fStatDepth) << "\t"
            << to_string(itr->stGF.stDCRangeCore.iStatDepthLen) << "\t"

            << FloatToStr(itr->stGF.stDCRangeLeft.fStatDepth) << "\t"
            << to_string(itr->stGF.stDCRangeLeft.iStatDepthLen) << "\t"

            << FloatToStr(itr->stGF.stDCRangeRight.fStatDepth) << "\t"
            << to_string(itr->stGF.stDCRangeRight.iStatDepthLen) << "\t"

            << FloatToStr(itr->stGF.stNoSvRangeLeft.fAvgDepth) << "\t"
            << FloatToStr(itr->stGF.stNoSvRangeRight.fAvgDepth) << "\t"

            << to_string(itr->stGF.iStartClipReadsNum) << "\t"
            << to_string(itr->stGF.iEndClipReadsNum) << "\t"
            << to_string(itr->stGF.iDRReadsNum) << "\t"
            << to_string(itr->stGF.iSpeciReadsNum) << endl;

        (*m_pOfs) << to_string(iIndex) << " -- \t"
             << FloatToStr(itr->stGF.stNoSvRangeLeft.fAvgDepth) << ", "
             << FloatToStr(itr->stGF.stDCRangeLeft.fStatDepth) << ", "
             << FloatToStr(itr->stGF.stDCRangeCore.fStatDepth) << ", "
             << FloatToStr(itr->stGF.stDCRangeRight.fStatDepth) << ", "
             << FloatToStr(itr->stGF.stNoSvRangeRight.fAvgDepth) << endl;

//        ofsDepthInfo << FloatToStr(itr->stGF.stNoSvRangeLeft.fAvgDepth) << ", "
//                     << FloatToStr(itr->stGF.stDCRangeLeft.fStatDepth) << ", "
//                     << FloatToStr(itr->stGF.stDCRangeCore.fStatDepth) << ", "
//                     << FloatToStr(itr->stGF.stDCRangeRight.fStatDepth) << ", "
//                     << FloatToStr(itr->stGF.stNoSvRangeRight.fAvgDepth) << endl;

        PrintCandiDepthDetail(ofsDepthInfo, itr, iIndex);
        iGoodNum++;
//        //For No Sv Left
//        ofsDepthInfo << itr->stGF.stNoSvRangeLeft.GetAvg() << " --- "
//                     << to_string(itr->GetAccuLen()) << " --- "
//                     << arryDTName[itr->stGF.stDCRangeCore.enDepthType] << " --- "
//                     << "(" << arryDTName[itr->stGF.stDCRangeLeft.enDepthType] << ", "
//                            << arryDTName[itr->stGF.stDCRangeRight.enDepthType] << ")" << endl;

//        //For DC Sv Left
//        ofsDepthInfo << itr->stGF.stDCRangeLeft.GetZeroPure() << ",  "
//                     << itr->stGF.stDCRangeLeft.GetZeroWeek() << ",  "
//                     << itr->stGF.stDCRangeLeft.GetZeroMix() << ",  "
//                     << itr->stGF.stDCRangeLeft.GetRepWeek() << ",  "
//                     << itr->stGF.stDCRangeLeft.GetRepStrong() << ",  "
//                     << itr->stGF.stDCRangeLeft.GetStat() << " ---  " << endl;

//        //For DC sv Core
//        ofsDepthInfo << itr->stGF.stDCRangeCore.GetZeroPure() << ",  "
//                     << itr->stGF.stDCRangeCore.GetZeroWeek() << ",  "
//                     << itr->stGF.stDCRangeCore.GetZeroMix() << ",  "
//                     << itr->stGF.stDCRangeCore.GetRepWeek() << ",  "
//                     << itr->stGF.stDCRangeCore.GetRepStrong() << ",  "
//                     << itr->stGF.stDCRangeCore.GetStat() << " ---  " << endl;

//        //For DC sv Right
//        ofsDepthInfo << itr->stGF.stDCRangeRight.GetZeroPure() << ",  "
//                     << itr->stGF.stDCRangeRight.GetZeroWeek() << ",  "
//                     << itr->stGF.stDCRangeRight.GetZeroMix() << ",  "
//                     << itr->stGF.stDCRangeRight.GetRepWeek() << ",  "
//                     << itr->stGF.stDCRangeRight.GetRepStrong() << ",  "
//                     << itr->stGF.stDCRangeRight.GetStat() << " ---  " << endl;

//        //For No sv Right
//        ofsDepthInfo << itr->stGF.stNoSvRangeRight.GetAvg() << endl << endl;
    }    

    for(int i=0; i<dtMax+1; i++)
    {
        ofsDepthInfo << arryDTName[i] << ": " << to_string(arryDT[i]) << endl;
    }

    ofs.close();
    (*m_pOfs) << "iAbandonNum: " << to_string(iAbandonNum) << endl;
    (*m_pOfs) << "NegativeNum: " << to_string(iGoodNum) << endl;
    (*m_pOfs) << endl;

////**********************************************************************************************
//    //Now we try to collect features for each group -->
//    for(vector<St_DRGroup>::iterator itr = vDRGroup.begin(); itr != vDRGroup.end(); itr++)
//    {
//        cout << endl << "---" << "(" << to_string(itr->iStart) << ", " << to_string(itr->iEnd) << ")" << endl;

//        itr->stGF.iStart = itr->iStart;
//        itr->stGF.iEnd = itr->iEnd;
//        itr->stGF.iSpeciReadsNum = itr->vRegReads.size();
//        itr->stGF.iDRReadsNum = itr->vDR.size();
//        //Try to get abnormal rang of depth --> Let's do it tomorrow
//        GetDepthFeatures(*itr, strBamFile);
//        cout << "\t GetDepthFeatures" << endl;

//        //Get Clip features -->
//        GetClipFeatures(*itr);
//        cout << "\t GetClipFeatures" << endl;
//    }
//    //<--

//    //Next Step --> make some additional filtering --> Go!!
//    cout << endl << "***************** Temp Filter by Features *******************" << endl;
//    cout << "Number of Groups: " << to_string(vDRGroup.size()) << endl;
//    for(vector<St_DRGroup>::iterator itr = vDRGroup.end() - 1; itr >= vDRGroup.begin(); itr--)
//    {
//        if(itr->stGF.fPartialDepth <= m_fAvgDepth * .3 ||
//           itr->stGF.fPartialDepth >= m_fAvgDepth * 1.5)
//        {}
//        else
//            vDRGroup.erase(itr);
//    }
//    cout << "\t After feature Filter: "  << to_string(vDRGroup.size()) << endl;
}

void ClsLearning::GetPositiveAndNegativeGroup(vector<St_DRGroup>& vDRGroup, vector<St_SV>& vSvDEL, //vector<St_SCRGroup>& vSCRGroup,
                                              vector<St_DRGroup>& vPositiveGroup, vector<St_DRGroup>& vNegativeGroup)
{
    unsigned int iStartDiff = 500;
    vPositiveGroup.clear();
    vNegativeGroup.clear();
    for(vector<St_DRGroup>::iterator itr = vDRGroup.begin(); itr != vDRGroup.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_SV>::iterator subItr = vSvDEL.begin(); subItr != vSvDEL.end(); subItr++)
        {
            if(subItr->iPos >=itr->iLeftBoundary && subItr->iPos <= itr->iRightBoundary)
            {
                if(abs(int(subItr->iPos - itr->iStart)) < iStartDiff ||
                   abs(int(subItr->iEnd - itr->iEnd) < iStartDiff))
                {
                    bFind = true;
                    itr->stSv = *subItr;
                    break;
                }
            }
        }
        if(bFind)
        {
            itr->enType = gtPositive;
            vPositiveGroup.push_back(*itr);
        }
        else
        {
            itr->enType = gtNegative;
            vNegativeGroup.push_back(*itr);
        }
    }
}

void ClsLearning::PrintCandiDepthDetail(ofstream& ofs, vector<St_DRGroup>::iterator itr, int iIndex)
{
    string arryDTName[dtMax+1] = {"ZeroPure", "ZeroWeek", "ZeroMix", "RepWeek", "RepStrong", "dtMax"};

    //For No Sv Left
    ofs << "<" << to_string(iIndex) << ">" << endl;
    ofs << itr->stGF.stNoSvRangeLeft.GetAvg() << " --- "
         << to_string(itr->GetAccuLen()) << " --- "
         << arryDTName[itr->stGF.stDCRangeCore.enDepthType] << " --- "
         << "(" << arryDTName[itr->stGF.stDCRangeLeft.enDepthType] << ", "
                << arryDTName[itr->stGF.stDCRangeRight.enDepthType] << ")" << endl;
    //For DC Sv Left
    ofs << itr->stGF.stDCRangeLeft.GetZeroPure() << ",  "
         << itr->stGF.stDCRangeLeft.GetZeroWeek() << ",  "
         << itr->stGF.stDCRangeLeft.GetZeroMix() << ",  "
         << itr->stGF.stDCRangeLeft.GetRepWeek() << ",  "
         << itr->stGF.stDCRangeLeft.GetRepStrong() << ",  "
         << itr->stGF.stDCRangeLeft.GetStat() << " ---  " << endl;

    //For DC sv Core
    ofs << itr->stGF.stDCRangeCore.GetZeroPure() << ",  "
         << itr->stGF.stDCRangeCore.GetZeroWeek() << ",  "
         << itr->stGF.stDCRangeCore.GetZeroMix() << ",  "
         << itr->stGF.stDCRangeCore.GetRepWeek() << ",  "
         << itr->stGF.stDCRangeCore.GetRepStrong() << ",  "
         << itr->stGF.stDCRangeCore.GetStat() << " ---  "
         << GetRatio(itr->stGF.stDCRangeCore.fStatLenRatio) << endl;

    //For DC sv Right
    ofs << itr->stGF.stDCRangeRight.GetZeroPure() << ",  "
         << itr->stGF.stDCRangeRight.GetZeroWeek() << ",  "
         << itr->stGF.stDCRangeRight.GetZeroMix() << ",  "
         << itr->stGF.stDCRangeRight.GetRepWeek() << ",  "
         << itr->stGF.stDCRangeRight.GetRepStrong() << ",  "
         << itr->stGF.stDCRangeRight.GetStat() << " ---  " << endl;

    //For No sv Right
    ofs << itr->stGF.stNoSvRangeRight.GetAvg() << endl << endl;
}

//This is for clustering --> Go!!
void ClsLearning::PrintClusterElements(vector<St_DRGroup>& vDRGroup, vector<St_SV>& vSvDEL)
{
    vector<St_DRGroup> vPositiveGroup;
    vector<St_DRGroup> vNegativeGroup;
    GetPositiveAndNegativeGroup(vDRGroup, vSvDEL, vPositiveGroup, vNegativeGroup);
    sort(vPositiveGroup.begin(), vPositiveGroup.end(), sort_positive_group_by_sv_small2large);

    ofstream ofs;
    ofs.open(("./Output/CluterPureZero_Chrom_" + m_strChromName + ".txt").c_str()); // 首先要normalization -->
    ofstream ofsLabel;
    ofsLabel.open(("./Output/CluterLabel_Chrom_" + m_strChromName + ".txt").c_str());
//    //Negative First, Positive Second
//    for(vector<St_DRGroup>::iterator itr = vNegativeGroup.begin(); itr != vNegativeGroup.end(); itr++)
//    {
////        if(itr->stGF.stDCRangeCore.iStatDepthLen < 10)
////            continue;

//        string strElement = GetElementValue(itr->stGF);
//        ofs << strElement << endl;

//        ofsLabel << "0, "; // << endl;
//    }
//    cout << "Negative Number: " << to_string(vNegativeGroup.size()) << endl;


//    string strPosLable = "";
//    for(vector<St_DRGroup>::iterator itr = vPositiveGroup.begin(); itr != vPositiveGroup.end(); itr++)
//    {
////        if(itr->stGF.stDCRangeCore.iStatDepthLen < 10)
////            continue;

//        string strElement = GetElementValue(itr->stGF);
//        ofs << strElement << endl;

//        strPosLable += "1, ";// << endl;
//    }
//    strPosLable = strPosLable.substr(0, strPosLable.length() - 2);
//    ofsLabel << strPosLable;

//// **********************************************

    //Negative First, Positive Second
    int iEligableSample = 0;
    for(vector<St_DRGroup>::iterator itr = vPositiveGroup.begin(); itr != vPositiveGroup.end(); itr++)
    {
        if(itr->stGF.stDCRangeCore.iStatDepthLen < 10)
        {}
        else
            iEligableSample++;
        //    continue;

        string strElement = GetElementValue(itr->stGF);
        ofs << strElement << endl;

        ofsLabel << "1, "; // << endl;
    }
    (*m_pOfs) << "Positive Number: " << to_string(iEligableSample) << endl;


    string strNegLable = "";
    iEligableSample = 0;
    for(vector<St_DRGroup>::iterator itr = vNegativeGroup.begin(); itr != vNegativeGroup.end(); itr++)
    {
        if(itr->stGF.stDCRangeCore.iStatDepthLen < 10)
        {}
        else
            iEligableSample++;
        //    continue;

        string strElement = GetElementValue(itr->stGF);
        ofs << strElement << endl;

        strNegLable += "0, ";// << endl;
    }
    strNegLable = strNegLable.substr(0, strNegLable.length() - 2);
    ofsLabel << strNegLable;


    //ofsLabel << endl;
    (*m_pOfs) << "Negative Number: " << to_string(iEligableSample) << endl;

    ofs.close();
    ofsLabel.close();
}

int ClsLearning::GetNormalizedValue(int iRealLen, int iSegLen, int iStdNormLen)
{
    float fRatio = (float)iRealLen / iSegLen;
    int iValue = iStdNormLen * fRatio;
    return iValue;
}

string ClsLearning::GetElementValue(St_GroupFeatures& stGF)
{
    string strElement = "";
//    //This is for the length based on Zero Pure for all five ranges only.
//    strElement = to_string(GetNormalizedValue(stGF.stNoSvRangeLeft.iZeroPureDepthLen, stGF.stNoSvRangeLeft.iSegLen)) + ", " +
//                 to_string(GetNormalizedValue(stGF.stDCRangeLeft.iZeroPureDepthLen, stGF.stDCRangeLeft.iSegLen)) + ", " +
//                 to_string(GetNormalizedValue(stGF.stDCRangeCore.iZeroPureDepthLen, stGF.stDCRangeCore.iSegLen)) + ", " +
//                 to_string(GetNormalizedValue(stGF.stDCRangeRight.iZeroPureDepthLen, stGF.stDCRangeRight.iSegLen)) + ", " +
//                 to_string(GetNormalizedValue(stGF.stNoSvRangeRight.iZeroPureDepthLen, stGF.stNoSvRangeRight.iSegLen));

//    //------->
//    //DC Range Left
//    strElement += to_string(GetNormalizedValue(stGF.stDCRangeLeft.iZeroPureDepthLen, stGF.stDCRangeLeft.iSegLen)) + ", " +
//                  to_string(GetNormalizedValue(stGF.stDCRangeLeft.iZeroWeekDepthLen, stGF.stDCRangeLeft.iSegLen)) + ", " +
//                  to_string(GetNormalizedValue(stGF.stDCRangeLeft.iZeroMixDepthLen, stGF.stDCRangeLeft.iSegLen)) + ", " +
//                  to_string(GetNormalizedValue(stGF.stDCRangeLeft.iRepWeekDepthLen, stGF.stDCRangeLeft.iSegLen)) + ", " +
//                  to_string(GetNormalizedValue(stGF.stDCRangeLeft.iRepStrongDepthLen, stGF.stDCRangeLeft.iSegLen)) + ", ";

//    //DC Range Core
//    strElement += to_string(GetNormalizedValue(stGF.stDCRangeCore.iZeroPureDepthLen, stGF.stDCRangeCore.iSegLen)) + ", " +
//                  to_string(GetNormalizedValue(stGF.stDCRangeCore.iZeroWeekDepthLen, stGF.stDCRangeCore.iSegLen)) + ", " +
//                  to_string(GetNormalizedValue(stGF.stDCRangeCore.iZeroMixDepthLen, stGF.stDCRangeCore.iSegLen)) + ", " +
//                  to_string(GetNormalizedValue(stGF.stDCRangeCore.iRepWeekDepthLen, stGF.stDCRangeCore.iSegLen)) + ", " +
//                  to_string(GetNormalizedValue(stGF.stDCRangeCore.iRepStrongDepthLen, stGF.stDCRangeCore.iSegLen)) + ", ";

//    //DC Range Right
//    strElement += to_string(GetNormalizedValue(stGF.stDCRangeRight.iZeroPureDepthLen, stGF.stDCRangeRight.iSegLen)) + ", " +
//                  to_string(GetNormalizedValue(stGF.stDCRangeRight.iZeroWeekDepthLen, stGF.stDCRangeRight.iSegLen)) + ", " +
//                  to_string(GetNormalizedValue(stGF.stDCRangeRight.iZeroMixDepthLen, stGF.stDCRangeRight.iSegLen)) + ", " +
//                  to_string(GetNormalizedValue(stGF.stDCRangeRight.iRepWeekDepthLen, stGF.stDCRangeRight.iSegLen)) + ", " +
//                  to_string(GetNormalizedValue(stGF.stDCRangeRight.iRepStrongDepthLen, stGF.stDCRangeRight.iSegLen));
//    //<------

//    //-> This is try 2
    strElement += to_string(GetNormalizedValue(stGF.stDCRangeCore.iZeroPureDepthLen, stGF.stDCRangeCore.iSegLen)) + ", " +
                  to_string(GetNormalizedValue(stGF.stDCRangeCore.iZeroWeekDepthLen, stGF.stDCRangeCore.iSegLen)) + ", " +
                  to_string(GetNormalizedValue(stGF.stDCRangeCore.iZeroMixDepthLen, stGF.stDCRangeCore.iSegLen)) + ", " +
                  to_string(GetNormalizedValue(stGF.stDCRangeCore.GetRegDepthLen(), stGF.stDCRangeCore.iSegLen));
//                  //to_string(GetNormalizedValue(stGF.stDCRangeCore.iStatDepthLen, stGF.stDCRangeCore.iSegLen));
//                  to_string(GetNormalizedValue(stGF.stDCRangeCore.iRepWeekDepthLen, stGF.stDCRangeCore.iSegLen)) + ", " +
//                  to_string(GetNormalizedValue(stGF.stDCRangeCore.iRepStrongDepthLen, stGF.stDCRangeCore.iSegLen));
//    //<--

//    //-->This is try 3
//    strElement += FloatToStr(stGF.stNoSvRangeLeft.fAvgDepth) + ", " +
//                  FloatToStr(stGF.stDCRangeLeft.fStatDepth) + ", " +
//                  FloatToStr(stGF.stDCRangeCore.fStatDepth) + ", " +
//                  FloatToStr(stGF.stDCRangeRight.fStatDepth) + ", " +
//                  FloatToStr(stGF.stNoSvRangeRight.fAvgDepth);
//    //<--

//    //Try 4
//    strElement = GetNormDepth(stGF);

    return strElement;
}

string ClsLearning::GetNormDepth(St_GroupFeatures& stGF)
{
    float fStd = m_fAvgDepth * 2;
    vector<float> vValue = {stGF.stNoSvRangeLeft.fAvgDepth, stGF.stDCRangeLeft.fStatDepth,
                            stGF.stDCRangeCore.fStatDepth, stGF.stDCRangeRight.fStatDepth,
                            stGF.stNoSvRangeRight.fAvgDepth};
    sort(vValue.begin(), vValue.end());
    float fMax = *(vValue.end() - 1);
    float fScale = (float)fStd / fMax;

    string strDepth = "";

    vector<float> vRawValue = {stGF.stNoSvRangeLeft.fAvgDepth, stGF.stDCRangeLeft.fStatDepth,
                               stGF.stDCRangeCore.fStatDepth, stGF.stDCRangeRight.fStatDepth,
                               stGF.stNoSvRangeRight.fAvgDepth};
    for(vector<float>::iterator itr = vRawValue.begin(); itr != vRawValue.end(); itr++)
    {
        strDepth += FloatToStr((*itr) * fScale) + ", ";
    }
    strDepth = strDepth.substr(0, strDepth.length() - 2);
    return strDepth;
}

void ClsLearning::ClusterByPython(vector<St_DRGroup>& vDRGroup, vector<St_SV>& vSvDEL,
                                  string strPythonCode, St_FinalResult& stFinalResult)
{
    vector<St_DRGroup> vPositiveGroup;
    vector<St_DRGroup> vNegativeGroup;

    //For Positive and Negative
    GetPositiveAndNegativeGroup(vDRGroup, vSvDEL, vPositiveGroup, vNegativeGroup);
    sort(vPositiveGroup.begin(), vPositiveGroup.end(), sort_positive_group_by_sv_small2large);

    //Save the input file
    //We have two files here:
    //(1) Samples (mandate)
    //(2) Lables: this is only for testing  right now, we will remove it in final version
    ofstream ofs;
    string strSampleList = "./Output/CluterPureZero_Chrom_" + m_strChromName + ".txt";
    ofs.open(strSampleList.c_str()); // 首先要normalization -->
    ofstream ofsLabel;
    string strLabelList = "./Output/CluterLabel_Chrom_" + m_strChromName + ".txt";
    ofsLabel.open(strLabelList.c_str());

    string strLable = "";
    for(vector<St_DRGroup>::iterator itr = vPositiveGroup.begin(); itr != vPositiveGroup.end(); itr++)
    {
        string strElement = GetElementValue(itr->stGF);
        ofs << strElement << endl;

        strLable += "1, ";
    }
    if(vNegativeGroup.empty())
        strLable = strLable.substr(0, strLable.length() - 2);
    ofsLabel << strLable;

    strLable = "";
    for(vector<St_DRGroup>::iterator itr = vNegativeGroup.begin(); itr != vNegativeGroup.end(); itr++)
    {
        string strElement = GetElementValue(itr->stGF);
        ofs << strElement << endl;

        strLable += "0, ";
    }
    if(!vNegativeGroup.empty())
        strLable = strLable.substr(0, strLable.length() - 2);
    ofsLabel << strLable;

    ofs.close();
    ofsLabel.close();

    //Run Cluster through Python
    string strAgu1 = strSampleList;
    string strAgu2 = strLabelList;
    string strAgu3 = m_strChromName;
    string strAgu4 = m_strSampleName;
    string strPyMain = strPythonCode; //"/home/xin/lxwg/WorkStudio/Prototype/SV_Draft/SvDelPainter/Learning/SvLearning/SvClutering.py";
    string strCmd = "python3 " + strPyMain + " "
                    + strAgu1 +  " "
                    + strAgu2 + " "
                    + strAgu3 + " "
                    + strAgu4;
    system(strCmd.c_str());
    (*m_pOfs) << "Python finished!" << endl;

    //Read Label File generated by python --> go 完成这个十二点半睡觉
    vector<int> vLabel;
    string strLabel = "./Pics/ClusterLables_Chrom_" + m_strChromName + ".txt";
    ifstream ifs;
    ifs.open(strLabel.c_str());
    string strLine = "";
    getline(ifs, strLine);
    int iGroupNum = 0;
    for(int i=0; i<strLine.length(); i++)
    {
        int iGroupIndex = atoi(strLine.substr(i, 1).c_str());
        vLabel.push_back(iGroupIndex);
       if(iGroupIndex > iGroupNum)
           iGroupNum = iGroupIndex;
    }
    if(iGroupNum > 0)
        iGroupNum++;
    (*m_pOfs) << to_string(vLabel.size()) << endl;
    ifs.close();

    //Choose correct Group -->
    if(iGroupNum == 0)
    {
        (*m_pOfs) << "Bad iGroupNum!" << endl;
        return;
    }
    //put sample in same group;
    vector<St_GroupCluster> vGroupCluster;
    vGroupCluster.resize(iGroupNum);

    vector<St_DRGroup> vTmpGroup;
    vTmpGroup.insert(vTmpGroup.end(), vPositiveGroup.begin(), vPositiveGroup.end());
    vTmpGroup.insert(vTmpGroup.end(), vNegativeGroup.begin(), vNegativeGroup.end());
    int i = 0;
    for(vector<St_DRGroup>::iterator itr = vTmpGroup.begin(); itr != vTmpGroup.end(); itr++, i++)
    {
        vGroupCluster[vLabel[i]].vGroup.push_back(*itr);
    }

    //Get some good back
    i = 0;
    for(vector<St_GroupCluster>::iterator itr = vGroupCluster.begin();
        itr != vGroupCluster.end(); itr++, i++)
    {
        //Calculate some key value to check if we belongs to bad sample -> Go tomorrow
        (*m_pOfs) << "Group " << to_string(i) << ": " << to_string(itr->vGroup.size()) << endl << " === " << endl;
        CalcFeatureOfEachCluster(*itr);
        DoFilterOfEachCluster(*itr, i);
        (*m_pOfs) << endl << "*******************" << endl << endl;
    }

    //Check How many positive and how many negative
    //Output the final number of positive and negative
    int iPositiveNum = 0;
    int iNegativeNum = 0;
    for(vector<St_GroupCluster>::iterator itr = vGroupCluster.begin();
        itr != vGroupCluster.end(); itr++)
    {
        iPositiveNum += itr->iPositive;
        iNegativeNum += itr->iNegative;
    }

    (*m_pOfs) << endl << "Final Positive: " << to_string(iPositiveNum) << endl;
    (*m_pOfs) << endl << "Final Negative: " << to_string(iNegativeNum) << endl;

    stFinalResult.iPositiveNum = iPositiveNum;
    stFinalResult.iNegativeNum = iNegativeNum;
    stFinalResult.iStdNum = vSvDEL.size();

    //Save Candidate
    string strResult = "./Output/Result_" + this->m_strChromName + ".txt";
    ofs.open(strResult.c_str());
    for(vector<St_GroupCluster>::iterator itr = vGroupCluster.begin();
        itr != vGroupCluster.end(); itr++)
    {
        for(vector<St_DRGroup>::iterator itrGroup=itr->vGroup.begin();
            itrGroup!=itr->vGroup.end(); itrGroup++)
        {
            ofs << this->m_strChromName << "\t"
                << to_string(itrGroup->GetStart()) << "\t"
                << to_string(itrGroup->iEnd) << "\t"
                << to_string(itrGroup->GetAccuLen()) << endl;
        }
    }
    ofs.close();
}

enum En_BestLenType{bltZeroPure, bltZeroWeek, bltZeroMix, bltMax};
void ClsLearning::DoFilterOfEachCluster(St_GroupCluster& stCluster, int iClusterIndex)
{
    //Confirm which one is bad cluster or good cluster -> Go!! -->
    bool bGoodCluster = false;
    vector<float> vRatio {stCluster.fAvgZeroPureLenHalfBestRatio, stCluster.fAvgZeroWeekLenHalfBestRatio, stCluster.fAvgZeroMixLenHalfBestRatio};
    sort(vRatio.begin(), vRatio.end());
    float f2SmallValueSum = vRatio[0] + vRatio[1];

    int iGoodLowNum = 0;
    int iGoodHighNum = 0;    

    if(stCluster.fAvgZeroMixLenHalfBestRatio > .5)
        iGoodHighNum++;
    else if(stCluster.fAvgZeroMixLenHalfBestRatio > .1)
        iGoodLowNum++;

    if(stCluster.fAvgZeroPureLenHalfBestRatio > .5)
        iGoodHighNum++;
    else if(stCluster.fAvgZeroPureLenHalfBestRatio > .1)
        iGoodLowNum++;

    if(stCluster.fAvgZeroWeekLenHalfBestRatio > .5)
        iGoodHighNum++;
    else if(stCluster.fAvgZeroWeekLenHalfBestRatio > .1)
        iGoodLowNum++;

    if(iGoodHighNum > 0)
    {
        if(*(vRatio.end()-1) > .7)  // It's good, 因为vRatio 是排序过的，因此最后一个就是最大值
            bGoodCluster = true;
        else if(f2SmallValueSum > .2)
            bGoodCluster = true;
    }
    else if(iGoodHighNum + iGoodLowNum == 3)
        bGoodCluster = true;
    else if(*(vRatio.end()-1) > (.5 * .9))  // we can change some --> many be all together > .7
    {
        if(vRatio[0] + vRatio[1] > .2)
            bGoodCluster = true;
    }
    //<--

    if(bGoodCluster)
    {
        //this is good cluster
        FilterGoodCluster(stCluster, iClusterIndex);
    }
    else
    {
        //This is bad cluster
        FilterBadCluster(stCluster, iClusterIndex);
    }
}

//Purpose: Kick some bad one out
void ClsLearning::FilterGoodCluster(St_GroupCluster& stCluster, int iClusterIndex)
{
    //****** kick the bad out ******
    float arryAvgHalfBestRatio[3] = {stCluster.fAvgZeroPureLenHalfBestRatio, stCluster.fAvgZeroWeekLenHalfBestRatio, stCluster.fAvgZeroMixLenHalfBestRatio};
    float arryAvgRatio[3] = {stCluster.fAvgZeroPureLenRatio, stCluster.fAvgZeroWeekLenRatio, stCluster.fAvgZeroMixLenRatio};

    //Know which type of half best value is the largest one
    En_BestLenType enBLT = bltZeroMix;
    if(arryAvgHalfBestRatio[enBLT] <= arryAvgHalfBestRatio[bltZeroWeek])
    {
        enBLT = bltZeroWeek;
        if(arryAvgHalfBestRatio[enBLT] <= arryAvgHalfBestRatio[bltZeroPure])
            enBLT = bltZeroPure;
    }
    else
    {
        if(arryAvgHalfBestRatio[enBLT] <= arryAvgHalfBestRatio[bltZeroPure])
            enBLT = bltZeroPure;
    }

    int iOrgSize = stCluster.vGroup.size();
    int iPositiveNum = 0;
    for(vector<St_DRGroup>::iterator itr = stCluster.vGroup.end() - 1; itr >= stCluster.vGroup.begin(); itr--)
    {
        bool bBad = false;
        switch(enBLT)
        {
            case bltZeroPure:
            {
                if(itr->stGF.stDCRangeCore.GetZeroPureDepthLenRatio() < stCluster.fAvgZeroPureLenHalfBestRatio * .5)
                    bBad = true;

                if(bBad)
                    stCluster.vGroup.erase(itr);

                break;
            }
            case bltZeroWeek:
            {
                if(itr->stGF.stDCRangeCore.GetZeroWeekDepthLenRatio() < stCluster.fAvgZeroWeekLenHalfBestRatio * .5)
                    bBad = true;

                if(bBad)
                    stCluster.vGroup.erase(itr);
                break;
            }
            case bltZeroMix:
            {
                if(itr->stGF.stDCRangeCore.GetZeroMixDepthLenRatio() < stCluster.fAvgZeroMixLenHalfBestRatio * .5)
                {
                    if(itr->stGF.stDCRangeCore.GetZeroPureDepthLenRatio() +
                       itr->stGF.stDCRangeCore.GetZeroWeekDepthLenRatio() > stCluster.fAvgZeroMixLenHalfBestRatio)
                    {}
                    else
                        bBad = true;
                }

                if(bBad)
                    stCluster.vGroup.erase(itr);

                break;
            }
            default:
                break;
        }

        //Coount Positive
        if(!bBad)
        {
            if(itr->enType == gtPositive)
                iPositiveNum++;
        }

    }

    stCluster.iPositive = iPositiveNum;
    stCluster.iNegative = stCluster.vGroup.size() - iPositiveNum;

    (*m_pOfs) << "Cluster_Good_Filter_Before: " << to_string(iOrgSize) << endl;
    (*m_pOfs) << "Cluster_Good_Filter_After : " << to_string(stCluster.vGroup.size()) << endl;
    (*m_pOfs) << "\t" << "Positive: " << to_string(iPositiveNum) << endl;
    (*m_pOfs) << "\t" << "Negative: " << to_string(stCluster.vGroup.size() - iPositiveNum) << endl;

//    //Now we get the max best half of enBLT --> Go!!
//    //Pick the good out
//    ofstream ofs;
//    ofs.open(("cluster_good_" + to_string(iClusterIndex)).c_str());
//    ofstream ofsLabel;
//    ofsLabel.open(("cluster_good_" + to_string(iClusterIndex) + "_label").c_str());

//    //--> Output Cluster -> Go
//    string strLable = "";
//    string arryLabel[2] = {"0", "1"};
//    for(vector<St_DRGroup>::iterator itr = stCluster.vGroup.begin(); itr != stCluster.vGroup.end(); itr++)
//    {
//        string strElement = GetElementValue(itr->stGF);
//        ofs << strElement << endl;
//        strLable += arryLabel[itr->enType] + ", ";// << endl;
//    }
//    if(!stCluster.vGroup.empty())
//        strLable = strLable.substr(0, strLable.length() - 2);
//    //<--
//    ofsLabel << strLable;

//    ofs.close();
//    ofsLabel.close();
}

//Purpose: get some good one back
void ClsLearning::FilterBadCluster(St_GroupCluster& stCluster, int iClusterIndex)
{
    //--> Get Max Type
    float arryAvgHalfBestRatio[3] = {stCluster.fAvgZeroPureLenHalfBestRatio, stCluster.fAvgZeroWeekLenHalfBestRatio, stCluster.fAvgZeroMixLenHalfBestRatio};
    float arryAvgRatio[3] = {stCluster.fAvgZeroPureLenRatio, stCluster.fAvgZeroWeekLenRatio, stCluster.fAvgZeroMixLenRatio};

    En_BestLenType enBLT = bltZeroMix;

    if(arryAvgHalfBestRatio[enBLT] <= arryAvgHalfBestRatio[bltZeroWeek])
    {
        enBLT = bltZeroWeek;        
        if(arryAvgHalfBestRatio[enBLT] <= arryAvgHalfBestRatio[bltZeroPure])        
            enBLT = bltZeroPure;
    }
    else
    {
        if(arryAvgHalfBestRatio[enBLT] <= arryAvgHalfBestRatio[bltZeroPure])
            enBLT = bltZeroPure;
    }
    //<--
    int iOrgSize = stCluster.vGroup.size();
    int iPositiveNum = 0;
    for(vector<St_DRGroup>::iterator itr = stCluster.vGroup.end() - 1; itr >= stCluster.vGroup.begin(); itr--)
    {
        bool bBad = true;
        switch(enBLT) //Best
        {
            case bltZeroPure:
            {                            
                //cout << "bltZeroPure";
                //if(itr->stGF.stDCRangeCore.GetZeroPureDepthLenRatio() > stCluster.fAvgZeroPureLenHalfBestRatio)
                if(itr->stGF.stDCRangeCore.GetZeroPureDepthLenRatio() > stCluster.fAvgZeroPureLenHalfBestSTDPlus)
                    bBad = false;

                if(bBad)
                    stCluster.vGroup.erase(itr);
                break;
            }
            case bltZeroWeek:
            {
                //cout << "bltZeroWeek";
                //if(itr->stGF.stDCRangeCore.GetZeroWeekDepthLenRatio() > stCluster.fAvgZeroWeekLenHalfBestRatio)
                if(itr->stGF.stDCRangeCore.GetZeroWeekDepthLenRatio() > stCluster.fAvgZeroWeekLenHalfBestSTDPlus)
                    bBad = false;

                if(bBad)
                    stCluster.vGroup.erase(itr);
                break;
            }
            case bltZeroMix:
            {
                /*if(stCluster.fAvgZeroMixLenHalfBestRatio >= .5)
                {
                    if(itr->stGF.stDCRangeCore.GetZeroMixDepthLenRatio() >= stCluster.fAvgZeroMixLenHalfBestSTDPlus)
                        bBad = false;
                }
                else */if(stCluster.fAvgZeroMixLenHalfBestRatio < .5 &&
                        stCluster.fAvgZeroPureLenHalfBestRatio < .1 &&
                        stCluster.fAvgZeroWeekLenHalfBestRatio < .1)
                {
                    //cout << "bltZeroMix_1_1";
                    //if(itr->stGF.stDCRangeCore.GetZeroMixDepthLenRatio() >= stCluster.fAvgZeroMixLenHalfBestRatio)
                    if(itr->stGF.stDCRangeCore.GetZeroMixDepthLenRatio() >= stCluster.fAvgZeroMixLenHalfBestSTDPlus)
                        bBad = false;
                    else
                    {
                        //cout << "bltZeroMix_1_2";
                        //if(itr->stGF.stDCRangeCore.GetZeroPureDepthLenRatio() + itr->stGF.stDCRangeCore.GetZeroWeekDepthLenRatio() > stCluster.fAvgZeroMixLenHalfBestRatio)
                        if(itr->stGF.stDCRangeCore.GetZeroPureDepthLenRatio() + itr->stGF.stDCRangeCore.GetZeroWeekDepthLenRatio() >= stCluster.fAvgZeroMixLenHalfBestSTDPlus)
                            bBad = false;
                    }
                }
                else if( stCluster.fAvgZeroMixLenHalfBestRatio < .5 &&
                        (stCluster.fAvgZeroPureLenHalfBestRatio > .1 || stCluster.fAvgZeroWeekLenHalfBestRatio > .1))
                {
                    //cout << "bltZeroMix_2_1";
                    if(itr->stGF.stDCRangeCore.GetZeroMixDepthLenRatio() >= stCluster.fAvgZeroMixLenHalfBestRatio)
                        bBad = false;
                    else
                    {
//                        //cout << "bltZeroMix_2_2";
//                        float fCurSecondLarge = itr->stGF.stDCRangeCore.GetZeroPureDepthLenRatio();
//                        float fSecondLarge = stCluster.fAvgZeroPureLenHalfBestRatio;
//                        if(fSecondLarge < stCluster.fAvgZeroWeekLenHalfBestRatio)
//                        {
//                            fSecondLarge = stCluster.fAvgZeroWeekLenHalfBestRatio;
//                            fCurSecondLarge = itr->stGF.stDCRangeCore.GetZeroWeekDepthLenRatio();
//                        }

//                        if(itr->stGF.stDCRangeCore.GetZeroMixDepthLenRatio() >= stCluster.fAvgZeroMixLenRatio &&
//                           fCurSecondLarge > fSecondLarge)
//                            bBad = false;

                        if(itr->stGF.stDCRangeCore.GetZeroPureDepthLenRatio() + itr->stGF.stDCRangeCore.GetZeroWeekDepthLenRatio() >= stCluster.fAvgZeroMixLenHalfBestSTDPlus)
                            bBad = false;
                    }
                }
                else
                {}

                if(bBad)
                    stCluster.vGroup.erase(itr);
                break;
            }
            default:
                break;
        }
        //cout << endl;
        if(!bBad)
        {
            if(itr->enType == gtPositive)
                iPositiveNum++;
        }
    }

    stCluster.iPositive = iPositiveNum;
    stCluster.iNegative = stCluster.vGroup.size() - iPositiveNum;

    (*m_pOfs) << "Cluster_Bad_Filter_Before: " << to_string(iOrgSize) << endl;
    (*m_pOfs) << "Cluster_Bad_Filter_After : " << to_string(stCluster.vGroup.size()) << endl;
    (*m_pOfs) << "\t" << "Positive: " << to_string(iPositiveNum) << endl;
    (*m_pOfs) << "\t" << "Negative: " << to_string(stCluster.vGroup.size() - iPositiveNum) << endl;

//    //Pick the good out
//    ofstream ofs;
//    ofs.open(("cluster_bad_" + to_string(iClusterIndex)).c_str());
//    ofstream ofsLabel;
//    ofsLabel.open(("cluster_bad_" + to_string(iClusterIndex) + "_label").c_str());

//    //--> Output Cluster -> Go
//    string strLable = "";
//    string arryLabel[2] = {"0", "1"};
//    for(vector<St_DRGroup>::iterator itr = stCluster.vGroup.begin(); itr != stCluster.vGroup.end(); itr++)
//    {
//        string strElement = GetElementValue(itr->stGF);
//        ofs << strElement << endl;
//        strLable += arryLabel[itr->enType] + ", ";// << endl;
//    }
//    if(!stCluster.vGroup.empty())
//        strLable = strLable.substr(0, strLable.length() - 2);
//    //<--
//    ofsLabel << strLable;

//    ofs.close();
//    ofsLabel.close();
}

void ClsLearning::CalcFeatureOfEachCluster(St_GroupCluster& stCluster)
{
    //Go!!
    float fZeroPureNormStd = floor(.25 * m_fAvgDepth);;
    if(fZeroPureNormStd < 1)
        fZeroPureNormStd = 1;

    int iZeroMixDepthMin = floor(.25 * m_fAvgDepth);
    int iZeroMixDepthMax = ceil(.5 * m_fAvgDepth);
    float fRealZeroMixDepth = .5 * m_fAvgDepth;
    int iZeroMixDepthNormStd = iZeroMixDepthMax - iZeroMixDepthMin;

    vector<float> vZeroPureRatio;
    vector<float> vZeroWeekRatio;
    vector<float> vZeroWeekRatioNorm;

    vector<float> vZeroMixRatio;
    vector<float> vZeroMixRatioNorm;

    string arrySampleType[gtMax] = {"0", "1"};
    string strTypeNegSum = "";
    int iNegSum = 0;
    string strTypePositiveSum = "";
    int iPositiveSum = 0;
    for(vector<St_DRGroup>::iterator itr = stCluster.vGroup.begin(); itr != stCluster.vGroup.end(); itr++)
    {
        //For Zero range
        vZeroPureRatio.push_back(itr->stGF.stDCRangeCore.GetZeroPureDepthLenRatio());
        vZeroWeekRatio.push_back(itr->stGF.stDCRangeCore.GetZeroWeekDepthLenRatio());
        vZeroWeekRatioNorm.push_back(itr->stGF.stDCRangeCore.GetZeroWeekDepthLenRatio() *
                                     ((fZeroPureNormStd - itr->stGF.stDCRangeCore.fZeroWeekDepth)/fZeroPureNormStd));

        //For Zero Mix range
        vZeroMixRatio.push_back(itr->stGF.stDCRangeCore.GetZeroMixDepthLenRatio());
        vZeroMixRatioNorm.push_back(itr->stGF.stDCRangeCore.GetZeroMixDepthLenRatio() *
                                    (iZeroMixDepthNormStd - abs(fRealZeroMixDepth - itr->stGF.stDCRangeCore.fZeroWeekDepth)/iZeroMixDepthNormStd));

        if(itr->enType == gtNegative)
        {
            strTypeNegSum += FloatToStr(itr->stGF.stDCRangeCore.GetZeroPureDepthLenRatio()) + " " +
                             FloatToStr(itr->stGF.stDCRangeCore.GetZeroWeekDepthLenRatio()) + " " +
                             FloatToStr(itr->stGF.stDCRangeCore.GetZeroMixDepthLenRatio()) + " --- " +
                             FloatToStr(itr->stGF.stDCRangeCore.fRepWeekDepth) + "<" +
                                   FloatToStr(itr->stGF.stDCRangeCore.GetRepWeekDepthLenRatio()) + ">, " +
                             FloatToStr(itr->stGF.stDCRangeCore.fRepStrongDepth) + "<" +
                                   FloatToStr(itr->stGF.stDCRangeCore.GetRepStrongDepthLenRatio()) + ">";
            strTypeNegSum += " =====> " + itr->stGF.stDCRangeCore.GetZeroPure() + ",  " +
                                    itr->stGF.stDCRangeCore.GetZeroWeek() + ",  " +
                                    itr->stGF.stDCRangeCore.GetZeroMix() + ",  " +
                                    itr->stGF.stDCRangeCore.GetRepWeek() + ",  " +
                                    itr->stGF.stDCRangeCore.GetRepStrong() + ",  " +
                                    itr->stGF.stDCRangeCore.GetStat() + " ---  " +
                                    GetRatio(itr->stGF.stDCRangeCore.fStatLenRatio) + "\n";
            iNegSum++;
        }
        else if(itr->enType == gtPositive)
        {
            strTypePositiveSum += FloatToStr(itr->stGF.stDCRangeCore.GetZeroPureDepthLenRatio()) + " " +
                                  FloatToStr(itr->stGF.stDCRangeCore.GetZeroWeekDepthLenRatio()) + " " +
                                  FloatToStr(itr->stGF.stDCRangeCore.GetZeroMixDepthLenRatio()) + " --- " +
                                  FloatToStr(itr->stGF.stDCRangeCore.fRepWeekDepth) + "<" +
                                             FloatToStr(itr->stGF.stDCRangeCore.GetRepWeekDepthLenRatio()) + ">, " +
                                  FloatToStr(itr->stGF.stDCRangeCore.fRepStrongDepth) + "<" +
                                             FloatToStr(itr->stGF.stDCRangeCore.GetRepStrongDepthLenRatio()) + ">";
            //For DC sv Core
           strTypePositiveSum += " =====> " + itr->stGF.stDCRangeCore.GetZeroPure() + ",  " +
                                        itr->stGF.stDCRangeCore.GetZeroWeek() + ",  " +
                                        itr->stGF.stDCRangeCore.GetZeroMix() + ",  " +
                                        itr->stGF.stDCRangeCore.GetRepWeek() + ",  " +
                                        itr->stGF.stDCRangeCore.GetRepStrong() + ",  " +
                                        itr->stGF.stDCRangeCore.GetStat() + " ---  " +
                                        GetRatio(itr->stGF.stDCRangeCore.fStatLenRatio) + "\n";
           iPositiveSum++;
        }
    }

    strTypeNegSum += to_string(iNegSum) + "\n";
    strTypePositiveSum += to_string(iPositiveSum) + "\n";

    //Print out the length ratio
    (*m_pOfs) << "vZeroPureRatio: " << endl;
    PrintFVector(vZeroPureRatio, stCluster.fAvgZeroPureLenRatio, stCluster.fAvgZeroPureLenHalfBestRatio,
                 stCluster.fAvgZeroPureSTDPlus, stCluster.fAvgZeroPureSTDMinus,
                 stCluster.fAvgZeroPureLenHalfBestSTDPlus, stCluster.fAvgZeroPureLenHalfBestSTDMinus);
    (*m_pOfs) << " === " << endl;

    (*m_pOfs) << "vZeroWeekRatio: " << endl;
    PrintFVector(vZeroWeekRatio, stCluster.fAvgZeroWeekLenRatio, stCluster.fAvgZeroWeekLenHalfBestRatio,
                 stCluster.fAvgZeroWeekSTDPlus, stCluster.fAvgZeroWeekSTDMinus,
                 stCluster.fAvgZeroWeekLenHalfBestSTDPlus, stCluster.fAvgZeroWeekLenHalfBestSTDMinus);
    (*m_pOfs) << " === " << endl;

    (*m_pOfs) << "vZeroWeekRatioNorm: " << endl;
    PrintFVector(vZeroWeekRatioNorm, stCluster.fNromAvgZeroWeekLenRatio, stCluster.fNormAvgZeroWeekLenHalfBestRatio,
                 stCluster.fNormAvgZeroWeekSTDPlus, stCluster.fNormAvgZeroWeekSTDMinus,
                 stCluster.fNormAvgZeroWeekLenHalfBestSTDPlus, stCluster.fNormAvgZeroWeekLenHalfBestSTDMinus);
    (*m_pOfs) << " === " << endl;

    (*m_pOfs) << "vZeroMixRatio: " << endl;
    PrintFVector(vZeroMixRatio, stCluster.fAvgZeroMixLenRatio, stCluster.fAvgZeroMixLenHalfBestRatio,
                 stCluster.fAvgZeroMixSTDPlus, stCluster.fAvgZeroMixSTDMinus,
                 stCluster.fAvgZeroMixLenHalfBestSTDPlus, stCluster.fAvgZeroMixLenHalfBestSTDMinus);
    (*m_pOfs) << " === " << endl;

    (*m_pOfs) << "vZeroMixRatioNorm: " << endl;
    PrintFVector(vZeroMixRatioNorm, stCluster.fNormAvgZeroMixLenRatio, stCluster.fNormAvgZeroMixLenHalfBestRatio,
                 stCluster.fNormAvgZeroMixSTDPlus, stCluster.fNormAvgZeroMixSTDMinus,
                 stCluster.fNormAvgZeroMixLenHalfBestSTDPlus, stCluster.fNormAvgZeroMixLenHalfBestSTDMinus);
    (*m_pOfs) << " === " << endl;

    (*m_pOfs) << "Negative: " << endl << strTypeNegSum;
    (*m_pOfs) << " === " << endl;

    (*m_pOfs) << "Positive: " << endl << strTypePositiveSum;
    (*m_pOfs) << " === " << endl;
}

//Print vector did two things
//(1) Print Vector
//(2) Calculat the average value for "fAvg" and "fAvgHalfBest"
void ClsLearning::PrintFVector(vector<float>& vValue, float& fAvg, float& fAvgHalfBest,
                               float& fAvgSTDPlus, float& fAvgSTDMinus,
                               float& fAvgHalfBestSTDPlus, float& fAvgHalfBestSTDMinus)
{
    float fSum = 0;
    for(vector<float>::iterator itr = vValue.begin(); itr != vValue.end(); itr++)
    {
        (*m_pOfs) << FloatToStr(*itr) << ", ";
        fSum += *itr;
    }
    (*m_pOfs) << endl;
    //Also Get Average value -->
    //(1) Average whole
    fAvg = (float)fSum / vValue.size();
    float fSTD = CalculateSD(vValue);
    fAvgSTDPlus = fAvg + fSTD;
    fAvgSTDMinus = fAvg - fSTD;

    (*m_pOfs) << "Avg          : " << FloatToStr(fAvg) << ", STD: " << FloatToStr(fSTD)
         << "  ---- +1STD, +2STD, +3STD: "
         << FloatToStr(fAvg + fSTD) << ", "
         << FloatToStr(fAvg + 2*fSTD) << ", "
         << FloatToStr(fAvg + 3*fSTD)
         << "  ---- -1STD, -2STD, -3STD: "
         << FloatToStr(fAvg - fSTD) << ", "
         << FloatToStr(fAvg - 2*fSTD) << ", "
         << FloatToStr(fAvg - 3*fSTD) << endl;
    //(2) Average sorted half --> Go!!
    sort(vValue.begin(), vValue.end());
    int iNum = 0;
    fSum = 0;
    vector<float> vHalfBest;
    for(vector<float>::iterator itr = vValue.end() - 1; itr >= vValue.begin(); itr--)
    {
        if(itr < vValue.begin() + (vValue.size() / 2))
            break;

        iNum++;
        fSum += *itr;
        vHalfBest.push_back(*itr);
    }

    fAvgHalfBest = (float)fSum / iNum;
    fSTD = CalculateSD(vHalfBest);
    fAvgHalfBestSTDPlus = fAvgHalfBest + fSTD;
    fAvgHalfBestSTDMinus = fAvgHalfBest - fSTD;

    (*m_pOfs) << "Avg Half Best: " << FloatToStr(fAvgHalfBest) << ", STD: " << FloatToStr(fSTD)
         << "  ---- +1STD, +2STD, +3STD: "
         << FloatToStr(fAvgHalfBest + fSTD) << ", "
         << FloatToStr(fAvgHalfBest + 2*fSTD) << ", "
         << FloatToStr(fAvgHalfBest + 3*fSTD)
         << "  ---- -1STD, -2STD, -3STD: "
         << FloatToStr(fAvgHalfBest - fSTD) << ", "
         << FloatToStr(fAvgHalfBest - 2*fSTD) << ", "
         << FloatToStr(fAvgHalfBest - 3*fSTD) << endl;
}

float ClsLearning::CalculateSD(vector<float>& vValue)
{
    if(vValue.empty())
        return 0;

    float fSum = 0;
    for(vector<float>::iterator itr = vValue.begin(); itr != vValue.end(); itr++)
    {
        fSum += *itr;
    }
    float fAvg = (float)fSum / vValue.size();

    float fSD = 0;
    for(vector<float>::iterator itr = vValue.begin(); itr != vValue.end(); itr++)
    {
        fSD += pow(*itr - fAvg, 2);
    }

    return sqrt(fSD / vValue.size());
}


