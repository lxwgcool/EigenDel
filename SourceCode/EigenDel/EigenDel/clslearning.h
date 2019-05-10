#ifndef CLSLEARNING_H
#define CLSLEARNING_H
#include "clssvdeldetect.h"

struct St_GroupCluster
{
    vector<St_DRGroup> vGroup;

    //Check the key value
    float fAvgZeroPureLenRatio;
    float fAvgZeroPureSTDPlus;
    float fAvgZeroPureSTDMinus;
    float fAvgZeroPureLenHalfBestRatio;
    float fAvgZeroPureLenHalfBestSTDPlus;
    float fAvgZeroPureLenHalfBestSTDMinus;


    float fAvgZeroWeekLenRatio;
    float fAvgZeroWeekSTDPlus;
    float fAvgZeroWeekSTDMinus;
    float fAvgZeroWeekLenHalfBestRatio;
    float fAvgZeroWeekLenHalfBestSTDPlus;
    float fAvgZeroWeekLenHalfBestSTDMinus;


    float fNromAvgZeroWeekLenRatio;
    float fNormAvgZeroWeekSTDPlus;
    float fNormAvgZeroWeekSTDMinus;
    float fNormAvgZeroWeekLenHalfBestRatio;
    float fNormAvgZeroWeekLenHalfBestSTDPlus;
    float fNormAvgZeroWeekLenHalfBestSTDMinus;


    float fAvgZeroMixLenRatio;
    float fAvgZeroMixSTDPlus;
    float fAvgZeroMixSTDMinus;
    float fAvgZeroMixLenHalfBestRatio;
    float fAvgZeroMixLenHalfBestSTDPlus;
    float fAvgZeroMixLenHalfBestSTDMinus;


    float fNormAvgZeroMixLenRatio;
    float fNormAvgZeroMixSTDPlus;
    float fNormAvgZeroMixSTDMinus;
    float fNormAvgZeroMixLenHalfBestRatio;
    float fNormAvgZeroMixLenHalfBestSTDPlus;
    float fNormAvgZeroMixLenHalfBestSTDMinus;


    int iPositive;
    int iNegative;

    //For statistic filter threshold

    St_GroupCluster()
    {
        Clear();
    }

    void Clear()
    {
        fAvgZeroPureLenRatio = 0;
        fAvgZeroPureLenHalfBestRatio = 0;
        fAvgZeroPureSTDPlus = 0;
        fAvgZeroPureSTDMinus = 0;
        fAvgZeroPureLenHalfBestSTDPlus = 0;
        fAvgZeroPureLenHalfBestSTDMinus = 0;

        fAvgZeroWeekLenRatio = 0;
        fAvgZeroWeekLenHalfBestRatio = 0;
        fAvgZeroWeekSTDPlus = 0;
        fAvgZeroWeekSTDMinus = 0;
        fAvgZeroWeekLenHalfBestSTDPlus = 0;
        fAvgZeroWeekLenHalfBestSTDMinus = 0;

        fNromAvgZeroWeekLenRatio = 0;
        fNormAvgZeroWeekLenHalfBestRatio = 0;
        fNormAvgZeroWeekSTDPlus = 0;
        fNormAvgZeroWeekSTDMinus = 0;
        fNormAvgZeroWeekLenHalfBestSTDPlus = 0;
        fNormAvgZeroWeekLenHalfBestSTDMinus = 0;

        fAvgZeroMixLenRatio = 0;
        fAvgZeroMixLenHalfBestRatio = 0;
        fAvgZeroMixSTDPlus = 0;
        fAvgZeroMixSTDMinus = 0;
        fAvgZeroMixLenHalfBestSTDPlus = 0;
        fAvgZeroMixLenHalfBestSTDMinus = 0;

        fNormAvgZeroMixLenRatio = 0;
        fNormAvgZeroMixLenHalfBestRatio = 0;
        fNormAvgZeroMixSTDPlus = 0;
        fNormAvgZeroMixSTDMinus = 0;
        fNormAvgZeroMixLenHalfBestSTDPlus = 0;
        fNormAvgZeroMixLenHalfBestSTDMinus = 0;

        iPositive = 0;
        iNegative = 0;
    }
};

class ClsLearning
{
public:
    ClsLearning();

public:    
    void PrintFeatures(vector<St_DRGroup>& vDRGroup, vector<St_SCRGroup>& vSCRGroup,
                        vector<St_SV>& vSvDEL);

    void PrintClusterElements(vector<St_DRGroup>& vDRGroup, vector<St_SV>& vSvDEL);

    void ClusterByPython(vector<St_DRGroup>& vDRGroup, vector<St_SV>& vSvDEL,
                         string strPythonCode, St_FinalResult& stFinalResult);

private:
    void GetPositiveAndNegativeGroup(vector<St_DRGroup>& vDRGroup, vector<St_SV>& vSvDEL, //vector<St_SCRGroup>& vSCRGroup,
                                     vector<St_DRGroup>& vPositiveGroup, vector<St_DRGroup>& vNegativeGroup);

    void PrintCandiDepthDetail(ofstream& ofs, vector<St_DRGroup>::iterator itr, int iIndex);

    int GetNormalizedValue(int iRealLen, int iSegLen, int iStdNormLen = 1000);
    string GetElementValue(St_GroupFeatures& stGF);
    string GetNormDepth(St_GroupFeatures& stGF);

    void CalcFeatureOfEachCluster(St_GroupCluster& stCluster);
    void DoFilterOfEachCluster(St_GroupCluster& stCluster, int iClusterIndex);

    void PrintFVector(vector<float>& vValue, float& fAvg, float& fAvgHalfBest,
                      float& fAvgSTDPlus, float& fAvgSTDMinus,
                      float& fAvgHalfBestSTDPlus, float& fAvgHalfBestSTDMinus);

    void FilterGoodCluster(St_GroupCluster& stCluster, int iClusterIndex);
    void FilterBadCluster(St_GroupCluster& stCluster, int iClusterIndex);

    float CalculateSD(vector<float>& vValue);

public:
    float GetAvgDepth();
    void SetAvgDepth(float fV1);
    void SetChromName(string strChromName);
    void SetSampleName(string strSampleName);
    void SetOfstream(ofstream& ofs);

private:
    float m_fAvgDepth;
    string m_strChromName;
    string m_strSampleName;
    ofstream* m_pOfs;
};

#endif // CLSLEARNING_H
