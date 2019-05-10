#ifndef CLSCOMPARISON_H
#define CLSCOMPARISON_H
#include <string>
#include <vector>
#include "../../../ShareLibrary/clsvcf1000genome.h"
using namespace std;

struct St_BreakPoint
{
    string strChromID;
    int iStart;
    int iEnd;
};

enum En_Software{swPindel=0, swSprites, swSvGASVpro, swCNVnator, swDelly, swLumpy, swEigenDel, swMax};

struct St_SvResult
{
    string strChrom;
    int iStdSv;
    int iSvFind;
    int iSvPositive; //From SvFind
    int iSvNegative; //From SvFind

    St_SvResult()
    {
        Clear();
    }

    St_SvResult(string strV1, int iV2, int iV3, int iV4, int iV5)
    {
        this->strChrom = strV1;
        this->iStdSv = iV2;
        this->iSvFind = iV3;
        this->iSvPositive = iV4;
        this->iSvNegative = iV5;
    }

    void Clear()
    {
        strChrom = "";
        iStdSv = 0;
        iSvFind = 0;
        iSvPositive = 0;
        iSvNegative = 0;
    }

    float CalcAccuracy()
    {
        if(iSvFind == 0)
            return -1;

        float fAccuracy = (float)iSvPositive / iSvFind;
        return fAccuracy;
    }

    float CalcSensitivity()
    {
        if(iStdSv == 0)
            return -1;

        float fSensitivity = (float)iSvPositive / iStdSv;
        return fSensitivity;
    }

    float CalcF1Score()
    {
        float fAccuracy = CalcAccuracy();
        float fSensitivity = CalcSensitivity();
        if(fAccuracy + fSensitivity == 0)
            return -1;

        float fF1 = (float)2 * fAccuracy * fSensitivity / (fAccuracy + fSensitivity);
        return fF1;
    }
};

struct St_ComparisonResult
{
    En_Software enSoftware;
    vector<St_SvResult> vSvResult; //Chrom by Chrom

    int GetTotalStdSv()
    {
        int iStdSvSum = 0;
        for(vector<St_SvResult>::iterator itr = vSvResult.begin(); itr != vSvResult.end(); itr++)
        {
            if(itr->iStdSv == 0)
                continue;

            iStdSvSum += itr->iStdSv;
        }
        return iStdSvSum;
    }

    int GetSvFind()
    {
        int iSvFind = 0;
        for(vector<St_SvResult>::iterator itr = vSvResult.begin(); itr != vSvResult.end(); itr++)
        {
            if(itr->iStdSv == 0)
                continue;
            iSvFind += itr->iSvFind;
        }
        return iSvFind;
    }

    int GetSvPositive()
    {
        int iSvPositive = 0;
        for(vector<St_SvResult>::iterator itr = vSvResult.begin(); itr != vSvResult.end(); itr++)
        {
            if(itr->iStdSv == 0)
                continue;

            iSvPositive += itr->iSvPositive;
        }
        return iSvPositive;
    }

    int GetSvNegative()
    {
        int iSvNegative = 0;
        for(vector<St_SvResult>::iterator itr = vSvResult.begin(); itr != vSvResult.end(); itr++)
        {
            if(itr->iStdSv == 0)
                continue;

            iSvNegative += itr->iSvNegative;
        }
        return iSvNegative;
    }

    float CalcAccuracy()
    {
        int iSvFind = GetSvFind();
        if(iSvFind == 0)
            return -1;

        int iSvPositive = GetSvPositive();
        float fAccuracy = (float)iSvPositive / iSvFind;
        return fAccuracy;
    }

    float CalcSensitivity()
    {
        int iSvPositive = GetSvPositive();
        int iStdSv = GetTotalStdSv();
        if(iStdSv == 0)
            return -1;

        float fSensitivity = (float)iSvPositive / iStdSv;
        return fSensitivity;
    }

    float CalcF1Score()
    {
        float fAccuracy = CalcAccuracy();
        float fSensitivity = CalcSensitivity();
        if(fAccuracy + fSensitivity == 0)
            return -1;

        float fF1 = (float)2 * fAccuracy * fSensitivity / (fAccuracy + fSensitivity);
        return fF1;
    }
};

class ClsComparison
{
public:
    ClsComparison();

public:    
    void ParsePindel(string strFilePath, vector<St_BreakPoint>& vBP);
    void ParseSprites(string strFilePath, vector<St_BreakPoint>& vBP);
    void ParseGASVpro(string strFilePath, vector<St_BreakPoint>& vBP);
    void ParseCNVnator(string strFilePath, vector<St_BreakPoint>& vBP);

    void CompareStdSVWithPindel(string strFilePath, vector<St_ChromSV>& vChromSvDEL);
    void CompareStdSvWithSprites(string strFilePath, vector<St_ChromSV>& vChromSvDEL);
    void CompareStdSVWithLumpy(string strFilePath, vector<St_ChromSV>& vChromSvDEL);
    void CompareStdSVWithDelly(string strFilePath, vector<St_ChromSV>& vChromSvDEL);   

    void CompareStdSVWithGASVpro(string strFilePath, vector<St_ChromSV>& vChromSvDEL);
    void CompareStdSVWithCNVnator(string strFilePath, vector<St_ChromSV>& vChromSvDEL);

public:
    void GetResultFromVcf(string strFilePath, vector<St_ChromSV>& vChromSvDEL,
                          St_ComparisonResult& stCR); //CR: comparison result
    void GetByCompareStd(vector<St_ChromSV>& vStdSvDEL,
                         vector<St_ChromSV> vSoftwareSvDEL,
                         St_ComparisonResult& stCR);
    void PrintResult(St_ComparisonResult& stCR);    
};

#endif // CLSCOMPARISON_H
