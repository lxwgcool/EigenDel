#ifndef CLSBLAST_H
#define CLSBLAST_H
#include <string>
#include <vector>
using namespace std;

struct St_BlastAlignUnit
{
    string strName;
    int iAlignNum;
    float fAlignPercent;
    int iStartPos; // This is for subject sequence --> subject is, for example, reference
    int iEndPos;   // This is for subject sequence

    int iStartPosQuery; // This is for query sequence --> query is, for example, reads
    int iEndPosQuery;   // This is for query sequence

    float fIdentifyRatio;
    float fGapRatio;

    St_BlastAlignUnit()
    {
        Clear();
    }

    void Clear()
    {
        strName = "";
        fAlignPercent = -1;
        iStartPos = -1;
        iEndPos = -1;
        iAlignNum = -1;
        fIdentifyRatio = 0;
        fGapRatio = 0;
        iStartPosQuery = -1;
        iEndPosQuery = -1;
    }
    bool IsRC() // If it is reverse complementary
    {
        if(iStartPos < 0 || iEndPos < 0)
            return false;
        if(iStartPos > iEndPos)
            return true;
        else
            return false;
    }
};

struct St_BlastResult
{
    vector<St_BlastAlignUnit> vBlastAlign;

    void Init()
    {
        vBlastAlign.clear();
    }
};

enum En_FastaType{ftQuery=0, ftRef, ftMax};

class ClsBlast
{
public:
    ClsBlast();

public:
    //Create Fa File
    string CreatFaFile(string strName, string& strSeq, En_FastaType enFastaType);

    //Make comparison by blastn
    string TwoFastaFileAlign(string strQuerySeqPath, string strTargetSeqPath,
                             int iSeqSizeType=0); //0: normal size type, 1: small size type
                                                //2: large size type

    St_BlastResult ParseResult(string strRstPath, vector<St_BlastResult>& vBlastResult,
                               bool bRecord);

    int TwoSeqBlast(string& strQueryPath, string& strSubjectPath, bool bShowFlank=false);

    void GetTopNResultsFromTwoSeqBlast(string& strQueryPath, string& strSubjectPath, vector<St_BlastAlignUnit>& vBlastAlign, int iTopN);

private:
    int GetBlastResult(string& strQueryPath, string& strSubjectPath, St_BlastResult& stBlastResult);
    vector<St_BlastResult> m_vBlastResultSum;
};

#endif // CLSBLAST_H
