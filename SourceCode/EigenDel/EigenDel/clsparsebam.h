#ifndef CLSPARSEBAM_H
#define CLSPARSEBAM_H

#include <string>
#include <vector>
#include "../../../ShareLibrary/clsbasealgorithm.h"
#include "../../../ShareLibrary/clsvcf1000genome.h"
#include "../../../ShareLibrary/clsfastareader.h"
#include "../../../ShareLibrary/clsblast.h"
//#include "clsdrawimage.h"
#include <map>

#include "../../../ShareLibrary/bamtools/include/api/BamReader.h"
using namespace BamTools;

using namespace std;

enum En_ClipType{ctSoft, ctHard, ctMax};
enum En_MapType{mtMat1End, mtMat2End, mtMat2Start, mtMat1Start, mtMax};

struct St_RegReads
{
    int iReadsMapPos;
    string strName;

    St_RegReads()
    {
        Clear();
    }

    void Clear()
    {
        iReadsMapPos = -1;
        strName = "";
    }
};

struct St_ClipReads
{
    int iReadsMapPos;
    string strName;
    string strQuerySeq;
    string strAlignSeq;
    vector<int> vPos;
    vector<string> vClipSeq;

    int iQueryMapPos;
    int iQueryMapLength;
    int iQueryClipPos;
    int iQueryClipLength;


    int iMateMapPos;
    string strMateName;
    string strMateQuerySeq;
    string strMateAlignSeq;
    vector<int> vMatePos;
    vector<string> vMateClipSeq;

    int iInsertSize;
    bool bFirstMate;
    bool bSecondMate;
    bool bMateMapped;

    St_ClipReads()
    {
        Clear();
    }

    void Clear()
    {
        iReadsMapPos = -1;
        strName = "";
        strQuerySeq = "";
        strAlignSeq = "";
        vPos.clear();
        vClipSeq.clear();

        iMateMapPos = -1;
        strMateName = "";
        strMateQuerySeq = "";
        strMateAlignSeq = "";
        vMatePos.clear();
        vMateClipSeq.clear();

        iInsertSize = -1;
        bFirstMate = false;
        bSecondMate = false;
        bMateMapped = false;

        iQueryMapPos = -1;
        iQueryMapLength = -1;
        iQueryClipPos = -1;
        iQueryClipLength = -1;
    }
};

//Let's get map length
struct St_BorderSCR //Border Soft-Cliped-Reads
{
    int iReadsMapPos;
    int iMatPos;
    int iClipPos;
    En_ClipPart enClipPart;
    string strClipSeq;
    string strQuerySeq;
    int iInsertSize;
    string strName;

    //The position info of reads
    int iQueryMapPos;
    int iQueryMapLength;
    int iQueryClipPos;
    int iQueryClipLength;

    //For debug
    bool bSvContribute;
    int iDiscCountNum;

    En_MapType enMapType;
    bool bFirstMate;
    bool bSecondMate;
    bool bMateMapped;

    bool bDrippedIntoDRGroup;

    St_BorderSCR()
    {
        Clear();
    }

    void Clear()
    {
        iReadsMapPos = -1;
        iMatPos = -1;
        iClipPos = -1;
        enClipPart = cpMax;
        strClipSeq = " ";
        strQuerySeq = " ";
        iInsertSize = -1;
        strName = " ";

        //-->For debug
        bSvContribute = false;
        iDiscCountNum = 0;
        //<--

        enMapType = mtMax;
        bFirstMate = false;
        bSecondMate = false;
        bMateMapped = false;

        bDrippedIntoDRGroup = false;

        iQueryMapPos = -1;
        iQueryMapLength = -1;
        iQueryClipPos = -1;
        iQueryClipLength = -1;
    }

    bool IsLeftPair()
    {
        if(iReadsMapPos < iMatPos)
            return true;
        else
            return false;
    }

    bool IsRightPair()
    {
        if(iReadsMapPos > iMatPos)
            return true;
        else
            return false;
    }

    En_ClipPart GetClipPart()
    {
        if(this->iReadsMapPos == this->iClipPos)
            return cpLeft;
        else
            return cpRight;
    }
};

struct St_ChromBorderSCR
{
    vector<St_BorderSCR> vSCR;
    string strChrom;

    St_ChromBorderSCR():strChrom("")
    {}

    void Clear()
    {
        vSCR.clear();
        strChrom = "";
    }
};

//-->For discondent reads
struct St_DiscordantReads //Actually, it contains two reads (those two reads are discordant)
{
    //For left reads
    int iReadsPos;  //left side reads
    string strReadsAlignSeq;
    string strReadsClipSeq;
    bool bClip;

    //For right reads
    int iMatePos; //right side reads
    string strMateAlignSeq;
    string strMateClipSeq;
    bool bMateClip;

    //Discordant info
    int iInsertSize;

    //The position info of reads
    int iQueryMapPos;
    int iQueryMapLength;
    int iQueryClipPos;
    int iQueryClipLength;
    string strClipSeq;

    int iMatQueryMapPos;
    int iMatQueryMapLength;
    int iMatQueryClipPos;
    int iMatQueryClipLength;
    string strMatClipSeq;

    St_DiscordantReads()
    {
        Clear();
    }

    void Clear()
    {
        //Left Reads
        iReadsPos = 0;
        strReadsAlignSeq = "";
        strReadsClipSeq = "";
        bClip = false;

        //Right reads
        iMatePos = 0;
        strMateAlignSeq = "";
        strMateClipSeq = "";
        bMateClip = false;

        //Discordant Info
        iInsertSize = 0;

        //This is for draw the pics
        iQueryMapPos = -1;
        iQueryMapLength = -1;
        iQueryClipPos = -1;
        iQueryClipLength = -1;
        strClipSeq = "";

        iMatQueryMapPos = -1;
        iMatQueryMapLength = -1;
        iMatQueryClipPos = -1;
        iMatQueryClipLength = -1;
        strMatClipSeq = "";
    }
};

struct St_ChromDiscordantReads
{
    vector<St_DiscordantReads> vDR;
    string strChrom;

    St_ChromDiscordantReads():strChrom("")
    {}

    void Clear()
    {
        vDR.clear();
        strChrom = "";
    }
};

//-->

class ClsParseBam
{
public:
    ClsParseBam();
    ~ClsParseBam();

public:
    void ReadBamFile(string strBamFilePath, vector<BamAlignment>& vAl);
    void CalcInsertSize(string strBamFilePath, string strPicard);
    void CalcAvgDepth(string strBamFilePath);

    //Get to types of reads
    void GetBorderSCR(string strBamFilePath, vector<St_Fasta>& vFasta,
                      vector<St_ChromBorderSCR>& vChromBorderSCR);//vector<St_BorderSCR>& vBorderSCR);
    void GetDiscordantReads(string strBamFilePath, vector<St_Fasta>& vFasta,
                            vector<St_ChromDiscordantReads>& vChromDR);//vector<St_DiscordantReads>& vDiscdRreads);

    //We do both SCR and DR at the same time  -->
    void GetSCRAndDR(string strBamFilePath, vector<St_Fasta>& vFasta,
                     vector<St_ChromBorderSCR>& vChromBorderSCR,
                     vector<St_ChromDiscordantReads>& vChromDR);

    void GetSCRByAl(BamAlignment& al,
                    vector<St_ClipReads>& vClipReads, St_ClipReads& stClipReads,
                    St_ChromBorderSCR& stChromSCR, St_BorderSCR& stBorderSCR,
                    vector<St_Fasta>& vFasta, vector<St_ChromBorderSCR>& vChromBorderSCR);
    void PrintTestingInfoSCR(vector<St_ClipReads>& vClipReads, vector<St_ChromBorderSCR>& vChromBorderSCR);

    void GetDRByAl(BamAlignment& al,
                   St_DiscordantReads& stDR, St_ChromDiscordantReads& stChromDR,
                   vector<St_Fasta>& vFasta, vector<St_ChromDiscordantReads>& vChromDR,
                   float& fMinThreshold, float& fMaxThreshold);
    void PrintTestingInfoDR(vector<St_ChromDiscordantReads>& vChromDR);
    //<--

public:
    void SetMeanInsertSize(float fV1);
    void SetStdDevInsertSize(float fV1);
    void SetAvgDepth(float fV1);

    float GetMeanInsertSize();
    float GetStdDevInsertSize();
    float GetAvgDepth();

    void SetOfstram(ofstream& ofs);
    void SetChromName(string strChromName);

private:
    float m_fMeanInsertSize;
    float m_fStdDevInsertSize;
    float m_fAvgDepth;
    ofstream* m_pOfs;
    string m_strChromName;
};

#endif // CLSPARSEBAM_H
