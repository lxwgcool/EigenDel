#ifndef CLSVCF1000GENOME_H
#define CLSVCF1000GENOME_H

#include <string>
#include <vector>
using namespace std;

//This is for the VCF file of Structural Variation in 1000 Genome Project.
struct St_Specimen // For specific person
{
    string strName;
    int iGT1; //genotype 1
    int iGT2; //genotype 2

    St_Specimen():strName(""),iGT1(-1),iGT2(-1)
    {}

    void Clear()
    {
        strName = "";
        iGT1 = -1;
        iGT2 = -1;
    }
};

enum En_SVType{svtDel, svtIns, svtMax};
struct St_SV
{
    string strChrom;
    int iPos; // beginning
    int iEnd;
    string strType;
    int iLen;
    string strRef;
    St_Specimen stSample;

    St_SV():strChrom(""), iPos(-1), iEnd(-1), strType(""), iLen(-1), strRef("")
    {}

    St_SV(string strV1, int iV2, int iV3)
    {
        this->strChrom = strV1;
        this->iPos = iV2;
        this->iEnd = iV3;
    }

    void Clear()
    {
        strChrom = "";
        iPos = -1;
        iEnd = -1;
        strType = "";
        iLen = -1;
        strRef = "";
        stSample.Clear();       
    }

    bool operator == (const St_SV& rhs) const
    {
        if(this->strChrom == rhs.strChrom &&
           this->iPos == rhs.iPos &&
           this->iEnd == rhs.iEnd &&
           this->strType == rhs.strType)
            return true;
        else
            return false;
    }

    int GetLen()
    {
        return this->iEnd - this->iPos;
    }
};

struct St_ChromSV
{
    vector<St_SV> vSv;
    string strChrom;

    St_ChromSV():strChrom("")
    {}



    void Clear()
    {
        vSv.clear();
        strChrom = "";
    }
};

class ClsVcf1000Genome
{
public:
    ClsVcf1000Genome();

public:
    int ParseVcf(string strVcfFile, string strSampleName="");
    void GetDeletion(vector<St_ChromSV>& vChromSvDEL, string strSampleName="", string strChrom="");

private:
    int GetColNumber(string& strSeq, string strSampleName); //Get Column Number
    void SpliteSeq(string& strSeq, char cDelimiter, vector<string>& vSeq);

private:
    vector<St_ChromSV> m_vChromSv;
};

#endif // CLSVCF1000GENOME_H
