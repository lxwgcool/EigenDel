#include "clsvcf1000genome.h"
#include <fstream>
#include <iostream>
#include <unistd.h>

ClsVcf1000Genome::ClsVcf1000Genome()
{

}

int ClsVcf1000Genome::ParseVcf(string strVcfFile, string strSampleName)
{
    if(::access(strVcfFile.c_str(), 0) != 0)
    {
        cout << "Error: VCF File does not existed!" << endl;
        return 1;
    }

    vector<St_SV> vSV;
    vSV.clear();
    St_SV stSV;
    ifstream ifsVcf;
    ifsVcf.open(strVcfFile.c_str());
    string strLine;
    int iSampleIndex = -1;
    vector<string> vSeq;
    while(!ifsVcf.eof())
    {
        getline(ifsVcf, strLine);

        //Skip the comment line
        if(strLine.substr(0, 1) == "#")
        {
            if(strSampleName != "" && strLine.substr(0, 6) == "#CHROM")
            {
                //Try to find the column index of specific person -->
                iSampleIndex = GetColNumber(strLine, strSampleName);                
                //<--
            }            
            continue;
        }

        stSV.Clear();        

        //The regular line
        string::size_type iStart = 0;
        string::size_type iOffset = 0;
        int iLen = 0;
        string strTmp = "";
        // 1. 1st column --> Chome Index
        iOffset = strLine.find('\t', iStart);
        iLen = iOffset - iStart;
        strTmp = strLine.substr(iStart, iLen);
        stSV.strChrom = strTmp; //atoi(strTmp.c_str());
        iStart = iOffset + 1;

        // 2: 2nd column --> Position
        iOffset = strLine.find('\t', iStart);
        iLen = iOffset - iStart;
        strTmp = strLine.substr(iStart, iLen);
        stSV.iPos = atoi(strTmp.c_str());
        iStart = iOffset + 1;

        // 3: 4th column --> Ref Seq
        iOffset = strLine.find('\t', iStart);
        iStart = iOffset + 1;
        iOffset = strLine.find('\t', iStart);
        iLen = iOffset - iStart;
        strTmp = strLine.substr(iStart, iLen);
        stSV.strRef = strTmp;
        iStart = iOffset + 1;

        // 4: 8th column --> Info
        for(int i=0; i<3; i++)
        {
            iOffset = strLine.find('\t', iStart);
            iStart = iOffset + 1;
        }
        iOffset = strLine.find('\t', iStart);
        iLen = iOffset - iStart;
        strTmp = strLine.substr(iStart, iLen);
          //(1) iLen -- SVLEN
        iStart = strTmp.find("SVLEN");
        if(iStart != string::npos) // find it
        {
            iStart += 5 + 1; // length (SVLEN) is 5, and "=" is 1
            iOffset = strTmp.find(";", iStart);
            if(iOffset == string::npos) // going to the end of current sequence
            {
                iOffset = strTmp.length();
            }
            iLen = iOffset - iStart;
            stSV.iLen = atoi(strTmp.substr(iStart, iLen).c_str());
        }
          //(2) SVType -- SVTYPE
        iStart = strTmp.find("SVTYPE");
        if(iStart != string::npos) // find it
        {
            iStart += 6 + 1; // length (SVTYPE) is 6, and "=" is 1
            iOffset = strTmp.find(";", iStart);
            if(iOffset == string::npos) // going to the end of current sequence
            {
                iOffset = strTmp.length();
            }
            iLen = iOffset - iStart;
            stSV.strType = strTmp.substr(iStart, iLen);
        }
          //(3) iEnd -- ;END=
        iStart = strTmp.find(";END=");        
        if(iStart != string::npos) // find it
        {
            iStart += 5; // length (;END=) is 5
            iOffset = strTmp.find(";", iStart);
            if(iOffset == string::npos) // going to the end of current sequence
            {
                iOffset = strTmp.length();
            }
            iLen = iOffset - iStart;
            stSV.iEnd = atoi(strTmp.substr(iStart, iLen).c_str());
        }

        //For the specific sample -->
        if(iSampleIndex != -1)
        {
            vSeq.clear();
            SpliteSeq(strLine, '\t', vSeq);

            if(vSeq.size() > iSampleIndex)
            {
                stSV.stSample.strName = strSampleName;
                //Check Genotype -->
                string strSeq = vSeq[iSampleIndex];
                //cout << "Genotype: " << strSeq << endl;
                if(strSeq == "." || strSeq == "0") // Why genotype sometimes is "." ??
                {
                    stSV.stSample.iGT1 = 0;
                    stSV.stSample.iGT2 = 0;
                }
                else if(strSeq == "1")
                {
                    stSV.stSample.iGT1 = 1;
                    stSV.stSample.iGT2 = 0;
                }
                else if(strSeq == "2")
                {
                    stSV.stSample.iGT1 = 1;
                    stSV.stSample.iGT2 = 1;
                }
                else
                {
                    stSV.stSample.iGT1 = atoi(strSeq.substr(0, 1).c_str());
                    stSV.stSample.iGT2 = atoi(strSeq.substr(2, 1).c_str());
                }
                //<--
            }
        }
        //<--

        //put SV into vector
        vSV.push_back(stSV);
    }
    ifsVcf.close();

    //Separate vSv to vGenoSv -->
    m_vChromSv.clear();
    St_ChromSV stChromSv;
    for(vector<St_SV>::iterator itr = vSV.begin(); itr != vSV.end(); itr++)
    {
        bool bFindChrom = false;
        for(vector<St_ChromSV>::iterator itrChromSv = m_vChromSv.begin(); itrChromSv != m_vChromSv.end(); itrChromSv++)
        {
            if(itrChromSv->strChrom == itr->strChrom)
            {
                itrChromSv->vSv.push_back(*itr);
                bFindChrom = true;
                break;
            }
        }

        if(!bFindChrom)
        {
            stChromSv.Clear();
            stChromSv.strChrom = itr->strChrom;
            stChromSv.vSv.push_back(*itr);
            m_vChromSv.push_back(stChromSv);
        }
    }
    vSV.clear(); //release heap
    //<--
    return 0;
}

void ClsVcf1000Genome::GetDeletion(vector<St_ChromSV>& vChromSvDEL, string strSampleName, string strChrom)
{
    if(m_vChromSv.empty())
        return;

    vChromSvDEL.clear();
    St_ChromSV stChromSvDEL;
    for(vector<St_ChromSV>::iterator itrChromSv = m_vChromSv.begin(); itrChromSv != m_vChromSv.end(); itrChromSv++)
    {
        if(strChrom == "") // This means we need to consider the whole geonome
        {}
        else
        {
            if(itrChromSv->strChrom != strChrom)
                continue;
        }

        stChromSvDEL.Clear();
        stChromSvDEL.strChrom = itrChromSv->strChrom;
        for(vector<St_SV>::iterator itr = itrChromSv->vSv.begin(); itr != itrChromSv->vSv.end(); itr++)
        {
            if(itr->strType.find("DEL") != string::npos) //find it
            {
                if(strSampleName != "")
                {
                    //Get the deletion only contained bystr SampleName
                    if(itr->stSample.strName == strSampleName &&
                       (itr->stSample.iGT1 == 1 || itr->stSample.iGT2 == 1))
                        stChromSvDEL.vSv.push_back(*itr);
                }
                else
                    stChromSvDEL.vSv.push_back(*itr);
            }
        }
        vChromSvDEL.push_back(stChromSvDEL);

        //Only need to record one chromosome
        if(strChrom != "")
            break;
    }
}

int ClsVcf1000Genome::GetColNumber(string& strSeq, string strSampleName) //Get Column Number
{
    vector<string> vSeq;
    SpliteSeq(strSeq, '\t', vSeq);
    int iIndex = 0;
    int iSampleIndex = -1;
    for(vector<string>::iterator itr = vSeq.begin(); itr != vSeq.end(); itr++, iIndex++)
    {
        if(*itr == strSampleName)
        {
            iSampleIndex = iIndex;
            break;
        }
    }
    return iSampleIndex;
}

void ClsVcf1000Genome::SpliteSeq(string& strSeq, char cDelimiter, vector<string>& vSeq)
{
    vSeq.clear();
    string::size_type iOffset = 0;
    string::size_type iPos = 0;
    int iLen = 0;
    string strSubStr = "";

    while(iPos != string::npos)
    {
        iPos = strSeq.find(cDelimiter, iOffset);
        if(iPos != string::npos)
        {
            iLen = iPos - iOffset;
            strSubStr = strSeq.substr(iOffset, iLen);
            vSeq.push_back(strSubStr);
            iOffset = iPos + 1;
        }
    }
    return;
}
