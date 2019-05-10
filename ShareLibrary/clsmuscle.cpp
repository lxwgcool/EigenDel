#include "clsmuscle.h"
#include "clsbasealgorithm.h"
#include <unistd.h>
#include <map>

ClsMuscle::ClsMuscle()
{
    m_strMusclePath = "../../../../ShareLibrary/Muscle/muscle3.8.31_i86linux64";
    if(access(m_strMusclePath.c_str(), 0) != 0)
        cout << "Init Muscle Failed!" << endl;
}

ClsMuscle::~ClsMuscle()
{
}

string ClsMuscle::Run(string strFaPath)
{
    if(access(strFaPath.c_str(), 0) != 0)
    {
        cout << "Aligned Source File doesn't exist!" << endl;
        return "";
    }
    //Run muscle
    string strOutput = "./result.txt";
    string strCommand = m_strMusclePath +
                        " -in " + strFaPath +
                        " -clw -quiet " +
                        " -out " + strOutput;
    system(strCommand.c_str());
    return strOutput;
}

string ClsMuscle::MergeAlignedSeq(string strSeqPath)
{
    if(access(strSeqPath.c_str(), 0) != 0)
    {
        cout << "MergeAlignedSeq File doesn't exist!" << endl;
        return "";
    }

    //Try to merge the aligned file -->
    ifstream ifsAlignedFile;
    ifsAlignedFile.open(strSeqPath.c_str(), ios::in);
    string strLine = "";
    getline(ifsAlignedFile, strLine); // --> Skip the first line
    vector<string> vSeq;
    string strTarget = "";
    string::size_type iIndex = 0;
    while(!ifsAlignedFile.eof())
    {
        getline(ifsAlignedFile, strLine);
        //Parse the line --> Go Notice: the blan is "space" -- >Go!!
        if(strLine == "")
        {
            iIndex = 0;
            continue;
        }
        //Now go to check the target sequence
        if(iIndex < vSeq.size())
        {
            strTarget = strLine.substr(16, strLine.size() - 16);
            vSeq[iIndex] += strTarget;
            iIndex++;
        }
        else
        {
            strTarget = strLine.substr(16, strLine.size() - 16);
            vSeq.push_back(strTarget);
            iIndex++;
        }
    }

    if(vSeq.empty())
        return "";

    //Check those sequence
    int iStrLen = vSeq[0].length();
    //(1) we first cut the starting posiiton and ending position
    //    *delete the head part and the ending part if do not support by enough reads
       /// Beginning -->
    int iStart = 0;
    map<char, int> mpChar;
    for(int i=0; i<iStrLen; i++)
    {
        if((*(vSeq.end() - 1))[i] == '*')
        {
            iStart = i;
            break;
        }
        else
        {
            mpChar.clear();
            for(vector<string>::iterator itr = vSeq.begin(); itr != vSeq.end() - 1; itr++) // do not check the status row
            {
                mpChar[(*itr)[i]]++;
            }
            int iMax = 0;
            char cMax;
            for(map<char, int>::iterator itr = mpChar.begin(); itr != mpChar.end(); itr++)
            {
                if(iMax <= itr->second)
                {
                    cMax = itr->first;
                    iMax = itr->second;
                }
            }
            if(cMax == '-')
                continue;
            else
            {
                iStart = i;
                break;
            }
        }
    }

       /// Ending -->
    int iEnd = 0;
    mpChar.clear();
    for(int i=iStrLen-1; i>=0; i--)
    {
        if((*(vSeq.end() - 1))[i] == '*')
        {
            iEnd = i;
            break;
        }
        else
        {
            mpChar.clear();
            for(vector<string>::iterator itr = vSeq.begin(); itr != vSeq.end() - 1; itr++) // do not check the status row
            {
                mpChar[(*itr)[i]]++;
            }
            int iMax = 0;
            char cMax;
            for(map<char, int>::iterator itr = mpChar.begin(); itr != mpChar.end(); itr++)
            {
                if(iMax <= itr->second)
                {
                    cMax = itr->first;
                    iMax = itr->second;
                }
            }
            if(cMax == '-')
                continue;
            else
            {
                iEnd = i;
                break;
            }
        }
    }

    //(2) get the most frequent sequence in the valide range
    string strSumSeq = "";
    for(int i=iStart; i<=iEnd; i++)
    {
        if((*(vSeq.end() - 1))[i] == '*')
        {
            strSumSeq += vSeq[0][i];
            continue;
        }
        else
        {
            mpChar.clear();
            for(vector<string>::iterator itr = vSeq.begin(); itr != vSeq.end() - 1; itr++) // do not check the status row
            {
                if((*itr)[i] != '-') // only record the valid sequence
                {
                    mpChar[(*itr)[i]]++;
                }

            }
            int iMax = 0;
            char cMax;
            for(map<char, int>::iterator itr = mpChar.begin(); itr != mpChar.end(); itr++)
            {
                //Keep the first maximum value
                if(iMax < itr->second)
                {
                    cMax = itr->first;
                    iMax = itr->second;
                }
            }
            strSumSeq += cMax;
        }
    }
    //Now we have the correct one
    ifsAlignedFile.close();
    return strSumSeq;
}


