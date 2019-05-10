#include "clsblast.h"
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include "clsbasealgorithm.h"
#include "clsfastareader.h"
#include <algorithm>

ClsBlast::ClsBlast()
{
    m_vBlastResultSum.clear();
}

int ClsBlast::TwoSeqBlast(string& strQueryPath, string& strSubjectPath, bool bShowFlank)
{
    St_BlastResult stBlastResult;
    m_vBlastResultSum.clear();
    GetBlastResult(strQueryPath, strSubjectPath, stBlastResult);

    //取得里面align最好的一组然后进行相应的记录
    //1: Get the max number of alignment
    int iMaxAlignNum = -1;
    float fIdentifyRatio = 0;
    float fGapRatio = 0;
    bool bRC = false; // RC: Reverse Complementary!!!
    vector<St_BlastAlignUnit>::iterator itrBestAlgn = stBlastResult.vBlastAlign.begin();
    for(vector<St_BlastAlignUnit>::iterator itr = stBlastResult.vBlastAlign.begin();
        itr != stBlastResult.vBlastAlign.end(); itr++)
    {
        if(iMaxAlignNum < itr->iAlignNum)
        {
            iMaxAlignNum = itr->iAlignNum;
            fIdentifyRatio = itr->fIdentifyRatio;
            fGapRatio = itr->fGapRatio;
            if(itr->iStartPos < itr->iEndPos)
                bRC = false;
            else
                bRC = true;
            itrBestAlgn = itr;
        }
    }
    cout << "The Max Alignment Number is: " << IntToStr(iMaxAlignNum) << endl;
    cout << "Identities Ratio (Total) is: " << GetRatio(fIdentifyRatio) << endl;
    cout << "Gap Ratio is               : " << GetRatio(fGapRatio) << endl;
    cout << "Is Reverse Complementary: ";
    if(bRC)
        cout << "YES" << endl;
    else
        cout << "NO" << endl;

    if(bShowFlank)
    {
        //Extract the left 50 characters and the right 50 characters  --> Go
        if(iMaxAlignNum > 0) // It means we have got it
        {
            ClsFastaReader* pFaReader = new ClsFastaReader();
            //--->
            vector<St_Fasta> vFasta;
            pFaReader->ReadFastaRegular(strSubjectPath, vFasta);
            for(vector<St_Fasta>::iterator itr = vFasta.begin(); itr != vFasta.end(); itr++)
            {
                if(itr->strName == itrBestAlgn->strName)
                {
                    //if(itr->strSeq.find('Y') != string::npos)
                    //    cout << "WE GOT YYYY!!!" << endl;
                    const int FLANKLEN = 50;
                    //Big <--> Right Extend  and Small <--> Left Extend
                    string strLeftFlank = "";
                    string strRightFlank = "";
                    if(bRC)
                    {
                        //This is the left flank --> Blast start from 1
                        int iStart = itrBestAlgn->iStartPos;
                        int iLen = itr->strSeq.length() - iStart > FLANKLEN ? FLANKLEN :
                                                                              itr->strSeq.length() - iStart;
                        if(iLen > 0)
                        {
                            // the reason of GetReverse is that this is the RC Case
                            strLeftFlank = GetReverse(itr->strSeq.substr(iStart, iLen));
                        }

                        //This is the right flank
                        iLen = itrBestAlgn->iEndPos - FLANKLEN > 0 ? FLANKLEN : 0;
                        if(iLen > 0)
                        {
                            iStart = itrBestAlgn->iEndPos - iLen;
                            strRightFlank = GetReverse(itr->strSeq.substr(iStart, iLen));
                        }
                    }
                    else // Normal Orientation
                    {
                        //For Left Flank
                        int iLen = itrBestAlgn->iStartPos - FLANKLEN > 0 ? FLANKLEN : 0;
                        int iStart = 0;
                        if(iLen > 0)
                        {
                            iStart = itrBestAlgn->iStartPos - iLen;
                            strLeftFlank = itr->strSeq.substr(iStart, iLen);
                        }

                        //For Right Flanl
                        iStart = itrBestAlgn->iEndPos;
                        iLen = itr->strSeq.length() - iStart > FLANKLEN ? FLANKLEN :
                                                                        itr->strSeq.length() - iStart;
                        strRightFlank = itr->strSeq.substr(iStart, iLen);
                    }
                    //Cout Left Flank and Right Flank
                    //Left Flank
                    cout << "Left Flank: -->" << endl;
                    cout << strLeftFlank << endl;
                    //Right Flank
                    cout << "Right Flank: -->" << endl;
                    cout << strRightFlank << endl;

                    break;
                }
            }
            delete pFaReader;
            pFaReader = NULL;
            //<---
        }
    }
    cout << "===============================" << endl;
    return iMaxAlignNum;
}

bool sort_blastAlignUnit_func(St_BlastAlignUnit stUnit1, St_BlastAlignUnit stUnit2)
{
    if(stUnit1.iAlignNum > stUnit2.iAlignNum)
        return true;
    else
        return false;
}

void ClsBlast::GetTopNResultsFromTwoSeqBlast(string& strQueryPath, string& strSubjectPath,
                                             vector<St_BlastAlignUnit>& vBlastAlign, int iTopN)
{
    St_BlastResult stBlastResult;
    m_vBlastResultSum.clear();
    int iResultCode = GetBlastResult(strQueryPath, strSubjectPath, stBlastResult);
    if(iResultCode != 0)
        return;

    if(stBlastResult.vBlastAlign.empty())
    {
        cout << "No Alignment Result." << endl;
        return;
    }

    //Sort Result --> From longest to shortest
    sort(stBlastResult.vBlastAlign.begin(), stBlastResult.vBlastAlign.end(), sort_blastAlignUnit_func);

    if(iTopN > stBlastResult.vBlastAlign.size())
        iTopN = stBlastResult.vBlastAlign.size();

    vBlastAlign.clear();
    vBlastAlign.insert(vBlastAlign.end(), stBlastResult.vBlastAlign.begin(), stBlastResult.vBlastAlign.begin()+iTopN);
}

int ClsBlast::GetBlastResult(string& strQueryPath, string& strSubjectPath, St_BlastResult& stBlastResult)
{
    //Check Path validation
    if(access(strQueryPath.c_str(), 0) != 0)
    {
        cout << "Error: \"QueryPath\" is not existed!" << endl;
        return -1;
    }
    if(access(strQueryPath.c_str(), 0) != 0)
    {
        cout << "Error: \"strSubjectPath\" is not existed!" << endl;
        return -1;
    }

    //Make Blast
    /// use default normal size
    string strBlastResultPath = TwoFastaFileAlign(strQueryPath, strSubjectPath, 1); /// Try to use small size
    stBlastResult = ParseResult(strBlastResultPath, m_vBlastResultSum, false);
    return 0;
}

string ClsBlast::CreatFaFile(string strName, string& strSeq, En_FastaType enFastaType)
{
    if(strSeq == "")
        return "";

    string strFolderPath = ::GetHigherFolderPath(::GetCurExeFolderPath()) +
                           "TempFile/";
    string strFileName = "";
    switch(enFastaType)
    {
        case ftQuery:
            strFileName = "querySeq.fa";
            break;
        case ftRef:
            strFileName = "refSeq.fa";
            break;
        default:
            strFileName = "querySeq.fa";
            break;
    }
    string strSingleSegFa = strFolderPath + strFileName; //"querySeq.fa";
    ofstream ofs;
    ofs.open(strSingleSegFa.c_str());
    if(!ofs.is_open())
    {
        cout << "create file failed" << endl;
        return "";
    }
    //Save fasta file
    ofs << ">" << strName << endl;
    ::DisplayString(ofs, strSeq, 100);
    ofs.close();
    return strSingleSegFa;
}

//Query Sequence: 这个里面我们默认值存储一条序列记录
//Target Sequence: 这个里面我们可以存储多条序列记录，从而进行比较方便的群体性分析
string ClsBlast::TwoFastaFileAlign(string strQuerySeqPath, string strTargetSeqPath, int iSeqSizeType)
{
    if( access(strQuerySeqPath.c_str(), 0) != 0 ||
        access(strTargetSeqPath.c_str(), 0) != 0 )
    {
        cout << "Error: one of input paths are not existed!" << endl;
        return "";
    }
    string strCmd = "";
    string strRootPath = GetHigherFolderPath(get_current_dir_name());
    string strBlastnPath = strRootPath + "../../../ShareLibrary/Blast/blastn";
    string strOutputPath = strRootPath + "../../../ShareLibrary/Blast/log.ini";
    switch(iSeqSizeType)
    {
        case 0: //Normal Size
            strCmd = strBlastnPath + " -query " + strQuerySeqPath + " -subject " + strTargetSeqPath +
                    " > " + strOutputPath;
            break;
        case 1: //Small Size
            strCmd = strBlastnPath + " -word_size 15 -evalue 100 " + " -dust no " + //"-gapopen 10 -gapextend 20 -reward 1 -dust no "
                     " -query " + strQuerySeqPath + " -subject " + strTargetSeqPath +
                     " > " + strOutputPath;
            break;
        case 2: //Large Size
            break;
        default:
            break;
    }
    //Run the the command line
    try
    {
        system(strCmd.c_str());
    }
    catch(...)
    {
        cout << "Error: command line failed to be executed!" << endl;
        return "";
    }
    return strOutputPath;
}

//这个parse的文件，是基于默认参数(主要是条目显示方面的)的生成文件的解析
St_BlastResult ClsBlast::ParseResult(string strRstPath, vector<St_BlastResult>& vBlastResult,
                                     bool bRecord) // Rst-->Result
{
    St_BlastAlignUnit stBlastUnit;
    St_BlastResult stBlastResult;

    //Parse the alignment result
    if( access(strRstPath.c_str(), 0) != 0)
    {
        cout << "Error: File Path is not existed!" << endl;
        return stBlastResult;
    }

    //let record every fragment, and then try to make the analysis for each fragment
    ifstream ifs;
    ifs.open(strRstPath.c_str());
    bool bFindHint = true;
    bool bFirstTimeFindSbjct = true;
    bool bFirstTimeFindQuery = true;
    string strLine = "";
    while(!ifs.eof())
    {
        getline(ifs, strLine);
        if(strLine.find("Subject") != string::npos)
        {
            int iStart = strLine.find(" ") + 1;
            int iCount = strLine.find(" ", iStart) - iStart;
            stBlastUnit.strName = strLine.substr(iStart, iCount);
            continue;
        }
        if(strLine.find("No hits found") != string::npos)
        {
            bFindHint = false;
            continue;
        }
        if(strLine.find("Identities") != string::npos)
        {
            //Identities is the new start of another part in the same reference sequence
            if(!bFirstTimeFindSbjct)
            {
                stBlastResult.vBlastAlign.push_back(stBlastUnit);
                stBlastUnit.Clear();
                bFirstTimeFindSbjct = true;
                bFirstTimeFindQuery = true;
            }
            //(1) Get Total Align Num
            int iStart = strLine.find('/') + 1;
            int iCount = strLine.find('(') - iStart - 1;
            stBlastUnit.iAlignNum = atoi(strLine.substr(iStart, iCount).c_str());

            //(2) Get Identities Percentage: -> substitution + indel
            iStart = strLine.find('(') + 1;
            iCount = strLine.find('%') - iStart;
            stBlastUnit.fIdentifyRatio = (float)atoi(strLine.substr(iStart, iCount).c_str()) / 100;

            //(3) Get Gap Percentage
            iStart = strLine.rfind('(') + 1;
            iCount = strLine.rfind('%') - iStart;
            stBlastUnit.fGapRatio = (float)atoi(strLine.substr(iStart, iCount).c_str()) / 100;

            continue;
        }
        if(strLine.find("Query ") != string::npos)
        {
            int iStart = strLine.find("Query") + 4 + 3; // the length of Sbjct
            int iStartCount = strLine.find(" ", iStart) - iStart;
            int iEndPos = strLine.rfind(" ") + 1;
            int iEndCount = strLine.length() - iEndPos;

            //这里减去1的原因是因为，他们的计数是从1开始的
            if(bFirstTimeFindQuery)
            {
                stBlastUnit.iStartPosQuery = atoi(strLine.substr(iStart, iStartCount).c_str())-1;
                bFirstTimeFindQuery = false;
                stBlastUnit.iEndPosQuery = atoi(strLine.substr(iEndPos, iEndCount).c_str())-1;
            }
            else // not the first one
                stBlastUnit.iEndPosQuery = atoi(strLine.substr(iEndPos, iEndCount).c_str())-1;
            continue;
        }
        if(strLine.find("Sbjct") != string::npos)
        {
            int iStart = strLine.find("Sbjct") + 4 + 3; // the length of Sbjct
            int iStartCount = strLine.find(" ", iStart) - iStart;
            int iEndPos = strLine.rfind(" ") + 1;
            int iEndCount = strLine.length() - iEndPos;

            //这里减去1的原因是因为，他们的计数是从1开始的
            if(bFirstTimeFindSbjct)
            {
                stBlastUnit.iStartPos = atoi(strLine.substr(iStart, iStartCount).c_str())-1;
                bFirstTimeFindSbjct = false;
                stBlastUnit.iEndPos = atoi(strLine.substr(iEndPos, iEndCount).c_str())-1;
            }
            else // not the first one
                stBlastUnit.iEndPos = atoi(strLine.substr(iEndPos, iEndCount).c_str())-1;
            continue;
        }
        if(strLine.find("Effective") != string::npos) // as new unit end -->end
        {
            if(bFindHint)
                stBlastResult.vBlastAlign.push_back(stBlastUnit);
            bFirstTimeFindSbjct = true;
            bFirstTimeFindQuery = true;
            bFindHint = true;
            stBlastUnit.Clear();
            continue;
        }
    }
    ifs.close();
    if(bRecord)
        vBlastResult.push_back(stBlastResult);

    return stBlastResult;
}
