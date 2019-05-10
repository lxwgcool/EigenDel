#include <iostream>
#include "clsparsebam.h"
#include "clsconfig.h"
#include <algorithm>
#include "../../../ShareLibrary/clsbwa.h"
#include "clssvdeldetect.h"
#include "clsdebug.h"
#include "clslearning.h"
#include <map>
#include <stdexcept>
#include <memory>
#include <array>
#include "clscomparison.h"
#include "unistd.h"
#include <algorithm>

using namespace std;

#define RUN_IN_CLUSTER
#define MULTITHREAD
#define GLOBALWHOLEGENO

//const char* SAMPLENAME = "NA12890";
//const int CHROM = 11;//20;
//Find Del Chrom By Chrom -->
struct St_FindSvDEL
{
    //vector<St_Fasta>* pVFasta;
    St_ChromSV* pSTDChromSV;

    St_FinalResult stFinalResult;
    string strChromName;
    string strSampleName;
    string strBamFilePath;
    string strPythonCodePath;    

    //--> this is for the case of using the statistic value of the whole genome
    //    (Not chromosome by chromosome respectively)
    float fAvgInsertSize;
    float fSTDInsertSize;
    float fAvgDepth;

    //--> This is for target reads
    St_ChromDiscordantReads* pChromDR;
    St_ChromBorderSCR* pChromSCR;

    bool bTargetReadsGood;

    St_FindSvDEL()
    {
        Reset();
    }

    void Init(//vector<St_Fasta>* pV1,
              St_ChromSV* pV2,
              string strV3, string strV4,
              string strV5,
              string strV6,              
              St_ChromBorderSCR* pV8, St_ChromDiscordantReads* pV9,
              bool bV10)
    {
        //pVFasta = pV1;
        pSTDChromSV = pV2;

        strChromName = strV3;
        strSampleName = strV4;
        strBamFilePath = strV5;
        strPythonCodePath = strV6;        

        pChromSCR = pV8;
        pChromDR = pV9;

        bTargetReadsGood = bV10;

        stFinalResult.Clear();
    }

    void Reset()
    {
        //pVFasta = NULL;
        pSTDChromSV = NULL;

        strChromName = "";
        strSampleName = "";
        //strBamFilePath = "";
        strPythonCodePath = "";       

        stFinalResult.Clear();

        fAvgInsertSize = -1;
        fSTDInsertSize = -1;
        fAvgDepth = -1;

        pChromSCR = NULL;
        pChromDR = NULL;

        bTargetReadsGood = false;
    }

    void SetGlobalGenoValue(float fV1, float fV2, float fV3)
    {
        this->fAvgInsertSize = fV1;
        this->fSTDInsertSize = fV2;
        this->fAvgDepth = fV3;
    }

    void PrintSCRAndDRInfo()
    {
        cout << "SCR: (" << pChromSCR->strChrom << ", " << to_string(pChromSCR->vSCR.size()) << ")"
             << "\t"
             << "DR: (" << pChromDR->strChrom << ", " << to_string(pChromDR->vDR.size()) << ")" << endl;
    }
};

void FindSvDel(string strConfigFile);
void Comparision(string strConfigFile);

void PrepareTargetReads(ClsParseBam* pClsParseBam, string strBamFile,
                        vector<St_Fasta>& vFasta, vector<St_ChromBorderSCR>& vChromBorderSCR,
                        vector<St_ChromDiscordantReads>& vChromDR);
                        //vector<St_BorderSCR>& vBorderSCR, vector<St_DiscordantReads>& vDiscdRreads);

void FindSvDelByChromosome(ClsParseBam* pClsParseBam, string strBamFile,
                           string strPythonCode,
                           string strChrom, string strSampleName,
                           vector<St_BorderSCR>& vBorderSCR, vector<St_DiscordantReads>& vDiscdRreads,
                           vector<St_SV>& vSvDEL, St_FinalResult& stFinalResult, ofstream& ofs);

void PrepareGroup(ClsSvDelDetect* pSvDelDetect, ClsDebug* pDebug,
                  vector<St_BorderSCR>& vBorderSCR, vector<St_DiscordantReads>& vDiscdRreads,
                  vector<St_SV>& vSvDEL, vector<St_DRGroup>& vDRGroup, vector<St_SCRGroup>& vSCRGroup,
                  string strBamFile, string strChrom);
void DrawImage();


void* FindSvDelForSingleChrom(void* pValue);
void FindDelByChrom(St_Config& stConfig, vector<St_Fasta>& vFasta, vector<St_ChromSV>& vChromSvDEL);
void SaveToFileEigenDelResult();
void PrintEigenDelResult(vector<St_FindSvDEL>& vFindSvDEL);

bool GetCurChromDRAndSCR(vector<St_ChromDiscordantReads>& vChromDR,
                         vector<St_ChromBorderSCR>& vChromBorderSCR,
                         string strSTDSvDelChromName,
                         int& iSCRIndex,
                         int& iDRIndex);


int main(int argc, char *argv[])
{           
    if(argc > 2)
    {
        cout << "To many arguments!" << endl;
        return 1;
    }

    if(argc < 2)
    {
        cout << "Please give the path of config file!" << endl;
        return 1;
    }

//    int iValue = ceil(.25 * 7);
//    iValue = floor(.25 * 7);
//    cout << "I am good" << endl;

    //--> For testing *****************
//    ClsSvDelDetect* pSvDelDetect = new ClsSvDelDetect();

//    St_DRGroup stDRGroup;
//    stDRGroup.iLeftBoundary = 16581704; //16582093 - 389;
//    stDRGroup.iRightBoundary = 16582131; //16582183 - 52;

//    stDRGroup.iStart = 16582093 - 300;
//    stDRGroup.iEnd = 16582183 - 74;

//    string strBamFile = "/home/xin/lxwg/WorkStudio/Prototype/SV_Draft/Data/SvBamFile/NA12878.chrom11.ILLUMINA.bwa.CEU.low_coverage.20121211.bam";
//    pSvDelDetect->GetDepthOfRange(stDRGroup, strBamFile);

//    delete pSvDelDetect;
//    pSvDelDetect = NULL;
    //<-- ************************

    FindSvDel(argv[1]);

    //Comparision(argv[1]);

    return 0;
}

bool sort_sv_small2large_func(St_SV stSv1, St_SV stSv2)
{
    if((stSv1.iEnd - stSv1.iPos) < (stSv2.iEnd - stSv2.iPos))
        return true;
    else
        return false;
};

void FindSvDel(string strConfigFile)
{
    //Read Config File
    ClsConfig* pConfig = new ClsConfig();
    St_Config stConfig;
    pConfig->ReadConfig(stConfig, strConfigFile.c_str());
    delete pConfig;
    pConfig = NULL;

    //--> Read Reference and get the corresponding ChromName and ChromName for whole genome --> Do tomorrow
    ClsFastaReader* pFastaReader = new ClsFastaReader();
    vector<St_Fasta> vFasta;
    pFastaReader->ReadFastaRegular(stConfig.strRef, vFasta, false); // Get the name list
    delete pFastaReader;
    pFastaReader = NULL;
    //<--

    //Read VCF File
    ClsVcf1000Genome* pClsVcf = new ClsVcf1000Genome();
    cout << "ParseVcf --> " << endl;
    pClsVcf->ParseVcf(stConfig.strVcf, stConfig.strSampleName);

    //vector<St_SV> vSvDEL; // only pick out deletion
    vector<St_ChromSV> vChromSvDEL;
    cout << "GetDeletion --> " << endl;
    pClsVcf->GetDeletion(vChromSvDEL, stConfig.strSampleName, stConfig.strChrom);
    //-->Out put Size:
    cout << "vChromSvDEL size: " << to_string(vChromSvDEL.size()) << endl;
    cout << "==== Detail ====" << endl;
    for(vector<St_ChromSV>::iterator itr = vChromSvDEL.begin(); itr != vChromSvDEL.end(); itr++)
    {
        cout << itr->strChrom << ": " << to_string(itr->vSv.size()) << endl;
    }
    cout << "================" << endl;

    delete pClsVcf;
    pClsVcf = NULL;

//    //Check the coverage --> ---------------------
//    cout << "Sort SV First" << endl;
//    sort(vSvDEL.begin(), vSvDEL.end(), sort_sv_small2large_func);
//    cout << "Sort Finished" << endl;
//    float fThreshold1 = 3;
//    int iNum1 = 0;
//    float fThreshold2 = 2;
//    int iNum2 = 0;
//    float fThreshold3 = 4;
//    int iNum3 = 0;


//    for(vector<St_SV>::iterator itr = vSvDEL.begin(); itr != vSvDEL.end(); itr++)
//    {
//        string strCmd = (string)"samtools depth -a -r " +
//                        "11:" + to_string(itr->iPos) + "-" + to_string(itr->iEnd) + " " + stConfig.strBamFile +
//                        " | awk '{sum+=$3} END {if(sum==\"\") {print 0} else {print sum/NR}}'";

//        string strResult = exec(strCmd.c_str());
//        float fCoverage = stof(strResult.c_str());

//        cout << strResult;

//        if(fCoverage <= fThreshold1)
//           iNum1++;

//        if(fCoverage <= fThreshold2)
//           iNum2++;

//        if(fCoverage <= fThreshold3)
//           iNum3++;


//        strCmd = (string)"samtools depth -a -r " +
//                 "11:" + to_string(itr->iPos - 300) + "-" + to_string(itr->iPos - 150) + " " + stConfig.strBamFile +
//                 " | awk '{sum+=$3} END {if(sum==\"\") {print 0} else {print sum/NR}}'";
//        cout << "\t" << exec(strCmd.c_str());

//        strCmd = (string)"samtools depth -a -r " +
//                 "11:" + to_string(itr->iEnd + 150) + "-" + to_string(itr->iEnd + 300) + " " + stConfig.strBamFile +
//                 " | awk '{sum+=$3} END {if(sum==\"\") {print 0} else {print sum/NR}}'";

//        cout << "\t" << exec(strCmd.c_str());
    //        cout << "\t" << "<" << to_string(itr->iPos) << ", " << to_string(itr->iEnd) << ">" << endl;
//    }
//    cout << IntToStr(iNum1) << " \t" << IntToStr(iNum2) << "\t" << IntToStr(iNum3) << endl;
//    //<-------------------
//    return;

#ifdef MULTITHREAD
    ////-------> Run Multiple Thread -->  ---------------------------------
    FindDelByChrom(stConfig, vFasta, vChromSvDEL);
    ////<-------
#else
    //Calc Insert Size
    ClsParseBam* pClsParseBam = new ClsParseBam();

#ifdef RUN_IN_CLUSTER
    pClsParseBam->CalcInsertSize(stConfig.strBamFile);
#else
    pClsParseBam->SetMeanInsertSize(410.026265);//(379);
    pClsParseBam->SetStdDevInsertSize(91.85553);//(38);
#endif

#ifdef RUN_IN_CLUSTER
    pClsParseBam->CalcAvgDepth(stConfig.strBamFile);
#else
    pClsParseBam->SetAvgDepth(7.0);
#endif

    cout << "Mean Insert Size         : " << FloatToStr(pClsParseBam->GetMeanInsertSize()) << endl;
    cout << "Std Deviation Insert Size: " << FloatToStr(pClsParseBam->GetStdDevInsertSize()) << endl;
    cout << "Avg Depth                : " << FloatToStr(pClsParseBam->GetAvgDepth()) << endl;

    //Prepare Target Reads: discordant reads and clipped reads
    vector<St_ChromBorderSCR> vChromBorderSCR;
    vector<St_ChromDiscordantReads> vChromDR;
    //vector<St_BorderSCR> vBorderSCR;
    //vector<St_DiscordantReads> vDiscdRreads;
    PrepareTargetReads(pClsParseBam, stConfig.strBamFile, vFasta, vChromBorderSCR, vChromDR);
    cout << "Finish PrepareTargetReads" << endl;

    //Get Result for each different chromosome respectively ->
    St_ChromDiscordantReads* pChromDR = NULL;
    St_ChromBorderSCR* pChromSCR = NULL;
    St_ChromSV* pChromSvDel = NULL;

    vector<St_FinalResult> vFinalResult;
    St_FinalResult stFinalResult;
    ofstream ofsAll;
    ofsAll.open("FindSvDEL_ALL.txt");
    for(vector<St_ChromDiscordantReads>::iterator itrChromDR = vChromDR.begin();
        itrChromDR != vChromDR.end(); itrChromDR++)
    {
        //(1)Get Discordant Reads
        string strChrom = itrChromDR->strChrom;
        pChromDR = &(*itrChromDR);

        //(2)Get Border Softclip Reads
        for(vector<St_ChromBorderSCR>::iterator itrChromSCR = vChromBorderSCR.begin();
            itrChromSCR != vChromBorderSCR.end(); itrChromSCR++)
        {
            if(itrChromSCR->strChrom == strChrom)
            {
                pChromSCR = &(*itrChromSCR);
                break;
            }
        }

        //(3)Get related standard SV Deletion -->
        for(vector<St_ChromSV>::iterator itrSvDel = vChromSvDEL.begin();
            itrSvDel != vChromSvDEL.end(); itrSvDel++)
        {
            if(itrSvDel->strChrom == strChrom)
            {
                pChromSvDel = &(*itrSvDel);
                break;
            }
        }

        if( pChromDR == NULL ||
            pChromSCR == NULL ||
            pChromSvDel == NULL)
        {}
        else
        {
            stFinalResult.strChrom = strChrom;
            cout << "Do FindSvDelByChromosome For -- " << strChrom << endl;
            FindSvDelByChromosome(pClsParseBam, stConfig.strBamFile,
                                  stConfig.strPythonCode, strChrom,
                                  pChromSCR->vSCR, pChromDR->vDR, pChromSvDel->vSv,
                                  stFinalResult, ofsAll);

            vFinalResult.push_back(stFinalResult);
            stFinalResult.Clear();
        }
        pChromDR = NULL;
        pChromSCR = NULL;
        pChromSvDel = NULL;
    }
    //<--
    ofsAll.close();

    cout << " ====== Final Sum ========" << endl;
    int iPosSum = 0;
    int iNegSum = 0;
    for(vector<St_FinalResult>::iterator itr = vFinalResult.begin();
        itr != vFinalResult.end(); itr++)
    {
        cout << itr->strChrom << ": " << to_string(itr->iPositiveNum)
             << "\t" << to_string(itr->iNegativeNum) << endl;

        iPosSum += itr->iPositiveNum;
        iNegSum += itr->iNegativeNum;
    }

    cout << "--------------" << endl;
    cout << "Total Positive: " << to_string(iPosSum) << endl;
    cout << "Total Negative: " << to_string(iNegSum) << endl;
    cout << "Total Sum     : " << to_string(iPosSum + iNegSum) << endl;

    //release the memory
    delete pClsParseBam;
    pClsParseBam = NULL;
#endif
}

//-->Do multiple thread
void FindDelByChrom(St_Config& stConfig, vector<St_Fasta>& vFasta, vector<St_ChromSV>& vChromSvDEL)
{
    //--> Collect Void Parameter for each thread
    vector<St_FindSvDEL> vFindSvDEL;
    vFindSvDEL.resize(vChromSvDEL.size());
    if(vFindSvDEL.empty())
        return;

    if(stConfig.strChrom != "")
    {
        if(vChromSvDEL.size() != 1)
            return;
    }

    //Create output folder
    string strCmd = "";
    //------->  For temporary file from EigenDel
    strCmd = "mkdir -p ./Output";
    ///1: Create folder
    system(strCmd.c_str());
    ///2: Clean folder
    strCmd = "rm -rf ./Output/*";
    system(strCmd.c_str());
    //------> This is for the file created by Python (sklearn)
    strCmd = "mkdir -p ./Pics";
    ///1: Create folder
    system(strCmd.c_str());
    ///2: Clean folder
    strCmd = "rm -rf ./Pics/*";
    system(strCmd.c_str());

//---> Check if we need to use the global value of the whole genome -->
#ifdef GLOBALWHOLEGENO
    ClsParseBam* pClsParseBam = new ClsParseBam();
    //1: Get global insert size
    pClsParseBam->CalcInsertSize(stConfig.strBamFile, stConfig.strPicard);
    cout << "Calc Insert Size Done!" << endl;
    //2: Get global depth
    pClsParseBam->CalcAvgDepth(stConfig.strBamFile);
    cout << "Calc Depth Done!" << endl;

//    //Debug
//    pClsParseBam->SetMeanInsertSize(410.033131);
//    pClsParseBam->SetStdDevInsertSize(92.985886);
//    pClsParseBam->SetAvgDepth(7.84);


    //3: extract those three types of value
    float fAvgInsertSize = pClsParseBam->GetMeanInsertSize();
    float fSTDInsertSize = pClsParseBam->GetStdDevInsertSize();
    float fAvgDepth = pClsParseBam->GetAvgDepth();

#endif
//<--- End    
    vector<St_ChromBorderSCR> vChromBorderSCR;
    vector<St_ChromDiscordantReads> vChromDR;
    PrepareTargetReads(pClsParseBam, stConfig.strBamFile, vFasta, vChromBorderSCR, vChromDR);
    cout << "vChromBorderSCR Size: " << to_string(vChromBorderSCR.size()) << endl;
    cout << "vChromDR Size       : " << to_string(vChromDR.size()) << endl;
    delete pClsParseBam;
    pClsParseBam = NULL;

    if(stConfig.strChrom != "") // i only want to do one chromosome -->
    {        
        int iSCRIndex = -1;
        int iDRIndex = -1;
        bool bGetReads = GetCurChromDRAndSCR(vChromDR, vChromBorderSCR, vChromSvDEL[0].strChrom, iSCRIndex, iDRIndex);
        if(!bGetReads)
        {
            cout << "Fail to find target reads!" << endl;
        }
        else
        {
            vFindSvDEL.begin()->Init(&vChromSvDEL[0],
                                     stConfig.strChrom, stConfig.strSampleName, stConfig.strBamFile,
                                     stConfig.strPythonCode,
                                     &vChromBorderSCR[iSCRIndex], &vChromDR[iDRIndex], bGetReads);
            #ifdef GLOBALWHOLEGENO
                vFindSvDEL.begin()->SetGlobalGenoValue(fAvgInsertSize, fSTDInsertSize, fAvgDepth);
            #endif
        }
    }
    else
    {
        vector<St_FindSvDEL>::iterator itrFindSvDEL = vFindSvDEL.begin();
        int i = 0;
        for(vector<St_ChromSV>::iterator itrChromSvDEL = vChromSvDEL.begin();
            itrChromSvDEL != vChromSvDEL.end(); itrChromSvDEL++, itrFindSvDEL++, i++)
        {
            //The Major job is getting Bam File -->
            string strChrom = itrChromSvDEL->strChrom;
            int iLen = stConfig.strBamFile.find(".bam");
            string strBamFile = stConfig.strBamFile.substr(0, iLen) + ".REF_" +
                                strChrom + ".bam";
            int iSCRIndex = -1;
            int iDRIndex = -1;
            bool bTargetReadsGood = GetCurChromDRAndSCR(vChromDR, vChromBorderSCR, strChrom, iSCRIndex, iDRIndex);
            if(!bTargetReadsGood)
            {
                cout << "Wrong: " << to_string(i) << endl;
            }
            else
            {
                itrFindSvDEL->Init(&(*itrChromSvDEL),
                                   strChrom, stConfig.strSampleName, strBamFile,
                                   stConfig.strPythonCode,
                                   &vChromBorderSCR[iSCRIndex], &vChromDR[iDRIndex], bTargetReadsGood);

                //itrFindSvDEL->PrintSCRAndDRInfo();

            #ifdef GLOBALWHOLEGENO
                itrFindSvDEL->SetGlobalGenoValue(fAvgInsertSize, fSTDInsertSize, fAvgDepth);
            #endif                
            }
        }
        cout << "Set itrFindSvDEL Finished" << endl;
    }

    int iThreadsNum = vChromSvDEL.size();
    if(iThreadsNum == 0)
    {
        cout << "OMG!! 0 THREAD!!!" << endl;
        return;
    }

    //Step 1: Define how may threads you want to use
    pthread_t threads[iThreadsNum];
    pthread_attr_t attr;

    //Step 2: Initialize and set thread joinable
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    //Step 3: Do the main body of multiple threads function
    for(int i = 0; i < iThreadsNum; i++)
    {
        //cout << "main() : creating thread, " << i << endl;
        if(!vFindSvDEL[i].bTargetReadsGood)
        {
            cout << "Thread " << to_string(i) << " Failed to be creared" << endl;
            continue;
        }

        int rc = pthread_create(&threads[i], NULL, FindSvDelForSingleChrom, (void *)&vFindSvDEL[i]);

        if(rc)
        {
            cout << "Error:unable to create thread," << rc << endl;
            exit(-1);
        }
    }

    //4: free attribute and wait for the other threads
    pthread_attr_destroy(&attr);

    void *status;
    //5: Join those threads together --> to make sure the remaining part of the code will be run only after those threads been finished
    cout << "Start Joint Threads -->" << endl;
    for(int i = 0; i < iThreadsNum; i++)
    {
        if(!vFindSvDEL[i].bTargetReadsGood)
            continue;

        int rc = pthread_join(threads[i], &status);
        if (rc)
        {
            cout << "Error:unable to join," << rc << endl;
            exit(-1);
        }

        //cout << "Main: completed thread id :" << i ;
        //cout << "  exiting with status :" << status << endl;
    }
    cout << "Joint Threads End <--" << endl;

    //6: Do the remaning thing    
//    cout << " ====== Final Sum ========" << endl;
//    int iPosSum = 0;
//    int iNegSum = 0;
//    for(vector<St_FindSvDEL>::iterator itr = vFindSvDEL.begin();
//        itr != vFindSvDEL.end(); itr++)
//    {
//        cout << itr->stFinalResult.strChrom << ": " << to_string(itr->stFinalResult.iPositiveNum)
//             << "\t" << to_string(itr->stFinalResult.iNegativeNum) << endl;

//        iPosSum += itr->stFinalResult.iPositiveNum;
//        iNegSum += itr->stFinalResult.iNegativeNum;
//    }

//    cout << "--------------" << endl;
//    cout << "Total Positive: " << to_string(iPosSum) << endl;
//    cout << "Total Negative: " << to_string(iNegSum) << endl;
//    cout << "Total Sum     : " << to_string(iPosSum + iNegSum) << endl;

//    cout << endl << "ALL SET!!!" << endl;

    SaveToFileEigenDelResult(); // This is for the final candidate found by EigenDel
    PrintEigenDelResult(vFindSvDEL); // This is for comparison result
}

bool GetCurChromDRAndSCR(vector<St_ChromDiscordantReads>& vChromDR,
                         vector<St_ChromBorderSCR>& vChromBorderSCR,
                         string strSTDSvDelChromName,
                         int& iSCRIndex,
                         int& iDRIndex)
{
    iSCRIndex = -1;
    iDRIndex = -1;

    //(1)Get Border Softclip Reads
    for(vector<St_ChromBorderSCR>::iterator itrChromSCR = vChromBorderSCR.begin();
        itrChromSCR != vChromBorderSCR.end(); itrChromSCR++)
    {
        if(itrChromSCR->strChrom == strSTDSvDelChromName)
        {
            iSCRIndex = itrChromSCR - vChromBorderSCR.begin();
            break;
        }
    }

    //(2)Get discordant reads
    for(vector<St_ChromDiscordantReads>::iterator itrChromDR = vChromDR.begin();
        itrChromDR != vChromDR.end(); itrChromDR++)
    {
        if(itrChromDR->strChrom == strSTDSvDelChromName)
        {
            iDRIndex = itrChromDR - vChromDR.begin();
            break;
        }
    }

    if( iSCRIndex == -1 ||
        iDRIndex == -1)
    {
        return false;
    }
    else
        return true;
}

void SaveToFileEigenDelResult()
{
    string strCmd = "cat ./Output/Result_* > ./Result_Del_Sum.txt";
    system(strCmd.c_str());
}

void PrintEigenDelResult(vector<St_FindSvDEL>& vFindSvDEL)
{   
    cout << " ====== Final Sum ========" << endl;
    int iPosSum = 0;
    int iNegSum = 0;


    cout << "--------------" << endl;
    cout << "Total Positive: " << to_string(iPosSum) << endl;
    cout << "Total Negative: " << to_string(iNegSum) << endl;
    cout << "Total Sum     : " << to_string(iPosSum + iNegSum) << endl;

    cout << endl << "ALL SET!!!" << endl;

    St_ComparisonResult stCR;
    stCR.enSoftware = swEigenDel;
    //1: Get Result
    for(vector<St_FindSvDEL>::iterator itr = vFindSvDEL.begin();
        itr != vFindSvDEL.end(); itr++)
    {
        stCR.vSvResult.push_back(St_SvResult(itr->stFinalResult.strChrom,
                                             itr->stFinalResult.iStdNum,
                                             itr->stFinalResult.iPositiveNum + itr->stFinalResult.iNegativeNum,
                                             itr->stFinalResult.iPositiveNum,
                                             itr->stFinalResult.iNegativeNum));
    }
    //2: Print Result
    ClsComparison* pCompare = new ClsComparison();
    pCompare->PrintResult(stCR);    
    delete pCompare;
    pCompare = NULL;
}

void* FindSvDelForSingleChrom(void* pValue)
{
    St_FindSvDEL* pFindSvDEL = (St_FindSvDEL*)pValue; // Then Do --> The remaining work later -->

    if(pFindSvDEL->strChromName == "")
        return NULL;

    if(pFindSvDEL->pChromDR == NULL || pFindSvDEL->pChromSCR == NULL)
        return NULL;

//    if(pFindSvDEL->strBamFilePath == "" ||
//       ::access(pFindSvDEL->strBamFilePath.c_str(), 0) != 0)
//    {
//        delete pFindSvDEL;
//        pFindSvDEL = NULL;
//        return NULL;
//    }

    //Create a new file to record all temporary results -->
    ofstream ofs;
    string strFileName = "./Output/SvDEL_Chrom_" + pFindSvDEL->strChromName + ".txt";
    ofs.open(strFileName.c_str());
    //<--

    //Calc Insert Size
    ClsParseBam* pClsParseBam = new ClsParseBam();
    pClsParseBam->SetOfstram(ofs);
    pClsParseBam->SetChromName(pFindSvDEL->strChromName);

#ifdef GLOBALWHOLEGENO
    pClsParseBam->SetMeanInsertSize(pFindSvDEL->fAvgInsertSize);
    pClsParseBam->SetStdDevInsertSize(pFindSvDEL->fSTDInsertSize);
    pClsParseBam->SetAvgDepth(pFindSvDEL->fAvgDepth);
#else
#ifdef RUN_IN_CLUSTER
    pClsParseBam->CalcInsertSize(pFindSvDEL->strBamFilePath);
#else
    pClsParseBam->SetMeanInsertSize(410.026265);//(379);
    pClsParseBam->SetStdDevInsertSize(91.85553);//(38);
#endif

#ifdef RUN_IN_CLUSTER
    pClsParseBam->CalcAvgDepth(pFindSvDEL->strBamFilePath);
#else
    pClsParseBam->SetAvgDepth(7.0);
#endif
#endif

    ofs << "Mean Insert Size         : " << FloatToStr(pClsParseBam->GetMeanInsertSize()) << endl;
    ofs << "Std Deviation Insert Size: " << FloatToStr(pClsParseBam->GetStdDevInsertSize()) << endl;
    ofs << "Avg Depth                : " << FloatToStr(pClsParseBam->GetAvgDepth()) << endl;

//    //Prepare Target Reads: discordant reads and clipped reads
//    vector<St_ChromBorderSCR> vChromBorderSCR;
//    vector<St_ChromDiscordantReads> vChromDR;
//    //vector<St_BorderSCR> vBorderSCR;
//    //vector<St_DiscordantReads> vDiscdRreads;
//    PrepareTargetReads(pClsParseBam, pFindSvDEL->strBamFilePath, *(pFindSvDEL->pVFasta), vChromBorderSCR, vChromDR);
//    ofs << "Finish PrepareTargetReads" << endl;

//    //Get Result for each different chromosome respectively ->
//    St_ChromDiscordantReads* pChromDR = NULL;
//    St_ChromBorderSCR* pChromSCR = NULL;
//    St_ChromSV* pChromSvDel = NULL;

//    for(vector<St_ChromDiscordantReads>::iterator itrChromDR = vChromDR.begin();
//        itrChromDR != vChromDR.end(); itrChromDR++)
//    {
//        //(1)Get Discordant Reads
//        string strChrom = itrChromDR->strChrom;
//        pChromDR = &(*itrChromDR);

//        //(2)Get Border Softclip Reads
//        for(vector<St_ChromBorderSCR>::iterator itrChromSCR = vChromBorderSCR.begin();
//            itrChromSCR != vChromBorderSCR.end(); itrChromSCR++)
//        {
//            if(itrChromSCR->strChrom == strChrom)
//            {
//                pChromSCR = &(*itrChromSCR);
//                break;
//            }
//        }

//        //(3)Get related standard SV Deletion -->
//        if(pFindSvDEL->pChromSV->strChrom == strChrom)
//        {
//            pChromSvDel = pFindSvDEL->pChromSV;
//        }

//        if( pChromDR == NULL ||
//            pChromSCR == NULL ||
//            pChromSvDel == NULL)
//        {}
//        else
//        {
//            pFindSvDEL->stFinalResult.strChrom = strChrom;
//            ofs << "Do FindSvDelByChromosome For -- " << strChrom << endl;
//            FindSvDelByChromosome(pClsParseBam, pFindSvDEL->strBamFilePath,
//                                  pFindSvDEL->strPythonCodePath,
//                                  strChrom, pFindSvDEL->strSampleName,
//                                  pChromSCR->vSCR, pChromDR->vDR, pChromSvDel->vSv,
//                                  pFindSvDEL->stFinalResult, ofs);
//        }
//        pChromDR = NULL;
//        pChromSCR = NULL;
//        pChromSvDel = NULL;
//    }
//    //<--

    string strChrom = pFindSvDEL->pSTDChromSV->strChrom;
    pFindSvDEL->stFinalResult.strChrom = strChrom;
    ofs << "Do FindSvDelByChromosome For -- " << strChrom << endl;
    FindSvDelByChromosome(pClsParseBam, pFindSvDEL->strBamFilePath,
                          pFindSvDEL->strPythonCodePath,
                          strChrom, pFindSvDEL->strSampleName,
                          pFindSvDEL->pChromSCR->vSCR, pFindSvDEL->pChromDR->vDR, pFindSvDEL->pSTDChromSV->vSv,
                          pFindSvDEL->stFinalResult, ofs);

    ofs.close();
    delete pClsParseBam;
    pClsParseBam = NULL;
    pthread_exit(NULL);
}

void FindSvDelByChromosome(ClsParseBam* pClsParseBam, string strBamFile,
                           string strPythonCode,
                           string strChrom, string strSampleName,
                           vector<St_BorderSCR>& vBorderSCR, vector<St_DiscordantReads>& vDiscdRreads,
                           vector<St_SV>& vSvDEL, St_FinalResult& stFinalResult, ofstream& ofs)
{
    //Prepare Potential Deletion Group
    ClsSvDelDetect* pSvDelDetect = new ClsSvDelDetect();    
    pSvDelDetect->SetMeanInsertSize(pClsParseBam->GetMeanInsertSize());//(379);
    pSvDelDetect->SetStdDevInsertSize(pClsParseBam->GetStdDevInsertSize());//(38);
    pSvDelDetect->SetAvgDepth(pClsParseBam->GetAvgDepth()); //(4.82903)
    pSvDelDetect->SetOfstream(ofs);
    pSvDelDetect->SetChromName(strChrom);

    ClsDebug* pDebug = new ClsDebug();
    pDebug->SetStrBamFile(strBamFile);
    pDebug->SetMeanInserSize(379);
    pDebug->SetOfstream(ofs);
    pDebug->SetChromName(strChrom);

    vector<St_DRGroup> vDRGroup;
    vector<St_SCRGroup> vSCRGroup;
    PrepareGroup(pSvDelDetect, pDebug, vBorderSCR, vDiscdRreads, vSvDEL,
                 vDRGroup, vSCRGroup, strBamFile, strChrom);

    ClsLearning* pLearning = new ClsLearning();
    pLearning->SetAvgDepth(pClsParseBam->GetAvgDepth());
    pLearning->SetOfstream(ofs);
    pLearning->SetChromName(strChrom);
    pLearning->SetSampleName(strSampleName);

    //Generate the data matrix for both positive and negative training data <-- Go!!!
    ///Positive: comes from the real SV
    ///Negative: comes from the groups i get but do not belongs to SV
    // --> Go!!!
    pLearning->PrintFeatures(vDRGroup, vSCRGroup, vSvDEL);
    pLearning->PrintClusterElements(vDRGroup, vSvDEL);
    pLearning->ClusterByPython(vDRGroup, vSvDEL, strPythonCode, stFinalResult);
    //-->

    delete pSvDelDetect;
    pSvDelDetect = NULL;

    delete pDebug;
    pDebug = NULL;

    delete pLearning;
    pLearning = NULL;
}

void PrepareTargetReads(ClsParseBam* pClsParseBam, string strBamFile, vector<St_Fasta>& vFasta,
                        vector<St_ChromBorderSCR>& vChromBorderSCR,
                        vector<St_ChromDiscordantReads>& vChromDR)
                        //vector<St_BorderSCR>& vBorderSCR, vector<St_DiscordantReads>& vDiscdRreads)
{
//    //1: Get Border Softclip reads
//    cout << "Prepare GetBorderSCR" << endl;
//    pClsParseBam->GetBorderSCR(strBamFile, vFasta, vChromBorderSCR);

//    //2: Get Discordant Reads
//    cout << "Prepare GetDiscordantReads" << endl;
//    pClsParseBam->GetDiscordantReads(strBamFile, vFasta, vChromDR);
// *******************************************************************
    pClsParseBam->GetSCRAndDR(strBamFile, vFasta, vChromBorderSCR, vChromDR);
}

void PrepareGroup(ClsSvDelDetect* pSvDelDetect, ClsDebug* pDebug,
                  vector<St_BorderSCR>& vBorderSCR, vector<St_DiscordantReads>& vDiscdRreads,
                  vector<St_SV>& vSvDEL, vector<St_DRGroup>& vDRGroup, vector<St_SCRGroup>& vSCRGroup,
                  string strBamFile, string strChrom) //Go tomorrow
{
    //Step 1: Prepare group by discordant date
    //vector<St_DRGroup> vDRGroup;
    pSvDelDetect->ClusterDiscrodantReadsDelly(vDRGroup, vDiscdRreads, vBorderSCR);

    //Step 2: Prepare group by Border soft clip reads --> Go Later
    //vector<St_SCRGroup> vSCRGroup;
    //pSvDelDetect->ClusterSCReads(vSCRGroup, vBorderSCR);

    //Use depth to make the additional filter -->  this is our next step --> Go!! --> Finish it before back home --> Go!!!

    pSvDelDetect->DepthFilter(vDRGroup, vSCRGroup, strBamFile, strChrom);

    //Step 4: Collect the additiona types of read for each group (for example: one mapped and one unmapped)    
    pSvDelDetect->GetSpeciReads(vDRGroup, strBamFile);

    //Collect Features    
    pSvDelDetect->CollectFeatures(vDRGroup, strBamFile, strChrom);

    //Step 2.1: Only Keep the Groups realted to Standard SV    
    pDebug->GetStdSvDelRelatedGroup(vDRGroup, vSCRGroup, vSvDEL, strChrom);
    pDebug->UpdateGroupBoundary(); //Only for display !!!    
}

void DrawImage()
{}

void Comparision(string strConfigFile)
{
    //Read Config File
    ClsConfig* pConfig = new ClsConfig();
    St_Config stConfig;
    pConfig->ReadConfig(stConfig, strConfigFile.c_str());
    delete pConfig;
    pConfig = NULL;

    //Read VCF File
    ClsVcf1000Genome* pClsVcf = new ClsVcf1000Genome();
    cout << "ParseVcf --> " << endl;
    pClsVcf->ParseVcf(stConfig.strVcf, stConfig.strSampleName);

    //vector<St_SV> vSvDEL; // only pick out deletion
    vector<St_ChromSV> vChromSvDEL;
    cout << "GetDeletion --> " << endl;
    pClsVcf->GetDeletion(vChromSvDEL, stConfig.strSampleName, stConfig.strChrom);
    cout << "vChromSvDEL size: " << to_string(vChromSvDEL.size()) << endl;
    cout << "==== Detail ====" << endl;
    int iTotal = 0;
    for(vector<St_ChromSV>::iterator itr = vChromSvDEL.begin(); itr != vChromSvDEL.end(); itr++)
    {
        cout << itr->strChrom << ": " << to_string(itr->vSv.size()) << endl;
        iTotal += itr->vSv.size();
    }
    cout << "----" << endl << "Total: " << to_string(iTotal) << endl;
    cout << "================" << endl;
    delete pClsVcf;
    pClsVcf = NULL;

    ClsComparison* pComparison = new ClsComparison();
    //pComparison->CompareStdSvWithSprites(stConfig.strSpritesPath, vChromSvDEL);
    pComparison->CompareStdSVWithPindel(stConfig.strPindelPath, vChromSvDEL);
    pComparison->CompareStdSVWithDelly(stConfig.strDellyPath, vChromSvDEL);
    pComparison->CompareStdSVWithLumpy(stConfig.strLumpyPath, vChromSvDEL);//vSvDEL);
    pComparison->CompareStdSVWithCNVnator(stConfig.strCNVnatorPath, vChromSvDEL);
    pComparison->CompareStdSVWithGASVpro(stConfig.strGASVproPath, vChromSvDEL);

    delete pComparison;
    pComparison = NULL;
}
