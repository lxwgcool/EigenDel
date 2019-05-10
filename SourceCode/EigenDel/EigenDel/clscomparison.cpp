#include "clscomparison.h"
#include <fstream>
#include "../../../ShareLibrary/clsbasealgorithm.h"

char* ArrySoftware[swMax] = {"Pindel", "Sprites", "GASVpro", "CNVnator", "Delly", "Lumpy", "EigenDel"};

ClsComparison::ClsComparison()
{
}

void ClsComparison::ParsePindel(string strFilePath, vector<St_BreakPoint>& vBP)
{
    fstream ifsPindel;
    ifsPindel.open(strFilePath.c_str());
    string strLine = "";
    vBP.clear();
    St_BreakPoint stBP;
    while(!ifsPindel.eof())
    {
        getline(ifsPindel, strLine);

        if(strLine.find("##########") == string::npos) // do not find
            continue;

        //If find it
        getline(ifsPindel, strLine);
        int iStart = 0;
        int iLen = 0;
        int iEnd = 0;

        //For Chrom ID
        iStart = strLine.find("ChrID");
        iStart += 5;
        iEnd = strLine.find('\t', iStart);
        iLen = iEnd - iStart;
        stBP.strChromID = strLine.substr(iStart, iLen);
        trim(stBP.strChromID);

        //For BP Start position
        iStart = strLine.find("BP", iEnd);
        iStart = iStart + 2;
        iEnd = strLine.find('\t', iStart);
        iLen = iEnd - iStart;
        stBP.iStart = atoi(strLine.substr(iStart, iLen).c_str());

        //For ending position
        iStart = iEnd + 1;
        iEnd = strLine.find('\t', iStart);
        iLen = iEnd - iStart;
        stBP.iEnd = atoi(strLine.substr(iStart, iLen).c_str());

        if(stBP.iStart > stBP.iEnd)
        {
            int iTmp = stBP.iStart;
            stBP.iStart = stBP.iEnd;
            stBP.iEnd = iTmp;
        }

        vBP.push_back(stBP);
    }
}

void ClsComparison::CompareStdSVWithPindel(string strFilePath, vector<St_ChromSV>& vChromSvDEL)
{
    //1: Parse Pindel Result
    vector<St_BreakPoint> vBP;
    ParsePindel(strFilePath, vBP);

    //2: Converet Pindel Result to vSoftwareSvDEL
    vector<St_ChromSV> vSoftwareSvDEL;
    St_ChromSV stChromSv;
    for(vector<St_BreakPoint>::iterator itr = vBP.begin(); itr != vBP.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_ChromSV>::iterator subItr = vSoftwareSvDEL.begin();
            subItr != vSoftwareSvDEL.end(); subItr++)
        {
            if(subItr->strChrom == itr->strChromID)
            {
                subItr->vSv.push_back(St_SV(itr->strChromID, itr->iStart, itr->iEnd));
                bFind = true;
                break;
            }
        }
        if(!bFind) // new one
        {
            stChromSv.Clear();
            stChromSv.strChrom = itr->strChromID;
            stChromSv.vSv.push_back(St_SV(itr->strChromID, itr->iStart, itr->iEnd));
            vSoftwareSvDEL.push_back(stChromSv);
        }
    }

    //3: Get info by compare with standard database
    St_ComparisonResult stCR;
    stCR.enSoftware = swPindel;
    GetByCompareStd(vChromSvDEL, vSoftwareSvDEL, stCR);

    //4: Print result
    PrintResult(stCR);

//    const int MAXDIFF = 340;

//    //Compare to check how many covered
//    int iMiss = 0;
//    for(vector<St_SV>::iterator itr = vSvDEL.begin(); itr != vSvDEL.end(); itr++)
//    {
//        bool bFind =  false;
//        for(vector<St_BreakPoint>::iterator subItr = vBP.begin(); subItr != vBP.end(); subItr++)
//        {
//            if(abs(itr->iPos - subItr->iStart) <= MAXDIFF ||
//               abs(itr->iEnd - subItr->iEnd) <= MAXDIFF)
//            {
//                bFind = true;
//            }
//        }
//        if(!bFind)
//        {
//            cout << "Miss SV: " << " --- " << "Len: " << to_string(itr->iEnd - itr->iPos) << " --- "
//                 << "(" << to_string(itr->iPos) << ", "
//                 << to_string(itr->iEnd) << ")" << endl;
//            iMiss++;
//        }
//    }
//    cout << endl << "-----------" << endl
//      << "Miss: " << to_string(iMiss) << endl << "-----------" << endl;
//    cout << "STD SV(DEL) Total   : " << to_string(vSvDEL.size()) << endl;
//    cout << "Pindel SV(DEL) Tatal: " << to_string(vBP.size()) << endl;
}

void ClsComparison::ParseSprites(string strFilePath, vector<St_BreakPoint>& vBP)
{
    fstream ifsSprites;
    ifsSprites.open(strFilePath.c_str());
    string strLine = "";
    vBP.clear();
    St_BreakPoint stBP;
    while(!ifsSprites.eof())
    {
        getline(ifsSprites, strLine);

        //Parse it:
        int iStart = 0;
        int iEnd = 0;
        int iLen = 0;

        //For chrom name: --> column 0
        iEnd = strLine.find('\t', iStart);
        iLen = iEnd - iStart;
        stBP.strChromID = strLine.substr(iStart, iLen);
        iStart = iEnd + 1;

        //For start position of DEL --> column 2
        iEnd = strLine.find('\t', iStart);
        iStart = iEnd + 1;
        iEnd = strLine.find('\t', iStart);
        iLen = iEnd - iStart;
        stBP.iStart = atoi(strLine.substr(iStart, iLen).c_str());
        iStart = iEnd + 1;

        //For end position of DEL --> column 4
        iEnd = strLine.find('\t', iStart);
        iStart = iEnd + 1;
        iEnd = strLine.find('\t', iStart);
        iLen = iEnd - iStart;
        stBP.iEnd = atoi(strLine.substr(iStart, iLen).c_str());

        if(stBP.iStart > stBP.iEnd)
        {
            int iTmp = stBP.iStart;
            stBP.iStart = stBP.iEnd;
            stBP.iEnd = iTmp;
        }
        vBP.push_back(stBP);
    }
}

void ClsComparison::ParseGASVpro(string strFilePath, vector<St_BreakPoint>& vBP)
{
    fstream ifsGASVpro;
    ifsGASVpro.open(strFilePath.c_str());
    string strLine = "";
    vBP.clear();
    St_BreakPoint stBP;
    while(!ifsGASVpro.eof())
    {
        getline(ifsGASVpro, strLine);

        //skip the first line
        if(strLine[0] == '#')
            continue;

        //Check if it belongs to deletion:
        int iTypePos = strLine.rfind('\t') + 1;
        string strType = strLine.substr(iTypePos, strLine.length()-iTypePos);
        if(strType != "D") //we only keep the type "D"
            continue;

        //Parse it:
        int iStart = 0;
        int iEnd = 0;
        int iLen = 0;

        //For chrom name: --> column 1 (column starts from 0)
        iEnd = strLine.find('\t', iStart);
        iStart = iEnd + 1;
        iEnd = strLine.find('\t', iStart);
        iLen = iEnd - iStart;
        stBP.strChromID = strLine.substr(iStart, iLen);
        //Convert 23 to X and 24 to Y
        if(stBP.strChromID == "23")
            stBP.strChromID = "X";
        else if(stBP.strChromID == "24")
            stBP.strChromID = "Y";
        iStart = iEnd + 1;

        //For start position of DEL --> column 2
        iEnd = strLine.find('\t', iStart);
        iLen = iEnd - iStart;
        string strTmp = strLine.substr(iStart, iLen);
        int iSplit= strTmp.find(',');
        stBP.iStart = atoi(strTmp.substr(0, iSplit).c_str());
        stBP.iEnd = atoi(strTmp.substr(iSplit+1, strTmp.length()-(iSplit+1)).c_str());

        if(stBP.iStart > stBP.iEnd)
        {
            int iTmp = stBP.iStart;
            stBP.iStart = stBP.iEnd;
            stBP.iEnd = iTmp;
        }

        vBP.push_back(stBP);
    }
}

void ClsComparison::ParseCNVnator(string strFilePath, vector<St_BreakPoint>& vBP)
{
    fstream ifsCNVnator;
    ifsCNVnator.open(strFilePath.c_str());
    string strLine = "";
    vBP.clear();
    St_BreakPoint stBP;
    while(!ifsCNVnator.eof())
    {
        getline(ifsCNVnator, strLine);

        //Parse it:
        int iStart = 0;
        int iEnd = 0;
        int iLen = 0;

        //Check if it belongs to deletion:
        iEnd = strLine.find('\t');
        iLen = iEnd - iStart;
        string strType = strLine.substr(iStart, iLen);
        if(strType != "deletion") //we only keep the type "deletion"
            continue;
        else
            iStart = iEnd + 1;

        //For chrom name: --> column 1 (column starts from 0)
        iEnd = strLine.find('\t', iStart);
        iLen = iEnd - iStart;
        string strTmp = strLine.substr(iStart, iLen);
        int iTmpSplit = strTmp.find(':');
        int iTmpStart = 0;
        int iTmpLen = iTmpSplit - iTmpStart;
        stBP.strChromID = strTmp.substr(iTmpStart, iTmpLen);

        iTmpStart = iTmpSplit + 1;
        iTmpSplit = strTmp.find('-');
        iTmpLen = iTmpSplit - iTmpStart;
        stBP.iStart = atoi(strTmp.substr(iTmpStart, iTmpLen).c_str());

        iTmpStart = iTmpSplit + 1;
        iTmpLen = strTmp.length() - iTmpStart;
        stBP.iEnd = atoi(strTmp.substr(iTmpStart, iTmpLen).c_str());

        if(stBP.iStart > stBP.iEnd)
        {
            int iTmp = stBP.iStart;
            stBP.iStart = stBP.iEnd;
            stBP.iEnd = iTmp;
        }

        vBP.push_back(stBP);
    }
}

void ClsComparison::CompareStdSvWithSprites(string strFilePath, vector<St_ChromSV>& vChromSvDEL)
{
    //1: Parse Sprites Result
    vector<St_BreakPoint> vBP;
    ParseSprites(strFilePath, vBP);

    //2: Converet Pindel Result to vSoftwareSvDEL
    vector<St_ChromSV> vSoftwareSvDEL;
    St_ChromSV stChromSv;
    for(vector<St_BreakPoint>::iterator itr = vBP.begin(); itr != vBP.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_ChromSV>::iterator subItr = vSoftwareSvDEL.begin();
            subItr != vSoftwareSvDEL.end(); subItr++)
        {
            if(subItr->strChrom == itr->strChromID)
            {
                subItr->vSv.push_back(St_SV(itr->strChromID, itr->iStart, itr->iEnd));
                bFind = true;
                break;
            }
        }
        if(!bFind) // new one
        {
            stChromSv.Clear();
            stChromSv.strChrom = itr->strChromID;
            stChromSv.vSv.push_back(St_SV(itr->strChromID, itr->iStart, itr->iEnd));
            vSoftwareSvDEL.push_back(stChromSv);
        }
    }

    //3: Get info by compare with standard database
    St_ComparisonResult stCR;
    stCR.enSoftware = swSprites;
    GetByCompareStd(vChromSvDEL, vSoftwareSvDEL, stCR);

    //4: Print result
    PrintResult(stCR);
}

void ClsComparison::CompareStdSVWithGASVpro(string strFilePath, vector<St_ChromSV>& vChromSvDEL)
{
    //1: Parse GASVpro Result
    vector<St_BreakPoint> vBP;
    ParseGASVpro(strFilePath, vBP);

    //2: Converet GASVpro Result to vSoftwareSvDEL
    vector<St_ChromSV> vSoftwareSvDEL;
    St_ChromSV stChromSv;
    for(vector<St_BreakPoint>::iterator itr = vBP.begin(); itr != vBP.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_ChromSV>::iterator subItr = vSoftwareSvDEL.begin();
            subItr != vSoftwareSvDEL.end(); subItr++)
        {
            if(subItr->strChrom == itr->strChromID)
            {
                subItr->vSv.push_back(St_SV(itr->strChromID, itr->iStart, itr->iEnd));
                bFind = true;
                break;
            }
        }
        if(!bFind) // new one
        {
            stChromSv.Clear();
            stChromSv.strChrom = itr->strChromID;
            stChromSv.vSv.push_back(St_SV(itr->strChromID, itr->iStart, itr->iEnd));
            vSoftwareSvDEL.push_back(stChromSv);
        }
    }

    //3: Get info by compare with standard database
    St_ComparisonResult stCR;
    stCR.enSoftware = swSvGASVpro;
    GetByCompareStd(vChromSvDEL, vSoftwareSvDEL, stCR);

    //4: Print result
    PrintResult(stCR);
}

void ClsComparison::CompareStdSVWithCNVnator(string strFilePath, vector<St_ChromSV>& vChromSvDEL)
{
    //1: Parse CNVnator Result
    vector<St_BreakPoint> vBP;
    ParseCNVnator(strFilePath, vBP);

    //2: Converet GASVpro Result to vSoftwareSvDEL
    vector<St_ChromSV> vSoftwareSvDEL;
    St_ChromSV stChromSv;
    for(vector<St_BreakPoint>::iterator itr = vBP.begin(); itr != vBP.end(); itr++)
    {
        bool bFind = false;
        for(vector<St_ChromSV>::iterator subItr = vSoftwareSvDEL.begin();
            subItr != vSoftwareSvDEL.end(); subItr++)
        {
            if(subItr->strChrom == itr->strChromID)
            {
                subItr->vSv.push_back(St_SV(itr->strChromID, itr->iStart, itr->iEnd));
                bFind = true;
                break;
            }
        }
        if(!bFind) // new one
        {
            stChromSv.Clear();
            stChromSv.strChrom = itr->strChromID;
            stChromSv.vSv.push_back(St_SV(itr->strChromID, itr->iStart, itr->iEnd));
            vSoftwareSvDEL.push_back(stChromSv);
        }
    }

    //3: Get info by compare with standard database
    St_ComparisonResult stCR;
    stCR.enSoftware = swCNVnator;
    GetByCompareStd(vChromSvDEL, vSoftwareSvDEL, stCR);

    //4: Print result
    PrintResult(stCR);
}

void ClsComparison::CompareStdSVWithLumpy(string strFilePath, vector<St_ChromSV>& vChromSvDEL)
{
//    ClsVcf1000Genome* pClsVcf = new ClsVcf1000Genome();
//    pClsVcf->ParseVcf(strFilePath);
//    vector<St_ChromSV> vChromSvLumpyDEL; // only pick out deletion
//    pClsVcf->GetDeletion(vChromSvLumpyDEL);
//    cout << "vSvDEL size: " << to_string(vChromSvLumpyDEL.size()) << endl;
//    delete pClsVcf;
//    pClsVcf = NULL;

//    //Compare
//    //Compare to check how many covered
//    const int MAXDIFF = 340; //insert length
//    int iMiss = 0;
//    int iHit = 0;
////    for(vector<St_SV>::iterator itr = vSvDEL.begin(); itr != vSvDEL.end(); itr++)
////    {
////       bool bFind =  false;
////       for(vector<St_SV>::iterator subItr = vSvLumpyDEL.begin(); subItr !=vSvLumpyDEL.end(); subItr++)
////       {
////           if(abs(itr->iPos - subItr->iPos) <= MAXDIFF ||
////              abs(itr->iEnd - subItr->iEnd) <= MAXDIFF)
////           {
////               bFind = true;
////           }
////       }
////       if(!bFind)
////       {
////           cout << "Miss SV: " << " --- " << "Len: " << to_string(itr->iEnd - itr->iPos) << " --- "
////                << "(" << to_string(itr->iPos) << ", "
////                << to_string(itr->iEnd) << ")" << endl;
////           iMiss++;
////       }
////       else
////           iHit++;
////    }

//    for(vector<St_ChromSV>::iterator itrChromSVLumpy = vChromSvLumpyDEL.begin();
//        itrChromSVLumpy !=vChromSvLumpyDEL.end(); itrChromSVLumpy++)
//    {
//        St_ChromSV* pChromSV = NULL;

//        for(vector<St_ChromSV>::iterator itrChromSV = vChromSvDEL.begin();
//            itrChromSV != vChromSvDEL.end(); itrChromSV++)
//        {
//            if(itrChromSVLumpy->strChrom == itrChromSV->strChrom)
//            {
//                pChromSV = &(*itrChromSV);
//                break;
//            }
//        }

//        if(pChromSV == NULL)
//            continue;

//        //Not NULL --> Now we get the correct chromosome --> Let do comparison
//        for(vector<St_SV>::iterator itrLumpySv = itrChromSVLumpy->vSv.begin();
//            itrLumpySv != itrChromSVLumpy->vSv.end(); itrLumpySv++)
//        {
//            bool bFind = false;
//            for(vector<St_SV>::iterator subItr = pChromSV->vSv.begin(); subItr != pChromSV->vSv.end(); subItr++)
//            {
//                if(abs(itrLumpySv->iPos - subItr->iPos) <= MAXDIFF ||
//                   abs(itrLumpySv->iEnd - subItr->iEnd) <= MAXDIFF)
//               {
//                   bFind = true;
//                   break;
//               }
//            }

//            if(bFind)
//                iHit++;
//            else
//                iMiss++;
//        }
//    }

//    int iTotalLumpySv = 0;
//    int iTotalSTDSv = 0;

//    for(vector<St_ChromSV>::iterator itrChromSVLumpy = vChromSvLumpyDEL.begin();
//        itrChromSVLumpy !=vChromSvLumpyDEL.end(); itrChromSVLumpy++)
//    {
//        iTotalLumpySv += itrChromSVLumpy->vSv.size();
//    }

//    for(vector<St_ChromSV>::iterator itrChromSV = vChromSvDEL.begin();
//        itrChromSV != vChromSvDEL.end(); itrChromSV++)
//    {
//        iTotalSTDSv += itrChromSV->vSv.size();
//    }

//    cout << endl << "-----------" << endl
//         << "Miss: " << to_string(iMiss) << endl << "-----------" << endl;
//    cout << "STD SV(DEL)   Total: " << to_string(iTotalSTDSv) << endl << endl;
//    cout << "Positive(DEL)      : " << to_string(iHit) << endl;
//    cout << "Negative(DEL)      : " << to_string(iMiss) << endl;
//    cout << "Lumpy SV(DEL) Tatal: " << to_string(iTotalLumpySv) << endl;

    St_ComparisonResult stCR;
    stCR.enSoftware = swLumpy;
    //1: Get Result
    GetResultFromVcf(strFilePath, vChromSvDEL, stCR);
    //2: Print Result
    PrintResult(stCR);
}

void ClsComparison::CompareStdSVWithDelly(string strFilePath, vector<St_ChromSV>& vChromSvDEL)
{
    St_ComparisonResult stCR;
    stCR.enSoftware = swDelly;
    //1: Get Result
    GetResultFromVcf(strFilePath, vChromSvDEL, stCR);
    //2: Print Result
    PrintResult(stCR);
}

void ClsComparison::GetResultFromVcf(string strFilePath,
                                     vector<St_ChromSV>& vStdSvDEL,
                                     St_ComparisonResult& stCR)
{
    //1: get SV DEL record in VCF
    ClsVcf1000Genome* pClsVcf = new ClsVcf1000Genome();
    pClsVcf->ParseVcf(strFilePath);
    vector<St_ChromSV> vSoftwareSvDEL; // only pick out deletion
    pClsVcf->GetDeletion(vSoftwareSvDEL);
    //cout << "vSvDEL size: " << to_string(vChromSvLumpyDEL.size()) << endl;
    delete pClsVcf;
    pClsVcf = NULL;

    //Compare with the standard SV and record the result of comparison
    GetByCompareStd(vStdSvDEL, vSoftwareSvDEL, stCR);
}

void ClsComparison::GetByCompareStd(vector<St_ChromSV>& vStdSvDEL,
                                    vector<St_ChromSV> vSoftwareSvDEL,
                                    St_ComparisonResult& stCR)
{
    const int MAXDIFF = 340; //insert length --> boundary diff tolerance
    St_SvResult stResult;

    for(vector<St_ChromSV>::iterator itrSoftwareSv = vSoftwareSvDEL.begin();
        itrSoftwareSv != vSoftwareSvDEL.end(); itrSoftwareSv++)
    {
        //(1): Record the detected SV message of current chromosome
        stResult.strChrom = itrSoftwareSv->strChrom;
        stResult.iSvFind = itrSoftwareSv->vSv.size();

        //(2): Find if we can find this chromosome in Standard SV database
        St_ChromSV* pChromSV = NULL;
        for(vector<St_ChromSV>::iterator itrChromSV = vStdSvDEL.begin();
            itrChromSV != vStdSvDEL.end(); itrChromSV++)
        {
            if(stResult.strChrom == itrChromSV->strChrom)
            {
                pChromSV = &(*itrChromSV);
                break;
            }
        }

        if(pChromSV == NULL)
        {
            //stCR.vSvResult.push_back(stResult);
            stResult.Clear();
            continue;
        }
        else
            stResult.iStdSv = pChromSV->vSv.size();

        //Not NULL --> Now we get the correct chromosome --> Let do comparison
        for(vector<St_SV>::iterator itrSv = itrSoftwareSv->vSv.begin();
            itrSv != itrSoftwareSv->vSv.end(); itrSv++)
        {
            bool bFind = false;
            for(vector<St_SV>::iterator subItr = pChromSV->vSv.begin();
                subItr != pChromSV->vSv.end(); subItr++)
            {
                if(abs(itrSv->iPos - subItr->iPos) <= MAXDIFF ||
                   abs(itrSv->iEnd - subItr->iEnd) <= MAXDIFF)
               {
                   bFind = true;
                   break;
               }
            }

            if(bFind)
                stResult.iSvPositive++;
            else
                stResult.iSvNegative++;
        }

        stCR.vSvResult.push_back(stResult);
        stResult.Clear();
    }
}

void ClsComparison::PrintResult(St_ComparisonResult& stCR)
{
    cout << endl << "==========="
         << ArrySoftware[stCR.enSoftware]
         << "===========" << endl;

    //1: Out Put result for each
    cout << ">>>>>> Result of each chromosome" << endl;
    cout << "Chrom --- " << "Positive --- " << "Negative --- " << "STD --- "
         << "Accuracy --- " << "Sensitivity --- " << "F1-Score" << endl;
    for(vector<St_SvResult>::iterator itr = stCR.vSvResult.begin();
        itr != stCR.vSvResult.end(); itr++)
    {
        if(itr->iStdSv == 0)
           continue;

        cout << itr->strChrom << "    "
             << to_string(itr->iSvPositive) << "   "
             << to_string(itr->iSvNegative) << "   "
             << to_string(itr->iStdSv) << " --- "
             << FloatToStr(itr->CalcAccuracy(), 4) << "   "
             << FloatToStr(itr->CalcSensitivity(), 4) << "   "
             << FloatToStr(itr->CalcF1Score(), 4) << endl;
    }
    //2: Output Statistic result
    cout << ">>>>>> Result of each chromosome" << endl;
    cout << "STD SV(DEL)   Total: " << to_string(stCR.GetTotalStdSv()) << endl;
    cout << ArrySoftware[stCR.enSoftware] << " Tatal: " << to_string(stCR.GetSvFind()) << endl;
    cout << "Positive(DEL)      : " << to_string(stCR.GetSvPositive()) << endl;
    cout << "Negative(DEL)      : " << to_string(stCR.GetSvNegative()) << endl;
    cout << "Accuracy           : " << FloatToStr(stCR.CalcAccuracy(), 4) << endl;
    cout << "Sensitivity        : " << FloatToStr(stCR.CalcSensitivity(), 4) << endl;
    cout << "F1-Score           : " << FloatToStr(stCR.CalcF1Score(), 4) << endl;
    cout << "--------------------" << endl;
}
