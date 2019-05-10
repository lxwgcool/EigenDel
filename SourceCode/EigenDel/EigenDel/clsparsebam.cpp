#include "clsparsebam.h"
#include <algorithm>

ClsParseBam::ClsParseBam(): m_fMeanInsertSize(-1.),m_fStdDevInsertSize(-1.),
                            m_fAvgDepth(-1.), m_pOfs(NULL), m_strChromName("")
{}

ClsParseBam::~ClsParseBam()
{}

void ClsParseBam::ReadBamFile(string strBamFilePath, vector<BamAlignment>& vAl)
{
    //Go to check sorted bam file
    //1: parse bam file & output the expected reads
    BamReader* pBamReader = new BamReader();
    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");
    BamAlignment al;

    while(pBamReader->GetNextAlignment(al)) //Get each alignment result
    {
        if(!al.IsMapped() || al.QueryBases == "")
            continue;

        if(al.QueryBases.find('N') != string::npos ||
           al.QueryBases.find('n') != string::npos)
            continue;

        vAl.push_back(al);
    }

    delete pBamReader;
    pBamReader = NULL;
}

void ClsParseBam::CalcInsertSize(string strBamFilePath, string strPicard)
{
    //here we use picard to calculate the average insert size and the insert size standard deviation    
    string strMatricsFile = "./Output/insert_size_metrics_Chrom_" + m_strChromName + ".txt";
    string strHistogramFile = "./Output/insert_size_histogram_Chrom_" + m_strChromName + ".pdf";
    string strBreifISFile = "./Output/brief_info_insert_size_Chrom_" + m_strChromName + ".txt";
    //string strPicardOutput = "./picard_chrom_" + m_strChromName + ".output";

    string strCmd = (string)"java -XX:ParallelGCThreads=24 -jar " + strPicard + " CollectInsertSizeMetrics " +
                    "I=" + strBamFilePath + " " +
                    "O=" + strMatricsFile + " " +
                    "H=" + strHistogramFile + " " +
                    "M=0.5" + " " +
                    "VALIDATION_STRINGENCY=LENIENT"; // +
                    //" > " + strPicardOutput;
    cout << "Read: " << strBamFilePath << endl;
    system(strCmd.c_str());

    if(m_pOfs == NULL)
        cout << "Read Finished: " << strBamFilePath << endl;
    else
        (*m_pOfs) << "Read Finished: " << strBamFilePath << endl;


    //Get Breif insert size result
    strCmd = "head -8 " + strMatricsFile +
             " | tail -n 2 | awk '{ for (i=1; i<=NF; i++) RtoC[i]= (RtoC[i]? RtoC[i] FS $i: $i) } \
              END{ for (i in RtoC) print RtoC[i] }' | awk '{print $1 \"\t\" $2}' > " + strBreifISFile;

    system(strCmd.c_str());

    //For Mean Insert Size:
    strCmd = "awk '{print $2}' " + strBreifISFile + " | sed -n '6p; 7q'";
    string strResult = exec(strCmd.c_str());
    m_fMeanInsertSize = stof(strResult.c_str());

    //For Std Dev Insert Size:
    strCmd = "awk '{print $2}' " + strBreifISFile + " | sed -n '7p; 8q'";
    strResult = exec(strCmd.c_str());
    m_fStdDevInsertSize = stof(strResult.c_str());

    //Parse brief result
    //--->
    // m_iMeanInsertSize = ***;
    // m_iStdDevInsertSize = ***;
    //<---
}

void ClsParseBam::CalcAvgDepth(string strBamFilePath)
{
    string strCmd = (string)"samtools depth " + strBamFilePath +
                    " | awk '{sum+=$3} END {if(sum==\"\") {print 0} else {print sum/NR}}'";
    string strResult = exec(strCmd.c_str());
    m_fAvgDepth = stof(strResult.c_str());
    cout << "Average Depth: " << FloatToStr(m_fAvgDepth) << endl;
}

void ClsParseBam::SetMeanInsertSize(float fV1)
{
    this->m_fMeanInsertSize = fV1;
}

void ClsParseBam::SetStdDevInsertSize(float fV1)
{
    this->m_fStdDevInsertSize = fV1;
}

void ClsParseBam::SetOfstram(ofstream& ofs)
{
    this->m_pOfs = &ofs;
}

void ClsParseBam::SetChromName(string strChromName)
{
    this->m_strChromName = strChromName;
}

float ClsParseBam::GetMeanInsertSize()
{
    return this->m_fMeanInsertSize;
}

float ClsParseBam::GetStdDevInsertSize()
{
    return this->m_fStdDevInsertSize;
}

void ClsParseBam::SetAvgDepth(float fV1)
{
    this->m_fAvgDepth = fV1;
}

float ClsParseBam::GetAvgDepth()
{
    return this->m_fAvgDepth;
}

bool sort_scrClipPos_func(St_BorderSCR stScr1, St_BorderSCR stScr2)
{
    if(stScr1.iClipPos < stScr2.iClipPos)
        return true;
    else
        return false;
}

void ClsParseBam::GetBorderSCR(string strBamFilePath, vector<St_Fasta>& vFasta,
                               vector<St_ChromBorderSCR>& vChromBorderSCR)//vector<St_BorderSCR>& vBorderSCR)
{
    vector<St_ClipReads> vClipReads;
    St_ClipReads stClipReads;

    St_ChromBorderSCR stChromSCR;

    //---->Initial St_BorderSCR
    //vBorderSCR.clear();
    St_BorderSCR stBorderSCR;
    //<----

    BamReader* pBamReader = new BamReader();
    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");
    BamAlignment al;

//    int iNotStartEnd  = 0;

    while(pBamReader->GetNextAlignment(al)) //Get each alignment result
    {
        if(!al.IsMapped() || al.QueryBases == "")
            continue;

        if(al.IsDuplicate() || al.IsFailedQC() || !al.IsPrimaryAlignment()) //We should keep the case: al.IsReverseStrand()
            continue;

        if(al.QueryBases.find('N') != string::npos ||
           al.QueryBases.find('n') != string::npos)
            continue;                    

//        if(al.Position == 132558874)
//        {
//            cout << "I am the target!" << endl;
//            if(al.IsFirstMate())
//                cout << "I am the first mate!" << endl;

//            if(al.IsSecondMate())
//                cout << "I am the second mate!" << endl;

//            if(al.IsDuplicate())
//                cout << "I am duplicate!" << endl;
//        }

        bool bFindClip = false;
        int iOffSet = 0;

//        //----> Just For testing
//        if(al.IsSecondMate())
//        {
//            string strMapSeq = "";
//            string strClipSeq = "";
//            bool bClipped = false;
//            int iClipPos = 0;
//            for(std::vector<CigarOp>::iterator itr = al.CigarData.begin(); itr != al.CigarData.end(); itr++)
//            {
//                switch(itr->Type)
//                {
//                    case 'M': // alignment match (can be a sequence match or mismatch)
//                        strMapSeq += al.QueryBases.substr(iOffSet, itr->Length);
//                        iOffSet += itr->Length;
//                        break;
//                    case 'I': // insertion to the reference
//                        break;
//                    case 'D': // deletion from the reference
//                    case 'N':  // skipped region from the reference
//                        break;
//                    case 'S':  // soft clipping (clipped sequences present in SEQ)
//                        strClipSeq += al.QueryBases.substr(iOffSet, itr->Length);
//                        iClipPos = iOffSet;
//                        iOffSet += itr->Length;
//                        bClipped = true;
//                        break;
//                    case 'H':  // hard clipping (clipped sequences NOT present in SEQ)
//                        break;
//                    case 'P': // padding (silent deletion from padded reference)
//                    case '=': // sequence match
//                    case 'X': // sequence mismatch
//                        break;
//                }
//            }
//            if(bClipped)
//            {
//                cout << "Map Position : " << to_string(al.Position) << endl;
//                cout << "Clip Position: " << to_string(iClipPos) << endl;
//                cout << "Map Seq      : " << strMapSeq << endl;
//                cout << "Clip Seq     : " << strClipSeq << endl;
//                cout << "Query Seq    : " << al.QueryBases << endl;
//                cout << endl;
//            }
//        }
//        //<----
        iOffSet = 0;

        stClipReads.Clear();
        string strMapSeq = al.AlignedBases;
        string strClipSeq = "";
        int iQueryClipPos = -1;
        int iQueryMapPos = -1;
        for(std::vector<CigarOp>::iterator itr = al.CigarData.begin(); itr != al.CigarData.end(); itr++)
        {
            switch(itr->Type)
            {
                case 'M': // alignment match (can be a sequence match or mismatch)
                    //if(strMapSeq == "")
                    //    strMapSeq += al.QueryBases.substr(iOffSet, itr->Length);
                    if(iQueryMapPos == -1)
                        iQueryMapPos = iOffSet;
                    iOffSet += itr->Length;                    
                    break;
                case 'I': // insertion to the reference
                    break;
                case 'D': // deletion from the reference
                case 'N':  // skipped region from the reference
                    break;
                case 'S':  // soft clipping (clipped sequences present in SEQ)
                    if(strClipSeq == "")
                        strClipSeq += al.QueryBases.substr(iOffSet, itr->Length);
                    stClipReads.vClipSeq.push_back(al.QueryBases.substr(iOffSet, itr->Length));
                    stClipReads.vPos.push_back(iOffSet);                    

                    if(iQueryClipPos == -1)
                        iQueryClipPos = iOffSet;
                    iOffSet += itr->Length;
                    bFindClip = true;

                    break;
                case 'H':  // hard clipping (clipped sequences NOT present in SEQ)
                    break;
                case 'P': // padding (silent deletion from padded reference)
                case '=': // sequence match
                case 'X': // sequence mismatch
                    break;
            }
        }

        if(bFindClip)
        {
            if(strMapSeq.length() < 15 || strClipSeq.length() < 15)
                continue;

            int iRealInsertSize = al.InsertSize;
            if(al.IsSecondMate() &&
               abs(iRealInsertSize) > 100*2)
                iRealInsertSize = iRealInsertSize * -1;

//            if(iRealInsertSize == 0 && al.IsMateMapped()) //  this is belongs to the special case
//            {
//                iRealInsertSize = al.MatePosition - al.Position;
//                if(iRealInsertSize > 0)
//                    iRealInsertSize += al.QueryBases.length();
//                else
//                    iRealInsertSize -= al.QueryBases.length();
//            }

            stClipReads.strName = al.Name;
            stClipReads.strQuerySeq = al.QueryBases;
            stClipReads.iReadsMapPos = al.Position;
            stClipReads.iInsertSize = iRealInsertSize;
            stClipReads.iMateMapPos = al.MatePosition;
            stClipReads.bFirstMate = al.IsFirstMate();
            stClipReads.bSecondMate = al.IsSecondMate();
            stClipReads.bMateMapped = al.IsMateMapped();

            stClipReads.iQueryMapPos = iQueryMapPos;
            stClipReads.iQueryMapLength = strMapSeq.length();
            stClipReads.iQueryClipPos = iQueryClipPos;
            stClipReads.iQueryClipLength = strClipSeq.length();

            //Save clipped info
            vClipReads.push_back(stClipReads);

            //The number of clip which do not happed in the beginning or the ending part of the reads -->
            if(stClipReads.iQueryClipPos < 3 ||
               stClipReads.iQueryClipPos > (al.QueryBases.length() - strClipSeq.length() - 3))
            {
                //Here save softclip border reads
                int iClipPosInRef = al.Position + stClipReads.iQueryClipPos;
                En_ClipPart enClipPart = cpMax;
                if(stClipReads.iQueryClipPos < 5)
                {
                    enClipPart = cpLeft;
                }
                else if(stClipReads.iQueryClipPos > (al.QueryBases.length() - strClipSeq.length() - 5))
                {
                    enClipPart = cpRight;
                }

                //Save Border reads
                stBorderSCR.enClipPart = enClipPart;
                stBorderSCR.iClipPos = iClipPosInRef;
                stBorderSCR.iReadsMapPos = stClipReads.iReadsMapPos;
                stBorderSCR.strClipSeq = strClipSeq;
                stBorderSCR.strQuerySeq = stClipReads.strQuerySeq;
                stBorderSCR.iInsertSize = stClipReads.iInsertSize;
                stBorderSCR.strName = stClipReads.strName;
                stBorderSCR.iMatPos = stClipReads.iMateMapPos;
                stBorderSCR.bFirstMate = stClipReads.bFirstMate;
                stBorderSCR.bSecondMate = stClipReads.bSecondMate;
                stBorderSCR.bMateMapped = stClipReads.bMateMapped;

                stBorderSCR.iQueryMapPos = stClipReads.iQueryMapPos;
                stBorderSCR.iQueryMapLength = stClipReads.iQueryMapLength;
                stBorderSCR.iQueryClipPos = stClipReads.iQueryClipPos;
                stBorderSCR.iQueryClipLength = stClipReads.iQueryClipLength;

                //--> Check where we need to put
                bool bFindChrom = false;
                string strChrom = vFasta[al.RefID].GetBrifName();
                for(vector<St_ChromBorderSCR>::iterator itrChromSCR = vChromBorderSCR.begin();
                    itrChromSCR != vChromBorderSCR.end(); itrChromSCR++)
                {
                    if(itrChromSCR->strChrom == strChrom)
                    {
                        itrChromSCR->vSCR.push_back(stBorderSCR);
                        bFindChrom = true;
                        break;
                    }
                }
                //Is a new one in a new chromosome
                if(!bFindChrom)
                {
                    stChromSCR.Clear();
                    stChromSCR.strChrom = strChrom;
                    stChromSCR.vSCR.push_back(stBorderSCR);
                    vChromBorderSCR.push_back(stChromSCR);
                }
                //<-------
                //vBorderSCR.push_back(stBorderSCR);
                stBorderSCR.Clear();
            }
        }
    }

    if(m_pOfs == NULL)
    {
        cout << "vClipReads Size: " << to_string(vClipReads.size()) << endl;
        cout << "vChromBorderSCR Size: " << to_string(vChromBorderSCR.size()) << endl;
        //-->Output the number SCR for each different chromosome
        cout << "======= Details =======" << endl;
        for(vector<St_ChromBorderSCR>::iterator itrChromSCR = vChromBorderSCR.begin();
            itrChromSCR != vChromBorderSCR.end(); itrChromSCR++)
        {
            //Sort them first
            sort(itrChromSCR->vSCR.begin(), itrChromSCR->vSCR.end(), sort_scrClipPos_func);
            cout << itrChromSCR->strChrom << ": " << to_string(itrChromSCR->vSCR.size()) << endl;
        }
        cout << "=======================" << endl;
    }
    else
    {
        (*m_pOfs) << "vClipReads Size: " << to_string(vClipReads.size()) << endl;
        (*m_pOfs) << "vChromBorderSCR Size: " << to_string(vChromBorderSCR.size()) << endl;
        //-->Output the number SCR for each different chromosome
        (*m_pOfs) << "======= Details =======" << endl;
        for(vector<St_ChromBorderSCR>::iterator itrChromSCR = vChromBorderSCR.begin();
            itrChromSCR != vChromBorderSCR.end(); itrChromSCR++)
        {
            //Sort them first
            sort(itrChromSCR->vSCR.begin(), itrChromSCR->vSCR.end(), sort_scrClipPos_func);
            (*m_pOfs) << itrChromSCR->strChrom << ": " << to_string(itrChromSCR->vSCR.size()) << endl;
        }
        (*m_pOfs) << "=======================" << endl;
    }
    //<----

//    for(vector<St_ClipReads>::iterator itr = vClipReads.begin(); itr != vClipReads.end(); itr++)
//    {
//        //Pick the clipped reads in the boundary
//        bool bBorder = false;

//        int iIndex = 0;
//        string strClipSeq = "";
//        int iClipPosInRef = -1;
//        for(vector<int>::iterator itrSCPos = itr->vPos.begin(); itrSCPos != itr->vPos.end();
//            itrSCPos++, iIndex++)
//        {
//            //Check if it is the last clip part from the right side
//            //int iRightOffSet = *itrSCPos - (itr->strQuerySeq.length() - itr->vClipSeq[iIndex].length());
//            if(abs(*itrSCPos) <= 10)  //Left side
//            {
//                enClipPart = cpLeft;
//                strClipSeq = itr->vClipSeq[iIndex];
//                iClipPosInRef = itr->iReadsMapPos + *itrSCPos;
//                bBorder = true;
//                break;
//            }
//            else //if(abs(iRightOffSet) <= 10) //  This is right side
//            {
//                enClipPart = cpRight;
//                strClipSeq = itr->vClipSeq[iIndex];
//                iClipPosInRef = itr->iReadsMapPos + itr->strQuerySeq.length() -
//                                itr->vClipSeq[iIndex].length();
//                bBorder = true;
//                break;
//            }
//        }

//        if(bBorder)
//        {
//            //Save Border reads
//            stBorderSCR.enClipPart = enClipPart;
//            stBorderSCR.iClipPos = iClipPosInRef;
//            stBorderSCR.iReadsMapPos = itr->iReadsMapPos;
//            stBorderSCR.strClipSeq = itr->vClipSeq[iIndex];
//            stBorderSCR.strQuerySeq = itr->strQuerySeq;
//            stBorderSCR.iInsertSize = itr->iInsertSize;
//            stBorderSCR.strName = itr->strName;
//            stBorderSCR.iMatPos = itr->iMateMapPos;
//            stBorderSCR.bFirstMate = itr->bFirstMate;
//            stBorderSCR.bSecondMate = itr->bSecondMate;
//            stBorderSCR.bMateMapped = itr->bMateMapped;
//            vBorderSCR.push_back(stBorderSCR);

//            stBorderSCR.Clear();
//        }
//    }

    //sort soft clip reads by clip position  from small to large
    if(m_pOfs == NULL)
    {
        cout << "Start_sort_borderSCR -->" << endl;
    }
    else
    {
        (*m_pOfs) << "Start_sort_borderSCR -->" << endl;
    }

//    for(vector<St_BorderSCR>::iterator itr = vBorderSCR.begin();  itr != vBorderSCR.end(); itr++)
//    {
//        cout << to_string(itr->iClipPos) << "\t \t";
//    }
//    cout << endl;
    //sort(vBorderSCR.begin(), vBorderSCR.end(), sort_scrClipPos_func);
    //cout << "Sorted vBorderSCR Size: " << to_string(vBorderSCR.size()) << endl;

    delete pBamReader;
    pBamReader = NULL;
}
//<----

void ClsParseBam::GetDiscordantReads(string strBamFilePath, vector<St_Fasta>& vFasta,
                                     vector<St_ChromDiscordantReads>& vChromDR)
{
    /* How to do it:
     * 1: Read Bam
     * 2: Check if mapped
     * 3: Check discordant --> absolution value (we allow inversion)
     * 4: Record itself and itsmate
     * 5: Check duplicate
     * 6: End
     */

    //Go!!!

    St_ChromDiscordantReads stChromDR;
    vChromDR.clear();

    //parse bam file
    BamReader* pBamReader = new BamReader();

    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");
    BamAlignment al;
    St_DiscordantReads stDR;

    //vDiscdRreads.clear();

    float fMinThreshold = m_fMeanInsertSize - 3 * m_fStdDevInsertSize;
    float fMaxThreshold = m_fMeanInsertSize + 3 * m_fStdDevInsertSize;

    if(m_fStdDevInsertSize < 50)
    {
        fMinThreshold = m_fMeanInsertSize - 6 * m_fStdDevInsertSize;
        fMaxThreshold = m_fMeanInsertSize + 6 * m_fStdDevInsertSize;
    }

    while(pBamReader->GetNextAlignment(al)) //Get each alignment result
    {
        if(!al.IsMapped() || !al.IsMateMapped() || al.QueryBases == "") //　we need both map
            continue;

        if(al.IsDuplicate() || al.IsFailedQC() || !al.IsPrimaryAlignment() || al.IsReverseStrand())
            continue;

        if(al.QueryBases.find('N') != string::npos ||
           al.QueryBases.find('n') != string::npos)
            continue;

        // skip concordant reads
        int iRealInsertSize = al.InsertSize;
        if(al.IsSecondMate() && abs(iRealInsertSize) > 100*2)
            iRealInsertSize = iRealInsertSize * -1;

//        if(iRealInsertSize == 0) //这个情况属于认为的把mate去掉了,但是在bam file里面mate position是有记录的
//        {
//            iRealInsertSize = al.MatePosition - al.Position;
//            if(iRealInsertSize > 0)
//                iRealInsertSize += al.QueryBases.length();
//            else
//                iRealInsertSize -= al.QueryBases.length();
//        }

        if(//abs(iRealInsertSize) > iMinThreshold &&
           abs(iRealInsertSize) < fMaxThreshold)
            continue;

        //skip the case with too large insertion size --> (this too large means abnormal) -->
        int iAbnormalLarge = 100000;
        if(abs(iRealInsertSize) > iAbnormalLarge)
            continue;
        //<--

        //-->Get the name of current reference (chromosome)
        string strChrom = vFasta[al.RefID].GetBrifName();

        //Now they are both mapped discordant reads
        if(//iRealInsertSize >= 0 ||
           al.MatePosition > al.Position) // current reads is left, the mate is right
        {
            //First check if it has been recorded -->

            //-->Record alignment result
            string strMapSeq = al.AlignedBases;
            string strClipSeq = "";
            int iQueryClipPos = -1;
            int iQueryMapPos = -1;
            int iOffSet = 0;
            //<--

            bool bFind = false;
            St_DiscordantReads* pDR = NULL;
            for(vector<St_ChromDiscordantReads>::iterator itrChromDR = vChromDR.begin();
                itrChromDR != vChromDR.end(); itrChromDR++)
            {
                if(itrChromDR->strChrom != strChrom)
                    continue;

                //vector<St_DiscordantReads>::iterator itrTarget = itrChromDR->vDR.begin();
                for(vector<St_DiscordantReads>::iterator itr = itrChromDR->vDR.begin();
                    itr != itrChromDR->vDR.end(); itr++)
                {
                    if(itr->iReadsPos == al.Position &&
                       itr->iMatePos == al.MatePosition)
                    {
                        pDR = &(*itr);
                        bFind = true;
                        break;
                    }
                }
                break;
            }

            if(bFind)
            {
                //update existed value
                if(pDR->strReadsAlignSeq == "")
                {
                    pDR->strReadsAlignSeq = al.AlignedBases;
                    //Check if it contains clip part -->
                    for(std::vector<CigarOp>::iterator itrCigar = al.CigarData.begin();
                        itrCigar != al.CigarData.end(); itrCigar++)
                    {
                        if(itrCigar->Type == 'M')
                        {
                            //if(strMapSeq == "")
                            //    strMapSeq += al.QueryBases.substr(iOffSet, itr->Length);
                            if(iQueryMapPos == -1)
                                iQueryMapPos = iOffSet;
                            iOffSet += itrCigar->Length;
                        }

                        if(itrCigar->Type == 'S')
                        {
                            if(strClipSeq == "")
                                strClipSeq += al.QueryBases.substr(iOffSet, itrCigar->Length);
                            if(iQueryClipPos == -1)
                                iQueryClipPos = iOffSet;
                            iOffSet += itrCigar->Length;
                            pDR->bClip = true;
                        }
                    }

                    //Update the mapping info
                    pDR->iQueryMapPos = iQueryMapPos;
                    pDR->iQueryMapLength = strMapSeq.length();
                    if(pDR->bClip)
                    {
                        pDR->iQueryClipPos = iQueryClipPos;
                        pDR->iQueryClipLength = strClipSeq.length();
                        pDR->strClipSeq = strClipSeq;
                    }
                    //<--
                }
            }
            else
            {
                //If it is new
                stDR.Clear();
                stDR.iReadsPos = al.Position;
                stDR.iMatePos = al.MatePosition;
                stDR.strReadsAlignSeq = al.AlignedBases;
                stDR.iInsertSize = iRealInsertSize;
                //Check if it contains clip part -->
                for(std::vector<CigarOp>::iterator itrCigar = al.CigarData.begin();
                    itrCigar != al.CigarData.end(); itrCigar++)
                {
                    if(itrCigar->Type == 'M')
                    {
                        //if(strMapSeq == "")
                        //    strMapSeq += al.QueryBases.substr(iOffSet, itr->Length);
                        if(iQueryMapPos == -1)
                            iQueryMapPos = iOffSet;
                        iOffSet += itrCigar->Length;
                    }

                    if(itrCigar->Type == 'S')
                    {
                        if(strClipSeq == "")
                            strClipSeq += al.QueryBases.substr(iOffSet, itrCigar->Length);
                        if(iQueryClipPos == -1)
                            iQueryClipPos = iOffSet;
                        iOffSet += itrCigar->Length;
                        stDR.bClip = true;
                    }
                }
                //Update the mapping info
                stDR.iQueryMapPos = iQueryMapPos;
                stDR.iQueryMapLength = strMapSeq.length();
                if(stDR.bClip)
                {
                    stDR.iQueryClipPos = iQueryClipPos;
                    stDR.iQueryClipLength = strClipSeq.length();
                    stDR.strClipSeq = strClipSeq;
                }
                //<--
                //<--
                bool bFindChromDR = false;
                for(vector<St_ChromDiscordantReads>::iterator itr = vChromDR.begin();
                    itr != vChromDR.end(); itr++)
                {
                    if(itr->strChrom == strChrom)
                    {
                        itr->vDR.push_back(stDR);
                        bFindChromDR = true;
                        break;
                    }
                }
                if(!bFindChromDR)
                {
                    stChromDR.Clear();
                    stChromDR.strChrom = strChrom;
                    stChromDR.vDR.push_back(stDR);
                    vChromDR.push_back(stChromDR);
                }
                //vDiscdRreads.push_back(stDR);
                //<--
            }
        }
        else //negative insert size !!!
        {           
            //-->Record alignment result
            string strMapSeq = al.AlignedBases;
            string strClipSeq = "";
            int iQueryClipPos = -1;
            int iQueryMapPos = -1;
            int iOffSet = 0;
            //<--

            //First check if it has been recorded -->
            bool bFind = false;

            St_DiscordantReads* pDR = NULL;
            for(vector<St_ChromDiscordantReads>::iterator itrChromDR = vChromDR.begin();
                itrChromDR != vChromDR.end(); itrChromDR++)
            {
                if(itrChromDR->strChrom != strChrom)
                    continue;

                //vector<St_DiscordantReads>::iterator itrTarget = itrChromDR->vDR.begin();
                for(vector<St_DiscordantReads>::iterator itr = itrChromDR->vDR.begin();
                    itr != itrChromDR->vDR.end(); itr++)
                {
                    if(itr->iReadsPos == al.MatePosition &&
                       itr->iMatePos == al.Position)
                    {
                        pDR = &(*itr);
                        bFind = true;
                        break;
                    }
                }
                break;
            }

            if(bFind)
            {
                //update existed value
                if(pDR->strMateAlignSeq == "")
                {
                    pDR->strMateAlignSeq = al.AlignedBases;
                    //Check if it contains clip part -->
                    for(std::vector<CigarOp>::iterator itrCigar = al.CigarData.begin();
                        itrCigar != al.CigarData.end(); itrCigar++)
                    {
                        if(itrCigar->Type == 'M')
                        {
                            //if(strMapSeq == "")
                            //    strMapSeq += al.QueryBases.substr(iOffSet, itr->Length);
                            if(iQueryMapPos == -1)
                                iQueryMapPos = iOffSet;
                            iOffSet += itrCigar->Length;
                        }

                        if(itrCigar->Type == 'S')
                        {
                            if(strClipSeq == "")
                                strClipSeq += al.QueryBases.substr(iOffSet, itrCigar->Length);
                            if(iQueryClipPos == -1)
                                iQueryClipPos = iOffSet;
                            iOffSet += itrCigar->Length;
                            pDR->bMateClip = true;
                        }
                    }

                    //Update the mapping info
                    pDR->iMatQueryMapPos = iQueryMapPos;
                    pDR->iMatQueryMapLength = strMapSeq.length();
                    if(pDR->bMateClip)
                    {
                        pDR->iMatQueryClipPos = iQueryClipPos;
                        pDR->iMatQueryClipLength = strClipSeq.length();
                        pDR->strMatClipSeq = strClipSeq;
                    }
                    //<--
                }
            }
            else
            {
                //If it is new
                stDR.Clear();
                stDR.iReadsPos = al.MatePosition;
                stDR.iMatePos = al.Position;
                stDR.strMateAlignSeq = al.AlignedBases;
                stDR.iInsertSize = iRealInsertSize; //abs(iRealInsertSize);


                for(std::vector<CigarOp>::iterator itrCigar = al.CigarData.begin();
                    itrCigar != al.CigarData.end(); itrCigar++)
                {
                    if(itrCigar->Type == 'M')
                    {
                        //if(strMapSeq == "")
                        //    strMapSeq += al.QueryBases.substr(iOffSet, itr->Length);
                        if(iQueryMapPos == -1)
                            iQueryMapPos = iOffSet;
                        iOffSet += itrCigar->Length;
                    }

                    if(itrCigar->Type == 'S')
                    {
                        if(strClipSeq == "")
                            strClipSeq += al.QueryBases.substr(iOffSet, itrCigar->Length);
                        if(iQueryClipPos == -1)
                            iQueryClipPos = iOffSet;
                        iOffSet += itrCigar->Length;
                        stDR.bMateClip = true;
                    }
                }

                //Update the mapping info
                stDR.iMatQueryMapPos = iQueryMapPos;
                stDR.iMatQueryMapLength = strMapSeq.length();
                if(stDR.bMateClip)
                {
                    stDR.iMatQueryClipPos = iQueryClipPos;
                    stDR.iMatQueryClipLength = strClipSeq.length();
                    stDR.strMatClipSeq = strClipSeq;
                }
                //<--

//                //Check if it contains clip part -->
//                for(std::vector<CigarOp>::iterator itrCigar = al.CigarData.begin();
//                    itrCigar != al.CigarData.end(); itrCigar++)
//                {
//                    if(itrCigar->Type == 'S')
//                    {
//                        stDR.bMateClip = true;
//                        break;
//                    }
//                }
//                //<--

                bool bFindChromDR = false;
                for(vector<St_ChromDiscordantReads>::iterator itr = vChromDR.begin();
                    itr != vChromDR.end(); itr++)
                {
                    if(itr->strChrom == strChrom)
                    {
                        itr->vDR.push_back(stDR);
                        bFindChromDR = true;
                        break;
                    }
                }
                if(!bFindChromDR)
                {
                    stChromDR.Clear();
                    stChromDR.strChrom = strChrom;
                    stChromDR.vDR.push_back(stDR);
                    vChromDR.push_back(stChromDR);
                }
                //vDiscdRreads.push_back(stDR);
            }
        }
    }

    //output the statistic of what we get
    if(m_pOfs == NULL)
    {
        cout << "vChromDR Size: " << to_string(vChromDR.size()) << endl;
        cout << "====== Details ======" << endl;
        for(vector<St_ChromDiscordantReads>::iterator itr = vChromDR.begin();
            itr != vChromDR.end(); itr++)
        {
            cout << itr->strChrom << ": " <<  to_string(itr->vDR.size()) << endl;
        }
        cout << "=====================" << endl;
    }
    else
    {
        (*m_pOfs) << "vChromDR Size: " << to_string(vChromDR.size()) << endl;
        (*m_pOfs) << "====== Details ======" << endl;
        for(vector<St_ChromDiscordantReads>::iterator itr = vChromDR.begin();
            itr != vChromDR.end(); itr++)
        {
            (*m_pOfs) << itr->strChrom << ": " <<  to_string(itr->vDR.size()) << endl;
        }
        (*m_pOfs) << "=====================" << endl;
    }

    //cout << "Discordant Reads Size: " << to_string(vDiscdRreads.size()) << endl;

    delete pBamReader;
    pBamReader = NULL;
}

void ClsParseBam::GetSCRAndDR(string strBamFilePath, vector<St_Fasta>& vFasta,
                              vector<St_ChromBorderSCR>& vChromBorderSCR,
                              vector<St_ChromDiscordantReads>& vChromDR)
{
    //We do them together --> check each al
    ///--> For Soft Clip Reads
    vector<St_ClipReads> vClipReads;
    St_ClipReads stClipReads;
    St_ChromBorderSCR stChromSCR;
    St_BorderSCR stBorderSCR;
    vChromBorderSCR.clear();

    ///--> For Discordant Reads
    St_DiscordantReads stDR;
    St_ChromDiscordantReads stChromDR;
    vChromDR.clear();
    float fMinThreshold = m_fMeanInsertSize - 3 * m_fStdDevInsertSize;
    float fMaxThreshold = m_fMeanInsertSize + 3 * m_fStdDevInsertSize;
    if(m_fStdDevInsertSize < 50)
    {
        fMinThreshold = m_fMeanInsertSize - 6 * m_fStdDevInsertSize;
        fMaxThreshold = m_fMeanInsertSize + 6 * m_fStdDevInsertSize;
    }

    //->Do loop
    BamReader* pBamReader = new BamReader();
    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");
    BamAlignment al;

    while(pBamReader->GetNextAlignment(al)) //Get each alignment result
    {
        GetSCRByAl(al, vClipReads, stClipReads, stChromSCR, stBorderSCR,
                   vFasta, vChromBorderSCR);

        GetDRByAl(al, stDR, stChromDR,
                  vFasta, vChromDR,
                  fMinThreshold, fMaxThreshold);
    }
    delete pBamReader;
    pBamReader = NULL;

    //PrintTestingInfoSCR(vClipReads, vChromBorderSCR);
    //PrintTestingInfoDR(vChromDR);
}

void ClsParseBam::GetSCRByAl(BamAlignment& al,
                             vector<St_ClipReads>& vClipReads, St_ClipReads& stClipReads,
                             St_ChromBorderSCR& stChromSCR, St_BorderSCR& stBorderSCR,
                             vector<St_Fasta>& vFasta, vector<St_ChromBorderSCR>& vChromBorderSCR)
{
    if(!al.IsMapped() || al.QueryBases == "")
        return;

    if(al.IsDuplicate() || al.IsFailedQC() || !al.IsPrimaryAlignment()) //We should keep the case: al.IsReverseStrand()
        return;

    if(al.QueryBases.find('N') != string::npos ||
       al.QueryBases.find('n') != string::npos)
        return;

//        if(al.Position == 132558874)
//        {
//            cout << "I am the target!" << endl;
//            if(al.IsFirstMate())
//                cout << "I am the first mate!" << endl;

//            if(al.IsSecondMate())
//                cout << "I am the second mate!" << endl;

//            if(al.IsDuplicate())
//                cout << "I am duplicate!" << endl;
//        }

    bool bFindClip = false;
    int iOffSet = 0;

//        //----> Just For testing
//        if(al.IsSecondMate())
//        {
//            string strMapSeq = "";
//            string strClipSeq = "";
//            bool bClipped = false;
//            int iClipPos = 0;
//            for(std::vector<CigarOp>::iterator itr = al.CigarData.begin(); itr != al.CigarData.end(); itr++)
//            {
//                switch(itr->Type)
//                {
//                    case 'M': // alignment match (can be a sequence match or mismatch)
//                        strMapSeq += al.QueryBases.substr(iOffSet, itr->Length);
//                        iOffSet += itr->Length;
//                        break;
//                    case 'I': // insertion to the reference
//                        break;
//                    case 'D': // deletion from the reference
//                    case 'N':  // skipped region from the reference
//                        break;
//                    case 'S':  // soft clipping (clipped sequences present in SEQ)
//                        strClipSeq += al.QueryBases.substr(iOffSet, itr->Length);
//                        iClipPos = iOffSet;
//                        iOffSet += itr->Length;
//                        bClipped = true;
//                        break;
//                    case 'H':  // hard clipping (clipped sequences NOT present in SEQ)
//                        break;
//                    case 'P': // padding (silent deletion from padded reference)
//                    case '=': // sequence match
//                    case 'X': // sequence mismatch
//                        break;
//                }
//            }
//            if(bClipped)
//            {
//                cout << "Map Position : " << to_string(al.Position) << endl;
//                cout << "Clip Position: " << to_string(iClipPos) << endl;
//                cout << "Map Seq      : " << strMapSeq << endl;
//                cout << "Clip Seq     : " << strClipSeq << endl;
//                cout << "Query Seq    : " << al.QueryBases << endl;
//                cout << endl;
//            }
//        }
//        //<----
    iOffSet = 0;

    stClipReads.Clear();
    string strMapSeq = al.AlignedBases;
    string strClipSeq = "";
    int iQueryClipPos = -1;
    int iQueryMapPos = -1;
    for(std::vector<CigarOp>::iterator itr = al.CigarData.begin(); itr != al.CigarData.end(); itr++)
    {
        switch(itr->Type)
        {
            case 'M': // alignment match (can be a sequence match or mismatch)
                //if(strMapSeq == "")
                //    strMapSeq += al.QueryBases.substr(iOffSet, itr->Length);
                if(iQueryMapPos == -1)
                    iQueryMapPos = iOffSet;
                iOffSet += itr->Length;
                break;
            case 'I': // insertion to the reference
                break;
            case 'D': // deletion from the reference
            case 'N':  // skipped region from the reference
                break;
            case 'S':  // soft clipping (clipped sequences present in SEQ)
                if(strClipSeq == "")
                    strClipSeq += al.QueryBases.substr(iOffSet, itr->Length);
                stClipReads.vClipSeq.push_back(al.QueryBases.substr(iOffSet, itr->Length));
                stClipReads.vPos.push_back(iOffSet);

                if(iQueryClipPos == -1)
                    iQueryClipPos = iOffSet;
                iOffSet += itr->Length;
                bFindClip = true;

                break;
            case 'H':  // hard clipping (clipped sequences NOT present in SEQ)
                break;
            case 'P': // padding (silent deletion from padded reference)
            case '=': // sequence match
            case 'X': // sequence mismatch
                break;
        }
    }

    if(bFindClip)
    {
        if(strMapSeq.length() < 15 || strClipSeq.length() < 15)
            return;

        int iRealInsertSize = al.InsertSize;
        if(al.IsSecondMate() &&
           abs(iRealInsertSize) > 100*2)
            iRealInsertSize = iRealInsertSize * -1;

//            if(iRealInsertSize == 0 && al.IsMateMapped()) //  this is belongs to the special case
//            {
//                iRealInsertSize = al.MatePosition - al.Position;
//                if(iRealInsertSize > 0)
//                    iRealInsertSize += al.QueryBases.length();
//                else
//                    iRealInsertSize -= al.QueryBases.length();
//            }

        stClipReads.strName = al.Name;
        stClipReads.strQuerySeq = al.QueryBases;
        stClipReads.iReadsMapPos = al.Position;
        stClipReads.iInsertSize = iRealInsertSize;
        stClipReads.iMateMapPos = al.MatePosition;
        stClipReads.bFirstMate = al.IsFirstMate();
        stClipReads.bSecondMate = al.IsSecondMate();
        stClipReads.bMateMapped = al.IsMateMapped();

        stClipReads.iQueryMapPos = iQueryMapPos;
        stClipReads.iQueryMapLength = strMapSeq.length();
        stClipReads.iQueryClipPos = iQueryClipPos;
        stClipReads.iQueryClipLength = strClipSeq.length();

        //Save clipped info
        vClipReads.push_back(stClipReads);

        //The number of clip which do not happed in the beginning or the ending part of the reads -->
        if(stClipReads.iQueryClipPos < 3 ||
           stClipReads.iQueryClipPos > (al.QueryBases.length() - strClipSeq.length() - 3))
        {
            //Here save softclip border reads
            int iClipPosInRef = al.Position + stClipReads.iQueryClipPos;
            En_ClipPart enClipPart = cpMax;
            if(stClipReads.iQueryClipPos < 5)
            {
                enClipPart = cpLeft;
            }
            else if(stClipReads.iQueryClipPos > (al.QueryBases.length() - strClipSeq.length() - 5))
            {
                enClipPart = cpRight;
            }

            //Save Border reads
            stBorderSCR.enClipPart = enClipPart;
            stBorderSCR.iClipPos = iClipPosInRef;
            stBorderSCR.iReadsMapPos = stClipReads.iReadsMapPos;
            stBorderSCR.strClipSeq = strClipSeq;
            stBorderSCR.strQuerySeq = stClipReads.strQuerySeq;
            stBorderSCR.iInsertSize = stClipReads.iInsertSize;
            stBorderSCR.strName = stClipReads.strName;
            stBorderSCR.iMatPos = stClipReads.iMateMapPos;
            stBorderSCR.bFirstMate = stClipReads.bFirstMate;
            stBorderSCR.bSecondMate = stClipReads.bSecondMate;
            stBorderSCR.bMateMapped = stClipReads.bMateMapped;

            stBorderSCR.iQueryMapPos = stClipReads.iQueryMapPos;
            stBorderSCR.iQueryMapLength = stClipReads.iQueryMapLength;
            stBorderSCR.iQueryClipPos = stClipReads.iQueryClipPos;
            stBorderSCR.iQueryClipLength = stClipReads.iQueryClipLength;

            //--> Check where we need to put
            bool bFindChrom = false;
            string strChrom = vFasta[al.RefID].GetBrifName();
            for(vector<St_ChromBorderSCR>::iterator itrChromSCR = vChromBorderSCR.begin();
                itrChromSCR != vChromBorderSCR.end(); itrChromSCR++)
            {
                if(itrChromSCR->strChrom == strChrom)
                {
                    itrChromSCR->vSCR.push_back(stBorderSCR);
                    bFindChrom = true;
                    break;
                }
            }
            //Is a new one in a new chromosome
            if(!bFindChrom)
            {
                stChromSCR.Clear();
                stChromSCR.strChrom = strChrom;
                stChromSCR.vSCR.push_back(stBorderSCR);
                vChromBorderSCR.push_back(stChromSCR);
            }
            //<-------
            //vBorderSCR.push_back(stBorderSCR);
            stBorderSCR.Clear();
        }
    }
}

void ClsParseBam::PrintTestingInfoSCR(vector<St_ClipReads>& vClipReads,
                                      vector<St_ChromBorderSCR>& vChromBorderSCR)
{
    if(m_pOfs == NULL)
    {
        cout << "vClipReads Size: " << to_string(vClipReads.size()) << endl;
        cout << "vChromBorderSCR Size: " << to_string(vChromBorderSCR.size()) << endl;
        //-->Output the number SCR for each different chromosome
        cout << "======= Details =======" << endl;
        for(vector<St_ChromBorderSCR>::iterator itrChromSCR = vChromBorderSCR.begin();
            itrChromSCR != vChromBorderSCR.end(); itrChromSCR++)
        {
            //Sort them first
            sort(itrChromSCR->vSCR.begin(), itrChromSCR->vSCR.end(), sort_scrClipPos_func);
            cout << itrChromSCR->strChrom << ": " << to_string(itrChromSCR->vSCR.size()) << endl;
        }
        cout << "=======================" << endl;
    }
    else
    {
        (*m_pOfs) << "vClipReads Size: " << to_string(vClipReads.size()) << endl;
        (*m_pOfs) << "vChromBorderSCR Size: " << to_string(vChromBorderSCR.size()) << endl;
        //-->Output the number SCR for each different chromosome
        (*m_pOfs) << "======= Details =======" << endl;
        for(vector<St_ChromBorderSCR>::iterator itrChromSCR = vChromBorderSCR.begin();
            itrChromSCR != vChromBorderSCR.end(); itrChromSCR++)
        {
            //Sort them first
            sort(itrChromSCR->vSCR.begin(), itrChromSCR->vSCR.end(), sort_scrClipPos_func);
            (*m_pOfs) << itrChromSCR->strChrom << ": " << to_string(itrChromSCR->vSCR.size()) << endl;
        }
        (*m_pOfs) << "=======================" << endl;
    }
    //<----

//    for(vector<St_ClipReads>::iterator itr = vClipReads.begin(); itr != vClipReads.end(); itr++)
//    {
//        //Pick the clipped reads in the boundary
//        bool bBorder = false;

//        int iIndex = 0;
//        string strClipSeq = "";
//        int iClipPosInRef = -1;
//        for(vector<int>::iterator itrSCPos = itr->vPos.begin(); itrSCPos != itr->vPos.end();
//            itrSCPos++, iIndex++)
//        {
//            //Check if it is the last clip part from the right side
//            //int iRightOffSet = *itrSCPos - (itr->strQuerySeq.length() - itr->vClipSeq[iIndex].length());
//            if(abs(*itrSCPos) <= 10)  //Left side
//            {
//                enClipPart = cpLeft;
//                strClipSeq = itr->vClipSeq[iIndex];
//                iClipPosInRef = itr->iReadsMapPos + *itrSCPos;
//                bBorder = true;
//                break;
//            }
//            else //if(abs(iRightOffSet) <= 10) //  This is right side
//            {
//                enClipPart = cpRight;
//                strClipSeq = itr->vClipSeq[iIndex];
//                iClipPosInRef = itr->iReadsMapPos + itr->strQuerySeq.length() -
//                                itr->vClipSeq[iIndex].length();
//                bBorder = true;
//                break;
//            }
//        }

//        if(bBorder)
//        {
//            //Save Border reads
//            stBorderSCR.enClipPart = enClipPart;
//            stBorderSCR.iClipPos = iClipPosInRef;
//            stBorderSCR.iReadsMapPos = itr->iReadsMapPos;
//            stBorderSCR.strClipSeq = itr->vClipSeq[iIndex];
//            stBorderSCR.strQuerySeq = itr->strQuerySeq;
//            stBorderSCR.iInsertSize = itr->iInsertSize;
//            stBorderSCR.strName = itr->strName;
//            stBorderSCR.iMatPos = itr->iMateMapPos;
//            stBorderSCR.bFirstMate = itr->bFirstMate;
//            stBorderSCR.bSecondMate = itr->bSecondMate;
//            stBorderSCR.bMateMapped = itr->bMateMapped;
//            vBorderSCR.push_back(stBorderSCR);

//            stBorderSCR.Clear();
//        }
//    }

    //sort soft clip reads by clip position  from small to large
    if(m_pOfs == NULL)
    {
        cout << "Start_sort_borderSCR -->" << endl;
    }
    else
    {
        (*m_pOfs) << "Start_sort_borderSCR -->" << endl;
    }

//    for(vector<St_BorderSCR>::iterator itr = vBorderSCR.begin();  itr != vBorderSCR.end(); itr++)
//    {
//        cout << to_string(itr->iClipPos) << "\t \t";
//    }
//    cout << endl;
    //sort(vBorderSCR.begin(), vBorderSCR.end(), sort_scrClipPos_func);
    //cout << "Sorted vBorderSCR Size: " << to_string(vBorderSCR.size()) << endl;
}

void ClsParseBam::GetDRByAl(BamAlignment& al,
                            St_DiscordantReads& stDR, St_ChromDiscordantReads& stChromDR,
                            vector<St_Fasta>& vFasta, vector<St_ChromDiscordantReads>& vChromDR,
                            float& fMinThreshold, float& fMaxThreshold)
{
    if(!al.IsMapped() || !al.IsMateMapped() || al.QueryBases == "") //　we need both map
        return;

    if(al.IsDuplicate() || al.IsFailedQC() || !al.IsPrimaryAlignment() || al.IsReverseStrand())
        return;

    if(al.QueryBases.find('N') != string::npos ||
       al.QueryBases.find('n') != string::npos)
        return;

    // skip concordant reads
    int iRealInsertSize = al.InsertSize;
    if(al.IsSecondMate() && abs(iRealInsertSize) > 100*2)
        iRealInsertSize = iRealInsertSize * -1;

//        if(iRealInsertSize == 0) //这个情况属于认为的把mate去掉了,但是在bam file里面mate position是有记录的
//        {
//            iRealInsertSize = al.MatePosition - al.Position;
//            if(iRealInsertSize > 0)
//                iRealInsertSize += al.QueryBases.length();
//            else
//                iRealInsertSize -= al.QueryBases.length();
//        }

    if(//abs(iRealInsertSize) > iMinThreshold &&
       abs(iRealInsertSize) < fMaxThreshold)
        return;

    //skip the case with too large insertion size --> (this too large means abnormal) -->
    int iAbnormalLarge = 100000;
    if(abs(iRealInsertSize) > iAbnormalLarge)
        return;
    //<--

    //-->Get the name of current reference (chromosome)
    string strChrom = vFasta[al.RefID].GetBrifName();

    //Now they are both mapped discordant reads
    if(//iRealInsertSize >= 0 ||
       al.MatePosition > al.Position) // current reads is left, the mate is right
    {
        //First check if it has been recorded -->

        //-->Record alignment result
        string strMapSeq = al.AlignedBases;
        string strClipSeq = "";
        int iQueryClipPos = -1;
        int iQueryMapPos = -1;
        int iOffSet = 0;
        //<--

        bool bFind = false;
        St_DiscordantReads* pDR = NULL;
        for(vector<St_ChromDiscordantReads>::iterator itrChromDR = vChromDR.begin();
            itrChromDR != vChromDR.end(); itrChromDR++)
        {
            if(itrChromDR->strChrom != strChrom)
                continue;

            //vector<St_DiscordantReads>::iterator itrTarget = itrChromDR->vDR.begin();
            for(vector<St_DiscordantReads>::iterator itr = itrChromDR->vDR.begin();
                itr != itrChromDR->vDR.end(); itr++)
            {
                if(itr->iReadsPos == al.Position &&
                   itr->iMatePos == al.MatePosition)
                {
                    pDR = &(*itr);
                    bFind = true;
                    break;
                }
            }
            break;
        }

        if(bFind)
        {
            //update existed value
            if(pDR->strReadsAlignSeq == "")
            {
                pDR->strReadsAlignSeq = al.AlignedBases;
                //Check if it contains clip part -->
                for(std::vector<CigarOp>::iterator itrCigar = al.CigarData.begin();
                    itrCigar != al.CigarData.end(); itrCigar++)
                {
                    if(itrCigar->Type == 'M')
                    {
                        //if(strMapSeq == "")
                        //    strMapSeq += al.QueryBases.substr(iOffSet, itr->Length);
                        if(iQueryMapPos == -1)
                            iQueryMapPos = iOffSet;
                        iOffSet += itrCigar->Length;
                    }

                    if(itrCigar->Type == 'S')
                    {
                        if(strClipSeq == "")
                            strClipSeq += al.QueryBases.substr(iOffSet, itrCigar->Length);
                        if(iQueryClipPos == -1)
                            iQueryClipPos = iOffSet;
                        iOffSet += itrCigar->Length;
                        pDR->bClip = true;
                    }
                }

                //Update the mapping info
                pDR->iQueryMapPos = iQueryMapPos;
                pDR->iQueryMapLength = strMapSeq.length();
                if(pDR->bClip)
                {
                    pDR->iQueryClipPos = iQueryClipPos;
                    pDR->iQueryClipLength = strClipSeq.length();
                    pDR->strClipSeq = strClipSeq;
                }
                //<--
            }
        }
        else
        {
            //If it is new
            stDR.Clear();
            stDR.iReadsPos = al.Position;
            stDR.iMatePos = al.MatePosition;
            stDR.strReadsAlignSeq = al.AlignedBases;
            stDR.iInsertSize = iRealInsertSize;
            //Check if it contains clip part -->
            for(std::vector<CigarOp>::iterator itrCigar = al.CigarData.begin();
                itrCigar != al.CigarData.end(); itrCigar++)
            {
                if(itrCigar->Type == 'M')
                {
                    //if(strMapSeq == "")
                    //    strMapSeq += al.QueryBases.substr(iOffSet, itr->Length);
                    if(iQueryMapPos == -1)
                        iQueryMapPos = iOffSet;
                    iOffSet += itrCigar->Length;
                }

                if(itrCigar->Type == 'S')
                {
                    if(strClipSeq == "")
                        strClipSeq += al.QueryBases.substr(iOffSet, itrCigar->Length);
                    if(iQueryClipPos == -1)
                        iQueryClipPos = iOffSet;
                    iOffSet += itrCigar->Length;
                    stDR.bClip = true;
                }
            }
            //Update the mapping info
            stDR.iQueryMapPos = iQueryMapPos;
            stDR.iQueryMapLength = strMapSeq.length();
            if(stDR.bClip)
            {
                stDR.iQueryClipPos = iQueryClipPos;
                stDR.iQueryClipLength = strClipSeq.length();
                stDR.strClipSeq = strClipSeq;
            }
            //<--
            //<--
            bool bFindChromDR = false;
            for(vector<St_ChromDiscordantReads>::iterator itr = vChromDR.begin();
                itr != vChromDR.end(); itr++)
            {
                if(itr->strChrom == strChrom)
                {
                    itr->vDR.push_back(stDR);
                    bFindChromDR = true;
                    break;
                }
            }
            if(!bFindChromDR)
            {
                stChromDR.Clear();
                stChromDR.strChrom = strChrom;
                stChromDR.vDR.push_back(stDR);
                vChromDR.push_back(stChromDR);
            }
            //vDiscdRreads.push_back(stDR);
            //<--
        }
    }
    else //negative insert size !!!
    {
        //-->Record alignment result
        string strMapSeq = al.AlignedBases;
        string strClipSeq = "";
        int iQueryClipPos = -1;
        int iQueryMapPos = -1;
        int iOffSet = 0;
        //<--

        //First check if it has been recorded -->
        bool bFind = false;

        St_DiscordantReads* pDR = NULL;
        for(vector<St_ChromDiscordantReads>::iterator itrChromDR = vChromDR.begin();
            itrChromDR != vChromDR.end(); itrChromDR++)
        {
            if(itrChromDR->strChrom != strChrom)
                continue;

            //vector<St_DiscordantReads>::iterator itrTarget = itrChromDR->vDR.begin();
            for(vector<St_DiscordantReads>::iterator itr = itrChromDR->vDR.begin();
                itr != itrChromDR->vDR.end(); itr++)
            {
                if(itr->iReadsPos == al.MatePosition &&
                   itr->iMatePos == al.Position)
                {
                    pDR = &(*itr);
                    bFind = true;
                    break;
                }
            }
            break;
        }

        if(bFind)
        {
            //update existed value
            if(pDR->strMateAlignSeq == "")
            {
                pDR->strMateAlignSeq = al.AlignedBases;
                //Check if it contains clip part -->
                for(std::vector<CigarOp>::iterator itrCigar = al.CigarData.begin();
                    itrCigar != al.CigarData.end(); itrCigar++)
                {
                    if(itrCigar->Type == 'M')
                    {
                        //if(strMapSeq == "")
                        //    strMapSeq += al.QueryBases.substr(iOffSet, itr->Length);
                        if(iQueryMapPos == -1)
                            iQueryMapPos = iOffSet;
                        iOffSet += itrCigar->Length;
                    }

                    if(itrCigar->Type == 'S')
                    {
                        if(strClipSeq == "")
                            strClipSeq += al.QueryBases.substr(iOffSet, itrCigar->Length);
                        if(iQueryClipPos == -1)
                            iQueryClipPos = iOffSet;
                        iOffSet += itrCigar->Length;
                        pDR->bMateClip = true;
                    }
                }

                //Update the mapping info
                pDR->iMatQueryMapPos = iQueryMapPos;
                pDR->iMatQueryMapLength = strMapSeq.length();
                if(pDR->bMateClip)
                {
                    pDR->iMatQueryClipPos = iQueryClipPos;
                    pDR->iMatQueryClipLength = strClipSeq.length();
                    pDR->strMatClipSeq = strClipSeq;
                }
                //<--
            }
        }
        else
        {
            //If it is new
            stDR.Clear();
            stDR.iReadsPos = al.MatePosition;
            stDR.iMatePos = al.Position;
            stDR.strMateAlignSeq = al.AlignedBases;
            stDR.iInsertSize = iRealInsertSize; //abs(iRealInsertSize);


            for(std::vector<CigarOp>::iterator itrCigar = al.CigarData.begin();
                itrCigar != al.CigarData.end(); itrCigar++)
            {
                if(itrCigar->Type == 'M')
                {
                    //if(strMapSeq == "")
                    //    strMapSeq += al.QueryBases.substr(iOffSet, itr->Length);
                    if(iQueryMapPos == -1)
                        iQueryMapPos = iOffSet;
                    iOffSet += itrCigar->Length;
                }

                if(itrCigar->Type == 'S')
                {
                    if(strClipSeq == "")
                        strClipSeq += al.QueryBases.substr(iOffSet, itrCigar->Length);
                    if(iQueryClipPos == -1)
                        iQueryClipPos = iOffSet;
                    iOffSet += itrCigar->Length;
                    stDR.bMateClip = true;
                }
            }

            //Update the mapping info
            stDR.iMatQueryMapPos = iQueryMapPos;
            stDR.iMatQueryMapLength = strMapSeq.length();
            if(stDR.bMateClip)
            {
                stDR.iMatQueryClipPos = iQueryClipPos;
                stDR.iMatQueryClipLength = strClipSeq.length();
                stDR.strMatClipSeq = strClipSeq;
            }
            //<--

//                //Check if it contains clip part -->
//                for(std::vector<CigarOp>::iterator itrCigar = al.CigarData.begin();
//                    itrCigar != al.CigarData.end(); itrCigar++)
//                {
//                    if(itrCigar->Type == 'S')
//                    {
//                        stDR.bMateClip = true;
//                        break;
//                    }
//                }
//                //<--

            bool bFindChromDR = false;
            for(vector<St_ChromDiscordantReads>::iterator itr = vChromDR.begin();
                itr != vChromDR.end(); itr++)
            {
                if(itr->strChrom == strChrom)
                {
                    itr->vDR.push_back(stDR);
                    bFindChromDR = true;
                    break;
                }
            }
            if(!bFindChromDR)
            {
                stChromDR.Clear();
                stChromDR.strChrom = strChrom;
                stChromDR.vDR.push_back(stDR);
                vChromDR.push_back(stChromDR);
            }
            //vDiscdRreads.push_back(stDR);
        }
    }
}

void ClsParseBam::PrintTestingInfoDR(vector<St_ChromDiscordantReads>& vChromDR)
{
    //output the statistic of what we get
    if(m_pOfs == NULL)
    {
        cout << "vChromDR Size: " << to_string(vChromDR.size()) << endl;
        cout << "====== Details ======" << endl;
        for(vector<St_ChromDiscordantReads>::iterator itr = vChromDR.begin();
            itr != vChromDR.end(); itr++)
        {
            cout << itr->strChrom << ": " <<  to_string(itr->vDR.size()) << endl;
        }
        cout << "=====================" << endl;

        //cout << "Discordant Reads Size: " << to_string(vDiscdRreads.size()) << endl;
    }
    else
    {
        (*m_pOfs) << "vChromDR Size: " << to_string(vChromDR.size()) << endl;
        (*m_pOfs) << "====== Details ======" << endl;
        for(vector<St_ChromDiscordantReads>::iterator itr = vChromDR.begin();
            itr != vChromDR.end(); itr++)
        {
            (*m_pOfs) << itr->strChrom << ": " <<  to_string(itr->vDR.size()) << endl;
        }
        (*m_pOfs) << "=====================" << endl;

        //(*m_pOfs) << "Discordant Reads Size: " << to_string(vDiscdRreads.size()) << endl;
    }
}
