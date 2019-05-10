#include "clsdebug.h"
#include <algorithm>

ClsDebug::ClsDebug()
{
    m_vSvSupportGroup.clear();
    m_strBamFile = "";
    m_strChromName = "";
    m_pOfs = NULL;
}

ClsDebug::~ClsDebug()
{
    m_vSvSupportGroup.clear();
}

void ClsDebug::SetChromName(string strChromName)
{
    this->m_strChromName = strChromName;
}

void ClsDebug::SetOfstream(ofstream& ofs)
{
    this->m_pOfs = &ofs;
}


bool sort_SV_len_func(St_SV stSv1, St_SV stSv2)
{
    if((stSv1.iEnd - stSv1.iPos) < (stSv2.iEnd - stSv2.iPos))
        return true;
    else
        return false;
}

void ClsDebug::GetStdSvDelRelatedGroup(vector<St_DRGroup>& vDRGroup, vector<St_SCRGroup>& vSCRGroup,
                                       vector<St_SV>& vSvDEL, string strChrom)
{
    for(vector<St_DRGroup>::iterator itr = vDRGroup.begin(); itr != vDRGroup.end(); itr++)
    {
        (*m_pOfs) << "( " << to_string(itr->iStart) << ",\t" << to_string(itr->iEnd) << ")"
             << " --- " << to_string(itr->vDR.size())
             << " --- " << to_string(itr->vDRRangedSCR.size())
             << " --- " << to_string(itr->iEnd - itr->iStart)
             << " >>> Depth Info: " << FloatToStr(itr->stGF.stDCRangeCore.fStatDepth /*fPartialDepth*/) << " -- Len: " << to_string(itr->stGF.stDCRangeCore.iStatDepthLen)
             << endl;
    }

    (*m_pOfs) <<"---------"<< endl << endl;

    (*m_pOfs) << "-----vSCRGroup Detail-----" << endl << endl;
    for(vector<St_SCRGroup>::iterator itr = vSCRGroup.begin(); itr != vSCRGroup.end(); itr++)
    {
        (*m_pOfs) << "( " << to_string(itr->iStart) << ",\t" << to_string(itr->iEnd) << ")" << " --- "
             << to_string(itr->vStartSCR.size() + itr->vEndSCR.size())
             << " --- " << to_string(itr->iEnd - itr->iStart)
             << endl;
    }

    (*m_pOfs) << endl << " ============== SV =================" << endl;
    //Check the deletion length of Standard SV
    sort(vSvDEL.begin(), vSvDEL.end(), sort_SV_len_func);
    for(vector<St_SV>::iterator itr = vSvDEL.begin(); itr != vSvDEL.end(); itr++)
    {
        (*m_pOfs) << "<" << to_string(itr->iPos) << ", " << to_string(itr->iEnd) << ">"
             << "\t ... SV_Len: " << to_string(itr->iEnd - itr->iPos) << endl;
    }
    (*m_pOfs) << "============================" << endl;

    //Check how many hit with SV
    m_vSvSupportGroup.clear();
    St_SvSupportGroup stSvSupportGroup;
    int iHitNum = 0;
    for(vector<St_SV>::iterator itr = vSvDEL.begin(); itr != vSvDEL.end(); itr++)
    {
        stSvSupportGroup.Clear();
        stSvSupportGroup.stSv = *itr;

        (*m_pOfs) << "<" << to_string(itr->iPos) << ", " << to_string(itr->iEnd) << ">"
             << " ... SV Len: " << to_string(itr->iEnd - itr->iPos) << endl;
        unsigned int iStartDiff = 500;
        int iRealStartDiff = -1;
        int iEndDiff = -1;
        int iRealEndDiff = -1;
        bool bFind = false;
        int iStart = -1;
        int iEnd = -1;
        int iLeftBoundary = -1;
        int iRightBoundary = -1;
        int iDRNum = 0;
        int iSCRNum = 0;

        //Try to find from discordant reads group
        for(vector<St_DRGroup>::iterator subItr = vDRGroup.begin(); subItr != vDRGroup.end(); subItr++)
        {
            //if(itr->iPos >= subItr->iStart-50 && itr->iPos <= subItr->iEnd+50)
            if(itr->iPos >= subItr->iLeftBoundary && itr->iPos <= subItr->iRightBoundary)
            {
                if(abs(int(itr->iPos - subItr->iStart)) < iStartDiff ||
                   abs(int(itr->iEnd - subItr->iEnd) < iStartDiff))
                {
                    iStart = subItr->iStart;
                    iEnd = subItr->iEnd;
                    iLeftBoundary = subItr->iLeftBoundary;
                    iRightBoundary = subItr->iRightBoundary;

                    iStartDiff = itr->iPos - subItr->iStart;
                    iRealStartDiff = itr->iPos - subItr->iLeftBoundary;
                    //if(iEndDiff > abs(int(itr->iEnd - subItr->iEnd)))
                    //{
                    iEndDiff = subItr->iEnd - itr->iEnd;
                    iRealEndDiff = subItr->iRightBoundary - itr->iEnd;
                    //}
                    bFind = true;
                    stSvSupportGroup.vDRGroup.push_back(*subItr);

                    iDRNum = subItr->vDR.size();
                    iSCRNum = subItr->vDRRangedSCR.size();
                }
            }
        }

        //try to find from clip reads group        
        for(vector<St_SCRGroup>::iterator subItr = vSCRGroup.begin(); subItr != vSCRGroup.end(); subItr++)
        {
            if(itr->iPos >= subItr->iStart-50 && itr->iPos <= subItr->iEnd+50)
            {
                if(abs(int(itr->iPos - subItr->iStart)) < iStartDiff)
                {
                    iStart = subItr->iStart;
                    iEnd = subItr->iEnd;
                    iStartDiff = abs(int(itr->iPos - subItr->iStart));
                    iRealStartDiff = itr->iPos - subItr->iStart;
                    if(iEndDiff > abs(int(itr->iEnd - subItr->iEnd)))
                    {
                        iEndDiff = abs(int(itr->iEnd - subItr->iEnd));
                        iRealEndDiff = subItr->iEnd - itr->iEnd;
                    }
                    bFind = true;
                    stSvSupportGroup.vSCRGroup.push_back(*subItr);
                }
            }
        }

        if(bFind)
        {
            iHitNum++;
            (*m_pOfs) << "\tFind --> Start_Diff & End_Diff: ("
                 << to_string(iStartDiff) << ", " << to_string(iEndDiff) << ")";

            (*m_pOfs) << " --- Real_Start_Diff & Read_End_Diff: ("
                 << to_string(iRealStartDiff) << ", " << to_string(iRealEndDiff) << ")" << endl;

            (*m_pOfs) << "\t\t --- Start_End: <" << to_string(iStart) << ", " << to_string(iEnd)
                 << "> --- Read_Start_End(LR Boundary): <"
                 << to_string(iLeftBoundary) << ", " << to_string(iRightBoundary) << ">" << endl;

            (*m_pOfs) << "\t\t --- Target Reads Counts: " << to_string(iDRNum) << " --- " << to_string(iSCRNum) << endl;
        }
        else
            (*m_pOfs) << "\t --- None ---" << endl;

        m_vSvSupportGroup.push_back(stSvSupportGroup);
    }

    (*m_pOfs) << "Total Number of SV: " << to_string(vSvDEL.size()) << endl;
    (*m_pOfs) << "Number of Hit Sv  : " << to_string(iHitNum) << endl;


    //Check the Coverage -->
    vector<St_DRGroup> vSumDRGroup;
    vector<St_SCRGroup> vSumSCRGroup;
    for(vector<St_SvSupportGroup>::iterator itr = m_vSvSupportGroup.begin(); itr != m_vSvSupportGroup.end(); itr++)
    {
        if(!itr->vDRGroup.empty())
        {
            vSumDRGroup.insert(vSumDRGroup.end(), itr->vDRGroup.begin(), itr->vDRGroup.end());
            (*m_pOfs) << FloatToStr(itr->vDRGroup.begin()->stGF.stDCRangeCore.fStatDepth)
                 << "\t --- " << "Len: " << to_string(itr->vDRGroup.begin()->stGF.stDCRangeCore.iStatDepthLen)
                 << " --- Group Length: " << to_string(itr->vDRGroup.begin()->iEnd - itr->vDRGroup.begin()->iStart)
                 << ", " << to_string(itr->vDRGroup.begin()->iRightBoundary - itr->vDRGroup.begin()->iLeftBoundary)
                 << " -- Real Sv Length: " << to_string(itr->stSv.iEnd - itr->stSv.iPos) << endl;

            (*m_pOfs) << "\t" << "<" << to_string(itr->vDRGroup.begin()->iStart) << ", "
                 << to_string(itr->vDRGroup.begin()->iEnd) << ">"
                 << "\t" << "<" << to_string(itr->vDRGroup.begin()->iLeftBoundary - 50)
                 << ", " << to_string(itr->vDRGroup.begin()->iRightBoundary + 50) << ">"
                 << "\t --- " << "Sv Info: <" << to_string(itr->stSv.iPos)
                 << ", " << to_string(itr->stSv.iEnd) << ">" << endl << endl;
        }

        if(!itr->vSCRGroup.empty())
            vSumSCRGroup.insert(vSumSCRGroup.end(), itr->vSCRGroup.begin(), itr->vSCRGroup.end());
    }

    (*m_pOfs) << endl << "vSumDRGroup Size: " << IntToStr(vSumDRGroup.size()) << endl;
//    for(vector<St_DRGroup>::iterator itr = vSumDRGroup.begin(); itr != vSumDRGroup.end(); itr++)
//    {
////        string strCmd = (string)"samtools depth -a -r " +
////                        "11:" + to_string(itr->iStart) + "-" + to_string(itr->iEnd) + " " + m_strBamFile +
////                        " | awk '{sum+=$3} END {if(sum==\"\") {print 0} else {print sum/NR}}'";

////        string strResult = exec(strCmd.c_str());

////        cout << strResult;
//        cout << to_string(itr->stGF.fPartialDepth) << "\t --- " << "Len: " << to_string(itr->stGF.iDepthRangLen) << endl;
//        cout << "\t" << "<" << to_string(itr->iStart) << ", " << to_string(itr->iEnd) << ">" << endl;
//    }

    (*m_pOfs) << endl << "vSumSCRGroup Size: " << IntToStr(vSumSCRGroup.size()) << endl;
    for(vector<St_SCRGroup>::iterator itr = vSumSCRGroup.begin(); itr != vSumSCRGroup.end(); itr++)
    {
        string strCmd = (string)"samtools depth -a -r " +
                        strChrom + ":" + to_string(itr->iStart) + "-" + to_string(itr->iEnd) + " " + m_strBamFile +
                        " | awk '{sum+=$3} END {if(sum==\"\") {print 0} else {print sum/NR}}'";

        string strResult = exec(strCmd.c_str());

        (*m_pOfs) << strResult;
        (*m_pOfs) << "\t" << "<" << to_string(itr->iStart) << ", " << to_string(itr->iEnd) << ">" << endl;
    }
    //-->
}

void ClsDebug::SetStrBamFile(string strV1)
{
    this->m_strBamFile = strV1;
}

void ClsDebug::SetMeanInserSize(int iV1)
{
    this->m_iMeanInsertSize = iV1;
}

void ClsDebug::UpdateGroupBoundary()
{
//    //Let do the real thing first
//    //1: update boundary by discordant reads only -->  No No No: should do it here-->
//    for(vector<St_SvSupportGroup>::iterator itr = m_vSvSupportGroup.begin(); itr != m_vSvSupportGroup.end(); itr++)
//    {
//        if(itr->vDRGroup.empty())
//            continue;

//        St_DRGroup& stDRGroup = *itr->vDRGroup.begin();

//        FilterBadDR(stDRGroup.vDR);

//        //ReSet Start and End position -->
//        bool bStartGood = false;
//        bool bEndGood = false;
//        for(vector<St_DiscordantReads>::iterator itrDR = stDRGroup.vDR.begin(); itrDR != stDRGroup.vDR.end(); itrDR++)
//        {
//            if(stDRGroup.iStart == itrDR->iReadsPos)
//                bStartGood = true;

//            if(stDRGroup.iEnd == itrDR->iMatePos)
//                bEndGood = true;
//        }
//        if(!bStartGood)
//            stDRGroup.iStart = stDRGroup.iLeftBoundary;

//        if(!bEndGood)
//            stDRGroup.iEnd = stDRGroup.iRightBoundary;
//        //<--
//    }
//    //<--

    (*m_pOfs) << endl << "************ UpdateGroupBoundary ************" << endl;
    for(vector<St_SvSupportGroup>::iterator itr = m_vSvSupportGroup.begin(); itr != m_vSvSupportGroup.end(); itr++)
    {
        (*m_pOfs) << "DR_Group_Num: " << to_string(itr->vDRGroup.size()) << endl;
        if(itr->vDRGroup.empty())
            continue;

        St_DRGroup& stDRGroup = *itr->vDRGroup.begin();

        //For Discordant reads group --> based on the softclip reads dropped into the group
        int iMinPosTailClip = -22222;
        int iMinPosTailClipMat = -22222;
        int iMinPosTailClipMatType = -1; // 1: first mat, 2: second mat

        int iMaxPosTailClip = -22222;
        int iMaxPosTailClipMat = -22222;
        int iMaxPosTailClipMatType = -1;

        int iMinPosHeadClip = -22222;
        int iMinPosHeadClipMat = -22222;
        int iMinPosHeadClipMatType = -1;

        int iMaxPosHeadClip = -22222;
        int iMaxPosHeadClipMat = -22222;
        int iMaxPosHeadClipMatType = -1;

        for(vector<St_BorderSCR>::iterator itr = stDRGroup.vDRRangedSCR.begin();
            itr != stDRGroup.vDRRangedSCR.end(); itr++)
        {
            //For tail clipped (right part clipped)
            if(itr->enClipPart == cpRight)
            {
                if(iMinPosTailClip == -22222)
                {
                    iMinPosTailClip = itr->iClipPos;
                    iMinPosTailClipMat = itr->iMatPos;
                    if(itr->bFirstMate)
                        iMinPosTailClipMatType = 2;
                    if(itr->bSecondMate)
                        iMinPosTailClipMatType = 1;
                }
                else if(iMinPosTailClip > itr->iClipPos)
                {
                    iMinPosTailClip = itr->iClipPos;
                    iMinPosTailClipMat = itr->iMatPos;
                    if(itr->bFirstMate)
                        iMinPosTailClipMatType = 2;
                    if(itr->bSecondMate)
                        iMinPosTailClipMatType = 1;
                }

                if(iMaxPosTailClip == -22222)
                {
                    iMaxPosTailClip = itr->iClipPos;
                    iMaxPosTailClipMat = itr->iMatPos;
                    if(itr->bFirstMate)
                        iMaxPosTailClipMatType = 2;
                    if(itr->bSecondMate)
                        iMaxPosTailClipMatType = 1;
                }
                else if(iMinPosTailClip < itr->iClipPos)
                {
                    iMaxPosTailClip = itr->iClipPos;
                    iMaxPosTailClipMat = itr->iMatPos;
                    if(itr->bFirstMate)
                        iMaxPosTailClipMatType = 2;
                    if(itr->bSecondMate)
                        iMaxPosTailClipMatType = 1;
                }
            }

            //For head clipped (right part clipped)
            if(itr->enClipPart == cpLeft)
            {
                if(iMinPosHeadClip == -22222)
                {
                    iMinPosHeadClip = itr->iClipPos;
                    iMinPosHeadClipMat = itr->iMatPos;
                    if(itr->bFirstMate)
                        iMinPosHeadClipMatType = 2; //since it is the mat, so it should be opposit than current type
                    if(itr->bSecondMate)
                        iMinPosHeadClipMatType = 1;
                }
                else if(iMinPosHeadClip > itr->iClipPos)
                {
                    iMinPosHeadClip = itr->iClipPos;
                    iMinPosHeadClipMat = itr->iMatPos;
                    if(itr->bFirstMate)
                        iMinPosHeadClipMatType = 2; //since it is the mat, so it should be opposit than current type
                    if(itr->bSecondMate)
                        iMinPosHeadClipMatType = 1;
                }

                if(iMaxPosHeadClip == -22222)
                {
                    iMaxPosHeadClip = itr->iClipPos;
                    iMaxPosHeadClipMat = itr->iMatPos;
                    if(itr->bFirstMate)
                        iMaxPosHeadClipMatType = 2; //since it is the mat, so it should be opposit than current type
                    if(itr->bSecondMate)
                        iMaxPosHeadClipMatType = 1;
                }
                else if(iMaxPosHeadClip < itr->iClipPos)
                {
                    iMaxPosHeadClip = itr->iClipPos;
                    iMaxPosHeadClipMat = itr->iMatPos;
                    if(itr->bFirstMate)
                        iMaxPosHeadClipMatType = 2; //since it is the mat, so it should be opposit than current type
                    if(itr->bSecondMate)
                        iMaxPosHeadClipMatType = 1;
                }
            }
        }

        //Output something
        int iStartDiff = itr->stSv.iPos - (stDRGroup.iStart + stDRGroup.iStartAlignLen);
        if(stDRGroup.bStartClipAdjusted)
            iStartDiff = itr->stSv.iPos - stDRGroup.iStart;

        int iEndDiff = stDRGroup.iEnd - itr->stSv.iEnd;        

        int iRealStartDiff = itr->stSv.iPos - (stDRGroup.iLeftBoundary + stDRGroup.iLeftBoundaryAlignLen);
        int iRealEndDiff = stDRGroup.iRightBoundary - itr->stSv.iEnd;

        int iSCRStartDiff = -22222;
        if(iMaxPosTailClip >= 0)
            iSCRStartDiff = itr->stSv.iPos - iMaxPosTailClip;
        int iSCREndDiff = -22222;
        if(iMinPosHeadClip >= 0)
            iSCREndDiff = iMinPosHeadClip - itr->stSv.iEnd;

        int iSCRRealStartDiff = -22222;
        if(iMinPosTailClip >= 0)
            iSCRRealStartDiff = itr->stSv.iPos - iMinPosTailClip;
        int iSCRRealEndDiff = -22222;
        if(iMaxPosHeadClip >= 0)
            iSCRRealEndDiff = iMaxPosHeadClip - itr->stSv.iEnd;

        (*m_pOfs) << "DR_Group --> SV Len: " << to_string(itr->stSv.GetLen()) << " --- Diff: ("
             << to_string(iStartDiff) << ", " << to_string(iEndDiff) << ")";

        (*m_pOfs) << " --- Real_Diff: ("
             << to_string(iRealStartDiff) << ", " << to_string(iRealEndDiff) << ")";

        (*m_pOfs) << " --- SCR_Diff: ("
             << to_string(iSCRStartDiff) << ", " << to_string(iSCREndDiff) << ")";

        (*m_pOfs) << " --- Real_SCR_Diff: ("
             << to_string(iSCRRealStartDiff) << ", " << to_string(iSCRRealEndDiff) << ")" << endl;

        //Cout Mate Diff
        int iMatDiff11 = -222222;
        if(iMaxPosTailClipMatType == 1)
            iMatDiff11 = stDRGroup.iLeftBoundary - iMaxPosTailClipMat;
        else if(iMaxPosTailClipMatType == 2)
            iMatDiff11 = stDRGroup.iRightBoundary - iMaxPosTailClipMat;

        int iMatDiff12 = -222222;
        if(iMinPosHeadClipMatType == 1)
            iMatDiff12 = stDRGroup.iLeftBoundary - iMinPosHeadClipMat;
        else if(iMinPosHeadClipMatType == 2)
            iMatDiff12 = stDRGroup.iRightBoundary - iMinPosHeadClipMat;

        int iMatDiff21 = -222222;
        if(iMinPosTailClipMatType == 1)
            iMatDiff21 = stDRGroup.iLeftBoundary - iMinPosTailClipMat;
        else if(iMinPosTailClipMatType == 2)
            iMatDiff21 = stDRGroup.iRightBoundary - iMinPosTailClipMat;

        int iMatDiff22 = -222222;
        if(iMaxPosHeadClipMatType == 1)
            iMatDiff22 = stDRGroup.iLeftBoundary - iMaxPosHeadClipMat;
        else if(iMaxPosHeadClipMatType == 2)
            iMatDiff22 = stDRGroup.iRightBoundary - iMaxPosHeadClipMat;

        (*m_pOfs) << "\t"
             << "(" << to_string(iMaxPosTailClipMatType) << ", " << to_string(iMinPosHeadClipMatType) << ")" << ", "
             << "<" << to_string(iMatDiff11) << ", " << to_string(iMatDiff12) << ">" << " ---- "
             << "(" << to_string(iMinPosTailClipMatType) << ", " << to_string(iMaxPosHeadClipMatType) << ")" << ", "
             << "<" << to_string(iMatDiff21) << ", " << to_string(iMatDiff22) << ">"
             << " - " << "StartAlignLen: " << to_string(stDRGroup.iStartAlignLen)
             << " - " << "LeftBoundAlignLen: " << to_string(stDRGroup.iLeftBoundaryAlignLen)
             << " - SV: <" << to_string(itr->stSv.iPos) << ", " << to_string(itr->stSv.iEnd) << ">"
             << ", " << to_string(itr->stSv.iEnd - itr->stSv.iPos) << endl;

        //output all the supported reads for the target SV
//        if((itr->stSv.iPos == 22470785 && itr->stSv.iEnd == 22472850) ||
//           (itr->stSv.iPos == 41818199 && itr->stSv.iEnd == 41819025))
        {
            //For DR group
            (*m_pOfs) << ">>>>>>>>>>>>> DR Group: " << to_string(stDRGroup.vDR.size()) << endl;
            for(vector<St_DiscordantReads>::iterator itrDR = stDRGroup.vDR.begin(); itrDR != stDRGroup.vDR.end(); itrDR++)
            {
                (*m_pOfs) << "< " << to_string(itrDR->iReadsPos) << ", " << to_string(itrDR->iMatePos) << ">"
                     << " --- " << to_string(itr->stSv.iPos - (itrDR->iReadsPos + (int)itrDR->strReadsAlignSeq.length()))
                                << ", " << to_string(itrDR->iMatePos - itr->stSv.iEnd);
                (*m_pOfs) << "\t ===> ClipInfo: Start(" << to_string(itrDR->iQueryMapLength) << ", "
                     << (itrDR->iQueryMapLength < 90 ? "T" : "F") << "), "
                     << "End( " << to_string(itrDR->iMatQueryMapLength) << ", "
                     << (itrDR->iMatQueryMapLength < 90 ? "T" : "F") << ")" << endl;
            }

            //For SCR dropped into current group
            (*m_pOfs) << " ---- " << endl << ">>>>>>>>>>>>> SCR Group: " << to_string(stDRGroup.vDRRangedSCR.size()) << endl;
            for(vector<St_BorderSCR>::iterator itrSCR = stDRGroup.vDRRangedSCR.begin();
                itrSCR != stDRGroup.vDRRangedSCR.end(); itrSCR++)
            {
                if(itrSCR->enClipPart == cpRight)
                {
                    (*m_pOfs) << "Tail Cut (Left Bound): < " << to_string(itrSCR->iReadsMapPos) << ">"
                         << " --- " << to_string(itr->stSv.iPos - itrSCR->iClipPos) << endl;

                    if(itrSCR->bFirstMate)                    
                        (*m_pOfs) << "\t Mate_Info: 2 > " << to_string(itrSCR->iMatPos - stDRGroup.iRightBoundary) << endl;
                    else                    
                        (*m_pOfs) << "\t Mate_Info: 1 > " << to_string(stDRGroup.iLeftBoundary - itrSCR->iMatPos) << endl;

                    //cout diff
                    (*m_pOfs) << "\t Diff(compare with START): " << to_string(itrSCR->iClipPos - stDRGroup.iStart) << endl;
                }

                if(itrSCR->enClipPart == cpLeft)
                {
                    (*m_pOfs) << "Head Cut (Right Bound): < " << to_string(itrSCR->iReadsMapPos) << ">"
                         << " --- " << to_string(itrSCR->iClipPos - itr->stSv.iEnd) << endl;

                    if(itrSCR->bFirstMate)                    
                        (*m_pOfs) << "\t Mate_Info: 2 > " << to_string(itrSCR->iMatPos - stDRGroup.iRightBoundary) << endl;
                    else                    
                        (*m_pOfs) << "\t Mate_Info: 1 > " << to_string(stDRGroup.iLeftBoundary - itrSCR->iMatPos) << endl;

                    (*m_pOfs) << "\t Diff(compare with END): " << to_string(itrSCR->iClipPos - stDRGroup.iEnd) << endl;
                }
                (*m_pOfs) << "\t " << "Insert_Size: " << to_string(itrSCR->iInsertSize)
                     << ", " << to_string(itrSCR->iMatPos) << endl;
            }
            (*m_pOfs) << "###################################" << endl;
        }
    }
}
