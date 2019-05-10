#include "clssvdeldetect.h"
#include <algorithm> //for sort
#include <cmath> //for ceil and floor

ClsSvDelDetect::ClsSvDelDetect():m_fMeanInsertSize(-1.), m_fStdDevInsertSize(-1.), m_fAvgDepth(-1.),
                                 m_strChromName(""), m_pOfs(NULL)
{}

ClsSvDelDetect::~ClsSvDelDetect()
{}

void ClsSvDelDetect::SetMeanInsertSize(float fV1)
{
    this->m_fMeanInsertSize = fV1;
}

void ClsSvDelDetect::SetStdDevInsertSize(float fV1)
{
    this->m_fStdDevInsertSize = fV1;
}

void ClsSvDelDetect::SetAvgDepth(float fV1)
{
    this->m_fAvgDepth = fV1;
}

void ClsSvDelDetect::SetChromName(string strChromName)
{
    this->m_strChromName = strChromName;
}

void ClsSvDelDetect::SetOfstream(ofstream& ofs)
{
    this->m_pOfs = &ofs;
}

bool sort_discordantReads_large_small_func(St_DiscordantReads stDR1, St_DiscordantReads stDR2)
{
    if(stDR1.iReadsPos >= stDR2.iReadsPos)
        return true;
    else
        return false;
}

//bool sort_discordantReads_mat_small_large_func(St_DiscordantReads stDR1, St_DiscordantReads stDR2)
//{
//    if(stDR1.iMatePos <= stDR2.iMatePos)
//        return true;
//    else
//        return false;
//}

void ClsSvDelDetect::ClusterDiscrodantReadsDelly(vector<St_DRGroup>& vDRGroup, vector<St_DiscordantReads>& vDiscdRreads,
                                                 vector<St_BorderSCR>& vBorderSCR)
{
    (*m_pOfs) << "m_fAvgDepth: " << FloatToStr(this->m_fAvgDepth) << endl;

    //1: Sort current discordant reads data se
    sort(vDiscdRreads.begin(), vDiscdRreads.end(), sort_discordantReads_large_small_func);
    //cout << "sort finished <--" << endl;

    //2: Make clustering
    if(vDiscdRreads.size() <= 1)
        return;

    //Filter Bad Discordant reads
    //FilterBadDR(vDiscdRreads);

    int iAllowedDiff = 101; // Now we set this diff as Reads Length
    St_DRGroup stDRGroup;
    int iSplitPos = 0;
    for(vector<St_DiscordantReads>::iterator itr = vDiscdRreads.end() - 1; itr >= vDiscdRreads.begin(); itr--)
    {
//        if(itr->iReadsPos <= iSplitPos)
//            continue;

        vector<St_DiscordantReads>::iterator tmpItr = itr;
        stDRGroup.vDR.push_back(*itr);
        stDRGroup.iEnd = itr->iMatePos;
        stDRGroup.iStart = itr->iReadsPos;
        stDRGroup.iLeftBoundary = itr->iReadsPos;
        stDRGroup.iRightBoundary = itr->iMatePos;
        stDRGroup.iStartAlignLen = itr->strReadsAlignSeq.length();
        stDRGroup.iLeftBoundaryAlignLen = itr->strReadsAlignSeq.length();

        for(vector<St_DiscordantReads>::iterator subItr = itr - 1; subItr >= vDiscdRreads.begin(); subItr--)
        {
            if(abs(int(subItr->iReadsPos - itr->iReadsPos)) < iAllowedDiff)
            {
                stDRGroup.vDR.push_back(*subItr);
                tmpItr = subItr;

                //Start should be the largest position, while, End should be the smallest position
                if(stDRGroup.iStart < tmpItr->iReadsPos)
                {
                    stDRGroup.iStart = tmpItr->iReadsPos;
                    stDRGroup.iStartAlignLen = tmpItr->strReadsAlignSeq.length();
                }

                if(stDRGroup.iEnd > tmpItr->iMatePos)
                    stDRGroup.iEnd = tmpItr->iMatePos;

                //For real boundary -->
                if(stDRGroup.iLeftBoundary > tmpItr->iReadsPos)
                {
                    stDRGroup.iLeftBoundary = tmpItr->iReadsPos;
                    stDRGroup.iLeftBoundaryAlignLen = tmpItr->strReadsAlignSeq.length();
                }

                if(stDRGroup.iRightBoundary < tmpItr->iMatePos)
                    stDRGroup.iRightBoundary = tmpItr->iMatePos;
                //<--
            }
            else
                break;
        }
        itr = tmpItr;

        //Add current group into group set
        if(stDRGroup.vDR.size() > 1) // erase isolated group (only contributed by 1 PE reads)
        {
            //we only keep the group contributed by more than one PE discordant reads
            vDRGroup.push_back(stDRGroup);
            iSplitPos = stDRGroup.iRightBoundary;
        }

        //Clean current temporary stDRGroup
        stDRGroup.Clear();
    }

    //Now we have done !!
    (*m_pOfs) << "The Number of DRGroup: " << to_string(vDRGroup.size()) << endl;

    //Combine the overlapped group -->
    //Advantage: 使得结果更加的整洁，组和组之间不存在任何的交集，比较独立
    //Disadvantage: 如果存在异常大的discordant reads, 而且这个discordant reads被归属到了某个组里面，那么这会导致在这个范围内的所有的组都被合并到一起了.
    (*m_pOfs) << endl <<  "---------" << endl;
    vector<St_DRGroup> vNewDRGroup;
    for(vector<St_DRGroup>::iterator itr = vDRGroup.begin(); itr != vDRGroup.end(); itr++)
    {
//        cout << "( " << to_string(itr->iStart) << ",\t" << to_string(itr->iEnd) << ")" << " --- " << to_string(itr->vDR.size()) << endl;

        stDRGroup.Clear();
        stDRGroup = *itr;

        bool bSaved = false;
        for(vector<St_DRGroup>::iterator subItr = itr+1; subItr < vDRGroup.end(); subItr++)
        {
            //dropped in the range
            if(subItr->iLeftBoundary >= stDRGroup.iLeftBoundary && subItr->iLeftBoundary <= stDRGroup.iRightBoundary)
            {
                stDRGroup.vDR.insert(stDRGroup.vDR.end(), subItr->vDR.begin(), subItr->vDR.end());
                //Update start and end
                if(subItr->iStart > stDRGroup.iStart)
                {
                    stDRGroup.iStart = subItr->iStart;
                    stDRGroup.iStartAlignLen = subItr->iStartAlignLen;
                }
                if(subItr->iEnd< stDRGroup.iEnd)
                    stDRGroup.iEnd = subItr->iEnd;

                //update real boundary --> we need to keep the shortest boundary range --> this is incorrect -->
                //Here --> 这个之前的逻辑错误是肯定需要修改的。LeftBoundary should be the real left boundary and RightBounday should be the real right boundary
                if(subItr->iLeftBoundary < stDRGroup.iLeftBoundary)
                {
                    stDRGroup.iLeftBoundary = subItr->iLeftBoundary;
                    stDRGroup.iLeftBoundaryAlignLen = subItr->iLeftBoundaryAlignLen;
                }
                if(subItr->iRightBoundary > stDRGroup.iRightBoundary)
                    stDRGroup.iRightBoundary = subItr->iRightBoundary;

                itr = subItr;

                if(subItr == vDRGroup.end())
                {
                     vNewDRGroup.push_back(stDRGroup);
                     bSaved =  true;
                     break;
                }

            }
            else
            {
                vNewDRGroup.push_back(stDRGroup);
                bSaved =  true;
                break;
            }
        }
//        if(!bSaved)
//            vNewDRGroup.push_back(stDRGroup);
    }

    (*m_pOfs) << "After merge: " << to_string(vNewDRGroup.size()) << endl;
    (*m_pOfs) << "vBorderSCR size: " << to_string(vBorderSCR.size()) << endl;

    //Filter the discordant reads which do not contain any split reads --> Go
    if(vBorderSCR.empty())
    {
        (*m_pOfs) << "vBorderSCR is empty" << endl;
        vDRGroup.clear();
        vDRGroup.insert(vDRGroup.end(), vNewDRGroup.begin(), vNewDRGroup.end());
        return;
    }

    //if vBorderSCR is not empty
    vector<St_BorderSCR>::iterator itrTmpScr = vBorderSCR.end() - 1;
    for(vector<St_DRGroup>::iterator itr = vNewDRGroup.end() - 1; itr >= vNewDRGroup.begin(); itr--) //iterator Group from large to small
    {
        for(vector<St_BorderSCR>::iterator subItr = itrTmpScr; subItr >= vBorderSCR.begin(); subItr--) //iterator softclip reads from large position to small as well
        {
            if(itr->iEnd >= subItr->iClipPos)
            {
                //if(itr->iStart + itr->iStartAlignLen <= subItr->iClipPos)
                if(itr->iStart <= subItr->iClipPos)
                {
                    //Keep it
                    itrTmpScr = subItr;
                    itr->vDRRangedSCR.push_back(*subItr);
                    subItr->bDrippedIntoDRGroup = true;
                    //Put all of SCR with in the range of current group into the data structure
                    for(vector<St_BorderSCR>::iterator subItr = itrTmpScr - 1; subItr >= vBorderSCR.begin(); subItr--)
                    {
                        if(itr->iLeftBoundary <= subItr->iClipPos)
                        {
                            itrTmpScr = subItr;
                            itr->vDRRangedSCR.push_back(*subItr);
                            subItr->bDrippedIntoDRGroup = true;
                        }
                        else
                            break;
                    }
                }
                else
                {
                    vNewDRGroup.erase(itr);
                    itrTmpScr = subItr;
                }
                break;
            }
        }

//        //This logic is so confused --> i rewrite this part:
//        //Notice: all the boundary should come from iLeftBoundary and iRightBoundary -->
//        for(vector<St_BorderSCR>::iterator subItr = itrTmpScr; subItr >= vBorderSCR.begin(); subItr--) //iterator softclip reads from large position to small as well
//        {
//            if(itr->iRightBoundary >= subItr->iClipPos)
//            {
//                if(itr->iLeftBoundary <= subItr->iClipPos)
//                {
//                    //Keep it
//                    itrTmpScr = subItr;
//                    itr->vDRRangedSCR.push_back(*subItr);
//                    subItr->bDrippedIntoDRGroup = true;
//                    //Put all of SCR with in the range of current group into the data structure
//                    for(vector<St_BorderSCR>::iterator subItr = itrTmpScr - 1; subItr >= vBorderSCR.begin(); subItr--)
//                    {
//                        if(itr->iLeftBoundary <= subItr->iClipPos)
//                        {
//                            itrTmpScr = subItr;
//                            itr->vDRRangedSCR.push_back(*subItr);
//                            subItr->bDrippedIntoDRGroup = true;
//                        }
//                        else
//                            break;
//                    }
//                }
//                else
//                {
//                    vNewDRGroup.erase(itr);
//                    itrTmpScr = subItr;
//                }
//                break;
//            }
//        }
    }
    (*m_pOfs) << "After clip reads filter: " << to_string(vNewDRGroup.size()) << endl;
    //<-- Now Group by discordant reads has been finished
    vDRGroup.clear();
    vDRGroup.insert(vDRGroup.end(), vNewDRGroup.begin(), vNewDRGroup.end());


    //Update the boundary of vDRGroup ******************************
    for(vector<St_DRGroup>::iterator itr = vDRGroup.begin(); itr != vDRGroup.end(); itr++)
    {
        FilterBadDR(itr->vDR);
        UpdateStartEndByDR(*itr); // Based on repeat DR reads
    }
    //<--

    //Update the boundary of vDRGroup by Softclip reads ******************************
    for(vector<St_DRGroup>::iterator itr = vDRGroup.end() - 1; itr >= vDRGroup.begin(); itr--)
    {
        //--> For debug
//        if(itr->iStart == 5106720 && itr->iEnd == 5107210)
        {
            UpdateBoundBySCR(*itr);
        }
//        else
//        {
//            vDRGroup.erase(itr);
//        }
        //<--
    }
    //<--
}

void ClsSvDelDetect::UpdateStartEndByDR(St_DRGroup& stDRGroup)
{
    if(stDRGroup.iLeftBoundary == 25564025 || stDRGroup.iRightBoundary == 25565420)
    {
        (*m_pOfs) << "I am the target" << endl;
    }

    //ReSet Start and End position -->
    bool bStartGood = false;
    bool bEndGood = false;

    //我们collect出现两次的　start point 和　ending point　然后以那个为参考值进行start和end的更新
    int iStartRepeatPos = -1;
    int iRepeatStartPosAlignLen = -1;
    bool bStartClip = false;

    int iEndRepeatPos = -1;

    map<int, int> mpStartPos;
    map<int, int> mpEndPos;

    for(vector<St_DiscordantReads>::iterator itrDR = stDRGroup.vDR.begin(); itrDR != stDRGroup.vDR.end(); itrDR++)
    {
        //--> This is the old irelavant logic
        if(stDRGroup.iStart == itrDR->iReadsPos)
            bStartGood = true;

        if(stDRGroup.iEnd == itrDR->iMatePos)
            bEndGood = true;
        //<--

        //-->
        //For start position
        if(mpStartPos.find(itrDR->iReadsPos + itrDR->iQueryMapLength) != mpStartPos.end()) // find it
        {
            if(itrDR->iReadsPos != mpStartPos[itrDR->iReadsPos + itrDR->iQueryMapLength])
                bStartClip = true;

            if(iStartRepeatPos < 0) // init value
            {
                iStartRepeatPos = itrDR->iReadsPos;
                iRepeatStartPosAlignLen = itrDR->iQueryMapLength;
            }
            else
            {
                if(itrDR->iReadsPos > iStartRepeatPos)
                {
                    iStartRepeatPos = itrDR->iReadsPos;
                    iRepeatStartPosAlignLen = itrDR->iQueryMapLength;
                }
            }
        }
        else
            mpStartPos[itrDR->iReadsPos + itrDR->iQueryMapLength] = itrDR->iReadsPos;

        //For end position
        if(mpEndPos.find(itrDR->iMatePos) != mpEndPos.end()) // find it
        {
            if(iEndRepeatPos < 0) // init value
                iEndRepeatPos = itrDR->iMatePos;
            else
            {
                if(itrDR->iMatePos < iEndRepeatPos)
                    iEndRepeatPos = itrDR->iMatePos;
            }
        }
        else
            mpEndPos[itrDR->iMatePos] = 1;
        //<--
    }
    if(!bStartGood)
        stDRGroup.iStart = stDRGroup.iLeftBoundary;

    if(!bEndGood)
        stDRGroup.iEnd = stDRGroup.iRightBoundary;

    //Update Boundary -->
    if(iStartRepeatPos > 0 && bStartClip)
    {
        stDRGroup.iStart = iStartRepeatPos;
        stDRGroup.iStartAlignLen = iRepeatStartPosAlignLen;
    }

    if(iEndRepeatPos > 0)
        stDRGroup.iEnd = iEndRepeatPos;

//    itr->iRepeatStartPos = iStartRepeatPos;
//    itr->iRepeatStartPosAlignLen = iRepeatStartPosAlignLen;
//    itr->iRepeatEndPos = iEndRepeatPos;
    //<--
}

bool sort_scrbound_count_large_to_small(St_SCRBound st1, St_SCRBound st2)
{
    if(st1.iCount > st2.iCount)
        return true;
    else
        return false;
}

bool sort_scrbound_iDiff_small_to_large(St_SCRBound st1, St_SCRBound st2)
{
    if(st1.iDiff < st1.iDiff)
        return true;
    else
        return false;
}

bool sort_scrbound_iDiff_closestmat_small_to_large(St_SCRBound st1, St_SCRBound st2)
{
    if(abs(st1.iClosestMat) < abs(st1.iClosestMat))
        return true;
    else
        return false;
}

void ClsSvDelDetect::UpdateBoundBySCR(St_DRGroup& stDRGroup)
{
    if(stDRGroup.iLeftBoundary == 25564025 || stDRGroup.iRightBoundary == 25565420)
    {
        (*m_pOfs) << "It's our target" <<endl;
    }

    //Get the maximum  appeared number
    vector<St_SCRBound> vLeftBoundSCR;
    vector<St_SCRBound> vRightBoundSCR;

    int iMaxDiff = 101;
    int iMaxClosestOffSet = 101;
    int iMaxDiffWeaker = 101 * 2.5;

    St_SCRBound stSCRBound;

    for(vector<St_BorderSCR>::iterator itr = stDRGroup.vDRRangedSCR.begin();
        itr != stDRGroup.vDRRangedSCR.end(); itr++)
    {
        //Escapte the reads with bad insert size -->
        if(abs(itr->iInsertSize) > (stDRGroup.iEnd - stDRGroup.iStart) * 10)
            continue;
        //<--

        //(1) this should contribute left boundary
        if(itr->enClipPart == cpRight)
        {
            bool bFind = false;
            for(vector<St_SCRBound>::iterator itrTmp = vLeftBoundSCR.begin(); itrTmp != vLeftBoundSCR.end(); itrTmp++)
            {
                if(abs(itr->iClipPos - itrTmp->iClipPos) <= 1)
                {
                    itrTmp->vSCR.push_back(*itr);
                    itrTmp->iCount++;

                    //Calc Clost mat
                    int iClosestMat = 111111;
                    if(itr->bFirstMate)
                    {
                        iClosestMat = itr->iMatPos - stDRGroup.iRightBoundary; // this is the second mate
                    }
                    else
                    {
                        iClosestMat = stDRGroup.iLeftBoundary - itr->iMatPos; // this is the frist mate
                    }

                    if(//iClosestMat > 0 &&
                       abs(itrTmp->iClosestMat) > abs(iClosestMat))
                    {
                        itrTmp->iClosestMat = iClosestMat;
                        itrTmp->iClipPos = itr->iClipPos;
                    }

                    bFind = true;
                    break;
                }
            }
            if(!bFind) //This is a new one
            {
                stSCRBound.Clear();
                stSCRBound.vSCR.push_back(*itr);
                stSCRBound.iClipPos = itr->iClipPos;
                stSCRBound.iCount++;
                stSCRBound.iDiff = abs(itr->iClipPos - (stDRGroup.iStart + stDRGroup.iStartAlignLen));

                //Calc Clost mat
                if(itr->bFirstMate)
                {
                    stSCRBound.iClosestMat = itr->iMatPos - stDRGroup.iRightBoundary;
                }
                else
                {
                    stSCRBound.iClosestMat = stDRGroup.iLeftBoundary - itr->iMatPos;
                }
//                if(stSCRBound.iClosestMat < 0)
//                    stSCRBound.iClosestMat = 100000;

                vLeftBoundSCR.push_back(stSCRBound);
            }
        }

        //(2) this shuld be the right boundary
        if(itr->enClipPart == cpLeft)
        {
            bool bFind = false;
            for(vector<St_SCRBound>::iterator itrTmp = vRightBoundSCR.begin(); itrTmp != vRightBoundSCR.end(); itrTmp++)
            {
                if(abs(itr->iClipPos - itrTmp->iClipPos) <= 1)
                {
                    itrTmp->vSCR.push_back(*itr);
                    itrTmp->iCount++;

                    //Calc Clost mat
                    int iClosestMat = 100000;
                    if(itr->bFirstMate)
                    {
                        iClosestMat = itr->iMatPos - stDRGroup.iRightBoundary; // this is the second mate
                    }
                    else
                    {
                        iClosestMat = stDRGroup.iLeftBoundary - itr->iMatPos; // this is the frist mate
                    }

                    if(abs(itrTmp->iClosestMat) > abs(iClosestMat))
                    {
                        itrTmp->iClosestMat = iClosestMat;
                        itrTmp->iClipPos = itr->iClipPos;
                    }

                    bFind = true;
                    break;
                }
            }
            if(!bFind) //This is a new one
            {
                stSCRBound.Clear();
                stSCRBound.vSCR.push_back(*itr);
                stSCRBound.iClipPos = itr->iClipPos;
                stSCRBound.iCount++;                
                stSCRBound.iDiff = abs(itr->iClipPos - stDRGroup.iEnd);

                //Calc Clost mat
                if(itr->bFirstMate)
                {
                    stSCRBound.iClosestMat = itr->iMatPos - stDRGroup.iRightBoundary;
                }
                else
                {
                    stSCRBound.iClosestMat = stDRGroup.iLeftBoundary - itr->iMatPos;
                }
//                if(stSCRBound.iClosestMat < 0)
//                    stSCRBound.iClosestMat = 100000;

                vRightBoundSCR.push_back(stSCRBound);
            }
        }
    }

    //--> Check if use repeat start or end position
    bool bLeftUpdated = false;
    bool bRightUpdated = false;
    UpdateStartEndBySameClipOrMate(stDRGroup, vLeftBoundSCR, vRightBoundSCR, bLeftUpdated, bRightUpdated);
    //<--

    //For left boundary SCR
    if(!vLeftBoundSCR.empty() && !bLeftUpdated)
    {
        sort(vLeftBoundSCR.begin(), vLeftBoundSCR.end(), sort_scrbound_count_large_to_small);
        vector<St_SCRBound> vTargetLeftBoundSCR;
        int iMaxCount = vLeftBoundSCR.begin()->iCount;
        if(iMaxCount > 1)
        {
            for(vector<St_SCRBound>::iterator itr = vLeftBoundSCR.begin(); itr != vLeftBoundSCR.end(); itr++)
            {
                if(itr->iCount >= iMaxCount)
                {
                    vTargetLeftBoundSCR.push_back(*itr);
                }
                else
                    break;
            }
        }
        if(!vTargetLeftBoundSCR.empty())
        {
            //check Small diff
            sort(vTargetLeftBoundSCR.begin(), vTargetLeftBoundSCR.end(), sort_scrbound_iDiff_small_to_large);
            if(vTargetLeftBoundSCR.begin()->iDiff < iMaxDiff && vTargetLeftBoundSCR.begin()->iClipPos != stDRGroup.iStart)
            {
                stDRGroup.iStart = vTargetLeftBoundSCR.begin()->iClipPos;
                stDRGroup.bStartClipAdjusted = true;
            }
            else
            {
                sort(vTargetLeftBoundSCR.begin(), vTargetLeftBoundSCR.end(), sort_scrbound_iDiff_closestmat_small_to_large);
                if( abs(vTargetLeftBoundSCR.begin()->iClosestMat) < iMaxClosestOffSet &&
                   vTargetLeftBoundSCR.begin()->iDiff < iMaxDiffWeaker &&
                   vTargetLeftBoundSCR.begin()->iClipPos != stDRGroup.iStart)
                {
                    stDRGroup.iStart = vTargetLeftBoundSCR.begin()->iClipPos;
                    stDRGroup.bStartClipAdjusted = true;
                }
            }
        }
    }

    //For right boundary SCR
    if(!vRightBoundSCR.empty() && !bRightUpdated)
    {
        sort(vRightBoundSCR.begin(), vRightBoundSCR.end(), sort_scrbound_count_large_to_small);
        vector<St_SCRBound> vTargetRightBoundSCR;
        int iMaxCount = vRightBoundSCR.begin()->iCount;
        if(iMaxCount > 1)
        {
            for(vector<St_SCRBound>::iterator itr = vRightBoundSCR.begin(); itr != vRightBoundSCR.end(); itr++)
            {
                if(itr->iCount >= iMaxCount)
                {
                    vTargetRightBoundSCR.push_back(*itr);
                }
                else
                    break;
            }
        }
        if(!vTargetRightBoundSCR.empty())
        {
            //check Small diff
            sort(vTargetRightBoundSCR.begin(), vTargetRightBoundSCR.end(), sort_scrbound_iDiff_small_to_large);
            if(vTargetRightBoundSCR.begin()->iDiff < iMaxDiff && vTargetRightBoundSCR.begin()->iClipPos != stDRGroup.iEnd)
            {
                stDRGroup.iEnd = vTargetRightBoundSCR.begin()->iClipPos;
                stDRGroup.bEndClipAdjusted = true;
            }
            else
            {
                sort(vTargetRightBoundSCR.begin(), vTargetRightBoundSCR.end(), sort_scrbound_iDiff_closestmat_small_to_large);
                if(abs(vTargetRightBoundSCR.begin()->iClosestMat) < iMaxClosestOffSet &&
                   vTargetRightBoundSCR.begin()->iDiff < iMaxDiffWeaker &&
                   vTargetRightBoundSCR.begin()->iClipPos != stDRGroup.iEnd)
                {
                    stDRGroup.iEnd = vTargetRightBoundSCR.begin()->iClipPos;
                    stDRGroup.bEndClipAdjusted = true;
                }
            }
        }
    }
}

float ClsSvDelDetect::GetMinValidIS()
{
    return m_fMeanInsertSize - 3 * m_fStdDevInsertSize;
}

float ClsSvDelDetect::GetMaxValidIS()
{
    return m_fMeanInsertSize + 3 * m_fStdDevInsertSize;
}

void ClsSvDelDetect::UpdateStartEndBySameClipOrMate(St_DRGroup& stDRGroup,
                                                    vector<St_SCRBound>& vLeftBoundSCR,
                                                    vector<St_SCRBound>& vRightBoundSCR,
                                                    bool& bLeftUpdated, bool& bRightUpdated)
{
    if(stDRGroup.iLeftBoundary == 12943496 || stDRGroup.iRightBoundary == 12944121)
        (*m_pOfs) << "I am target!" << endl;

    //For left part
    int iMaxDiff = 100 * 3; // triple times of  reads length
    int iValidClipPos = -1;

    int iMatDiffMax = 100 * 2.5;
    int iClipDiffStart = 100 * 2.5;
    int iClipDiffEnd = 100 * 2.5;
    bool bValidByMate = false;
    int iClipPosByMate = -1;
    for(vector<St_SCRBound>::iterator itr = vLeftBoundSCR.begin(); itr != vLeftBoundSCR.end(); itr++)
    {
        //-->
        for(vector<St_BorderSCR>::iterator itrSCR = itr->vSCR.begin(); itrSCR != itr->vSCR.end(); itrSCR++)
        {
            if(itrSCR->bSecondMate)
            {
                if(abs(itrSCR->iInsertSize) < this->GetMinValidIS() ||
                   abs(itrSCR->iInsertSize) > this->GetMaxValidIS())
                    continue;

                if(itrSCR->iInsertSize > 0)
                {
                    if( //stDRGroup.iLeftBoundary - itrSCR->iMatPos > 0 &&
                        abs(stDRGroup.iLeftBoundary - itrSCR->iMatPos) < iMatDiffMax &&
                        itrSCR->iClipPos - (stDRGroup.iStart + stDRGroup.iStartAlignLen) > 0 &&
                        itrSCR->iClipPos - (stDRGroup.iStart + stDRGroup.iStartAlignLen) < iClipDiffStart)
                    {
                        if(iClipPosByMate < 0)
                            iClipPosByMate = itrSCR->iClipPos;
                        else if(iClipPosByMate < itrSCR->iClipPos) // need large one
                        {
                            iClipPosByMate = itrSCR->iClipPos;
                        }
                        bValidByMate = true;
                    }
                }
                else // 反向了
                {
                    if( abs(stDRGroup.iRightBoundary - itrSCR->iMatPos) < iMatDiffMax &&
                        itrSCR->iClipPos - (stDRGroup.iStart + stDRGroup.iStartAlignLen) > 0 &&
                        itrSCR->iClipPos - (stDRGroup.iStart + stDRGroup.iStartAlignLen) < iClipDiffStart)
                    {
                        if(iClipPosByMate < 0)
                            iClipPosByMate = itrSCR->iClipPos;
                        else if(iClipPosByMate < itrSCR->iClipPos) // need large one
                        {
                            iClipPosByMate = itrSCR->iClipPos;
                        }
                        bValidByMate = true;
                    }
                }
            }

            if(itrSCR->bFirstMate && itrSCR->iInsertSize < 0)  //这里好像缺少了 insert size大于 0的情况
            {
                if( //stDRGroup.iLeftBoundary - itrSCR->iMatPos > 0 &&
                    abs(stDRGroup.iLeftBoundary - itrSCR->iMatPos) < iMatDiffMax &&
                    itrSCR->iClipPos - (stDRGroup.iStart + stDRGroup.iStartAlignLen) > 0 &&
                    itrSCR->iClipPos - (stDRGroup.iStart + stDRGroup.iStartAlignLen) < iClipDiffStart)
                {
                    if(iClipPosByMate < 0)
                        iClipPosByMate = itrSCR->iClipPos;
                    else if(iClipPosByMate < itrSCR->iClipPos) // need large one
                    {
                        iClipPosByMate = itrSCR->iClipPos;
                    }
                    bValidByMate = true;
                }
            }
        }
        //<--

        if(itr->vSCR.size() > 1)
        {
            int iDiff = itr->iClipPos - stDRGroup.iStart;

            if(iDiff>0 && iDiff<iMaxDiff)
            {
                if(iValidClipPos < 0)
                    iValidClipPos = itr->iClipPos;
                else
                {
                    if(itr->iClipPos < iValidClipPos)
                        iValidClipPos = itr->iClipPos;
                }
            }
        }
    }
    if(iValidClipPos > 0)
    {
        stDRGroup.iStart = iValidClipPos;
        stDRGroup.bStartClipAdjusted = true;
        bLeftUpdated = true;
    }
    else
    {
        //logic 2: Check mate
        if(bValidByMate)
        {
            stDRGroup.iStart = iClipPosByMate;
            stDRGroup.bStartClipAdjusted = true;
            bLeftUpdated = true;
        }
    }

    //For right part
    iValidClipPos = -1;

    bValidByMate = false;
    iClipPosByMate = -1;
    for(vector<St_SCRBound>::iterator itr = vRightBoundSCR.begin(); itr != vRightBoundSCR.end(); itr++)
    {
        //-->
        for(vector<St_BorderSCR>::iterator itrSCR = itr->vSCR.begin(); itrSCR != itr->vSCR.end(); itrSCR++)
        {
            if(itrSCR->bFirstMate)
            {
                if(abs(itrSCR->iInsertSize) < this->GetMinValidIS() ||
                   abs(itrSCR->iInsertSize) > this->GetMaxValidIS())
                    continue;

                if(itrSCR->iInsertSize > 0)
                {
                    if( //itrSCR->iMatPos - stDRGroup.iRightBoundary > 0 &&
                        abs(itrSCR->iMatPos - stDRGroup.iRightBoundary) < iMatDiffMax &&
                        stDRGroup.iEnd - itrSCR->iClipPos > 0 &&
                        stDRGroup.iEnd - itrSCR->iClipPos < iClipDiffEnd)
                    {
                        if(iClipPosByMate < 0)
                            iClipPosByMate = itrSCR->iClipPos;
                        else if(iClipPosByMate > itrSCR->iClipPos) // need large one
                        {
                            iClipPosByMate = itrSCR->iClipPos;
                        }
                        bValidByMate = true;
                    }
                }
                else //反向的case
                {
                    if( abs(itrSCR->iMatPos - stDRGroup.iLeftBoundary) < iMatDiffMax &&
                        stDRGroup.iEnd - itrSCR->iClipPos > 0 &&
                        stDRGroup.iEnd - itrSCR->iClipPos < iClipDiffEnd)
                    {
                        if(iClipPosByMate < 0)
                            iClipPosByMate = itrSCR->iClipPos;
                        else if(iClipPosByMate > itrSCR->iClipPos) // need large one
                        {
                            iClipPosByMate = itrSCR->iClipPos;
                        }
                        bValidByMate = true;
                    }
                }
            }

            if(itrSCR->bSecondMate && itrSCR->iInsertSize < 0)
            {
                if( abs(itrSCR->iMatPos - stDRGroup.iRightBoundary) < iMatDiffMax &&
                    stDRGroup.iEnd - itrSCR->iClipPos > 0 &&
                    stDRGroup.iEnd - itrSCR->iClipPos < iClipDiffEnd)
                {
                    if(iClipPosByMate < 0)
                        iClipPosByMate = itrSCR->iClipPos;
                    else if(iClipPosByMate > itrSCR->iClipPos) // need large one
                    {
                        iClipPosByMate = itrSCR->iClipPos;
                    }
                    bValidByMate = true;
                }
            }
        }
        //<--

        if(itr->vSCR.size() > 1)
        {
            int iDiff = stDRGroup.iEnd - itr->iClipPos;

            if(iDiff>0 && iDiff<iMaxDiff)
            {
                if(iValidClipPos < 0)
                    iValidClipPos = itr->iClipPos;
                else
                {
                    if(itr->iClipPos > iValidClipPos)
                        iValidClipPos = itr->iClipPos;
                }
            }
        }
    }
    if(iValidClipPos > 0)
    {
        stDRGroup.iEnd = iValidClipPos;
        stDRGroup.bEndClipAdjusted = true;
        bRightUpdated = true;
    }
    else
    {
        //logic 2: Check mate
        if(bValidByMate)
        {
            stDRGroup.iEnd = iClipPosByMate;
            stDRGroup.bEndClipAdjusted = true;
            bRightUpdated = true;
        }
    }
}

void ClsSvDelDetect::FilterBadDR(vector<St_DiscordantReads>& vDiscdRreads)
{
    for(vector<St_DiscordantReads>::iterator itr = vDiscdRreads.end() -1; itr >= vDiscdRreads.begin(); itr--)
    {
        bool bGood = false;
        if(itr == vDiscdRreads.end() - 1) // this is the first one
        {
            if(itr - 2 >= vDiscdRreads.begin()) // the next and the next next is valid
            {
                if(abs(itr->iReadsPos - (itr-1)->iReadsPos) < m_fMeanInsertSize ||
                   abs(itr->iReadsPos - (itr-2)->iReadsPos) < m_fMeanInsertSize)
                {
                    if(abs(itr->iMatePos - (itr-1)->iMatePos) < m_fMeanInsertSize ||
                       abs(itr->iMatePos - (itr-2)->iMatePos) < m_fMeanInsertSize)
                    {
                        bGood = true;
                    }
                }
            }
            else if(itr - 1 >= vDiscdRreads.begin())
            {
                if(abs(itr->iReadsPos - (itr-1)->iReadsPos) < m_fMeanInsertSize &&
                   abs(itr->iMatePos - (itr-1)->iMatePos) < m_fMeanInsertSize)
                       bGood = true;
            }
            else
                bGood = true;

        }
        else if(itr == vDiscdRreads.begin()) // this is the last one
        {
            if(itr + 2 <= vDiscdRreads.end() -1) // the previous and the previous previous is valid
            {
                if(abs(itr->iReadsPos - (itr+1)->iReadsPos) < m_fMeanInsertSize ||
                   abs(itr->iReadsPos - (itr+2)->iReadsPos) < m_fMeanInsertSize )
                {
                    if(abs(itr->iMatePos - (itr+1)->iMatePos) < m_fMeanInsertSize ||
                       abs(itr->iMatePos - (itr+2)->iMatePos) < m_fMeanInsertSize)
                    {
                        bGood = true;
                    }
                }
            }
            else
            {
                bGood = true;
            }
        }
        else // this is the middle one (definitely have one before and one after)
        {
            if(abs(itr->iReadsPos - (itr-1)->iReadsPos) < m_fMeanInsertSize ||
               abs(itr->iReadsPos - (itr+1)->iReadsPos) < m_fMeanInsertSize )
            {
                if(abs(itr->iMatePos - (itr-1)->iMatePos) < m_fMeanInsertSize ||
                   abs(itr->iMatePos - (itr+1)->iMatePos) < m_fMeanInsertSize)
                {
                    bGood = true;
                }
            }
        }
        if(!bGood)
            vDiscdRreads.erase(itr);
    }
}

bool sort_scrClipPos_dup_func(St_BorderSCR stScr1, St_BorderSCR stScr2) // duplicate: same sort function already identified in clsparsebam.cpp
{
    if(stScr1.iClipPos < stScr2.iClipPos)
        return true;
    else
        return false;
}

bool sort_scrReadsMapPos_dup_func(St_BorderSCR stScr1, St_BorderSCR stScr2) // duplicate: same sort function already identified in clsparsebam.cpp
{
    if(stScr1.iReadsMapPos < stScr2.iReadsMapPos)
        return true;
    else
        return false;
}

void ClsSvDelDetect::ClusterSCReads(vector<St_SCRGroup>& vSCRGroup, vector<St_BorderSCR>& vBorderSCR)
{
    //Here --> we can add additional situation to make a better cluser --> Go!! *************************************
    vSCRGroup.clear();

    vector<St_BorderSCR> vNewSCR;
    UpgradeSCR(vBorderSCR, vNewSCR);

    //Use target clip reads for clustering complementary --> now we start to make clustering
    //Step 1: Sort first: from small to large based on clip position
    sort(vNewSCR.begin(), vNewSCR.end(), sort_scrReadsMapPos_dup_func);

    //step 2: Start grouping
    int iClipDiff = this->m_fMeanInsertSize + 3 * this->m_fStdDevInsertSize - 2 * 101; // the maximum length of
    int iInsetSizeDiff = 200; //this->m_fMeanInsertSize + 3 * this->m_fStdDevInsertSize; //3 * this->m_fStdDevInsertSize;

    ///1: we first group the boundary reads --> Go
    vector<St_GroupBound> vGroupBound;
    St_GroupBound stGroupBound;

    for(vector<St_BorderSCR>::iterator itr = vNewSCR.begin(); itr < vNewSCR.end(); itr++)
    {
        stGroupBound.Clear();
        stGroupBound.vSCR.push_back(*itr);
        //En_ClipPart enClipPart = itr->GetClipPart();

        for(vector<St_BorderSCR>::iterator subItr = itr+1; subItr < vNewSCR.end(); subItr++)
        {
            if(abs(itr->iClipPos - subItr->iClipPos) <= iClipDiff &&
               abs(abs(itr->iInsertSize) - abs(subItr->iInsertSize)) < iInsetSizeDiff)
                    //&& subItr->GetClipPart() == enClipPart
            {
                stGroupBound.vSCR.push_back(*subItr);                
                itr = subItr;
            }
            else
            {                
                break;
            }
        }
        vGroupBound.push_back(stGroupBound);
    }


    //Filter isolate clip
    for(vector<St_GroupBound>::iterator itr = vGroupBound.end() - 1; itr >= vGroupBound.begin(); itr--)
    {
        if(itr->vSCR.size() < 4)
            vGroupBound.erase(itr);
    }
    (*m_pOfs) << "After Filter: " << to_string(vGroupBound.size()) << endl << endl;

    GetSCRGroup(vGroupBound, vSCRGroup);       
}

void ClsSvDelDetect::UpgradeSCR(vector<St_BorderSCR>& vBorderSCR, vector<St_BorderSCR>& vNewSCR)
{
    //Only keep 4 types of softclip and set the type value for each specific one
    vNewSCR.clear();
    int iMat1End = 0;
    int iMat2End = 0;
    int iMat1Start = 0;
    int iMat2Start = 0;

    int iDrippedOutofDRGroup = 0;
    for(vector<St_BorderSCR>::iterator itr = vBorderSCR.begin(); itr != vBorderSCR.end(); itr++)
    {
        if(itr->bDrippedIntoDRGroup) // We do not consider the soft clip reads which dropped into the range of discorcdant reads group
            continue;
        else
            iDrippedOutofDRGroup++;        

        //Type 1
        //Type 3
        if(itr->bFirstMate)
        {
            //Type 1
            if(itr->GetClipPart() == cpRight &&
               itr->bMateMapped &&
               itr->iReadsMapPos < itr->iClipPos && //  be cut off here
               abs(itr->iInsertSize) > m_fMeanInsertSize &&
               abs(itr->iInsertSize) < m_fMeanInsertSize + 3 * m_fStdDevInsertSize) //small deletion
            {
                itr->enMapType = mtMat1End;
                vNewSCR.push_back(*itr);
                iMat1End++;
            }

            //Type 2
            if(itr->GetClipPart() == cpLeft &&
               itr->bMateMapped &&
               itr->iReadsMapPos == itr->iClipPos &&
               abs(itr->iInsertSize) > m_fMeanInsertSize - 3*m_fStdDevInsertSize &&
               abs(itr->iInsertSize) < m_fMeanInsertSize + 3*m_fStdDevInsertSize)
            {
                itr->enMapType = mtMat1Start;
                vNewSCR.push_back(*itr);
                iMat1Start++;
            }
        }
        //Type 2
        //Type 4
        else if(itr->bSecondMate)
        {
            //Type 2
            if(itr->GetClipPart() == cpLeft &&
               itr->bMateMapped &&
               itr->iReadsMapPos == itr->iClipPos &&
               abs(itr->iInsertSize) > m_fMeanInsertSize &&
               abs(itr->iInsertSize) < m_fMeanInsertSize + 3*m_fStdDevInsertSize) //small deletion
            {
                itr->enMapType = mtMat2Start;
                vNewSCR.push_back(*itr);
                iMat2Start++;
            }

            //Type 4
            if(itr->GetClipPart() == cpRight &&
               itr->bMateMapped &&
               itr->iReadsMapPos < itr->iClipPos &&
               abs(itr->iInsertSize) > m_fMeanInsertSize - 3*m_fStdDevInsertSize &&
               abs(itr->iInsertSize) < m_fMeanInsertSize + 3*m_fStdDevInsertSize)
            {
                itr->enMapType = mtMat2End;
                vNewSCR.push_back(*itr);
                iMat2End++;
            }
        }
        else{}
    }
    (*m_pOfs) << "vBorderSCR Size: " << to_string(vBorderSCR.size()) << endl;
    (*m_pOfs) << "vNewSCR Size   : " << to_string(vNewSCR.size()) << endl;
    (*m_pOfs) << "\t iMate1End: " << to_string(iMat1End) << endl;
    (*m_pOfs) << "\t iMate1Start: " << to_string(iMat1Start) << endl;
    (*m_pOfs) << "\t iMate2End: " << to_string(iMat2End) << endl;
    (*m_pOfs) << "\t iMate2Start: " << to_string(iMat2Start) << endl;
    (*m_pOfs) << "iDrippedOutofDRGroup: " << to_string(iDrippedOutofDRGroup) << endl;
}

void ClsSvDelDetect::GetSCRGroup(vector<St_GroupBound>& vGroupBound, vector<St_SCRGroup>& vSCRGroup)
{
    vSCRGroup.clear();

    //--> Go!!
    //1: get all of left side reads and right side reads from the boundary data set
    St_SCRGroup stSCRGroup;
    for(vector<St_GroupBound>::iterator itr = vGroupBound.begin(); itr != vGroupBound.end(); itr++)
    {
        stSCRGroup.Clear();
        if(GetGroupFromBoundary(itr->vSCR, stSCRGroup))
        {
            vSCRGroup.push_back(stSCRGroup);
        }
    }
    (*m_pOfs) << "vSCRGroup: " << IntToStr(vSCRGroup.size()) << endl;
}

bool ClsSvDelDetect::GetGroupFromBoundary(vector<St_BorderSCR>& vSCR, St_SCRGroup& stSCRGroup) // record all of reads which made contribution -->
{
    if(vSCR.size() <= 2)
        return false;

    //Get real left boundary
    //First separate all the reads into 2 groups
    vector<St_BorderSCR> vStartBound;
    vector<St_BorderSCR> vEndBound;
    for(vector<St_BorderSCR>::iterator itr = vSCR.begin(); itr != vSCR.end(); itr++)
    {
        switch(itr->enMapType)
        {
            case mtMat1End:
                vStartBound.push_back(*itr);
                break;
            case mtMat1Start:
                vEndBound.push_back(*itr);
                break;
            case mtMat2End:
                vStartBound.push_back(*itr);
                break;
            case mtMat2Start:
                vEndBound.push_back(*itr);
                break;
            default:
                break;
        }
    }

    (*m_pOfs) << " --> Start G From Boundary" << endl;
    (*m_pOfs) << "vStartBound_1: " << IntToStr(vStartBound.size()) << endl;
    (*m_pOfs) << "vEndBound_1  : " << IntToStr(vEndBound.size()) << endl;

    if(vEndBound.empty() || vStartBound.empty())
    {
        (*m_pOfs)  << "Direct BAD!" << endl;
        (*m_pOfs)  << " <---- End" << endl << endl;
        return false;
    }

    //Filter the abnormal clip reads
//    sort(vStartBound.begin(), vStartBound.end(), sort_scrClipPos_func);
//    sort(vEndBound.begin(), vEndBound.end(), sort_scrClipPos_func);
    (*m_pOfs) << "sort both StartBound and EndBound Finished!" << endl;
    int iIndex = vStartBound.size()/2;
    int iMidClipPos = vStartBound[iIndex].iClipPos;
    int iOffSet = 100;
    for(vector<St_BorderSCR>::iterator itr = vStartBound.end() - 1; itr >=vStartBound.begin(); itr--)
    {
        if(itr->iClipPos >= iMidClipPos + iOffSet ||
           itr->iClipPos <= iMidClipPos - iOffSet)
            vStartBound.erase(itr);
    }


    iIndex = vEndBound.size()/2;
    if(vEndBound.size() % 2 == 0)
        iIndex--;
    iMidClipPos = vEndBound[iIndex].iClipPos;
    for(vector<St_BorderSCR>::iterator itr = vEndBound.end() - 1; itr >=vEndBound.begin(); itr--)
    {
        if(itr->iClipPos >= iMidClipPos + iOffSet ||
           itr->iClipPos <= iMidClipPos - iOffSet)
            vEndBound.erase(itr);
    }

    if(vStartBound.size() < 2 || vEndBound.size() < 2)//if(vStartBound.size() + vEndBound.size() <= 2)
    {
        (*m_pOfs)  << "After Filter BAD!" << endl;
        (*m_pOfs)  << " <---- End" << endl << endl;
        return false;
    }

    //Get Clip Pos
    //1: Start
    iIndex = vStartBound.size()/2;
    stSCRGroup.iStart = vStartBound[iIndex].iClipPos;//(vStartBound.end() - 1)->iClipPos;
    stSCRGroup.iLeftBoundary = vStartBound.begin()->iClipPos;

    //2: End
    iIndex = vEndBound.size()/2;
    //if(vEndBound.size() % 2 == 0)
    //    iIndex--;
    stSCRGroup.iEnd = vEndBound[iIndex].iClipPos;//vEndBound.begin()->iClipPos;
    stSCRGroup.iRightBoundary = (vEndBound.end() - 1)->iClipPos;

    (*m_pOfs) << "stSCRGroup.iStart: " << IntToStr(stSCRGroup.iStart) << endl;
    (*m_pOfs) << "stSCRGroup.iEnd  : " << IntToStr(stSCRGroup.iEnd) << endl;

    //Only keep size of small deletion
    if(stSCRGroup.iEnd > stSCRGroup.iStart &&
       abs(stSCRGroup.iEnd - stSCRGroup.iStart) >= 40 &&
       abs(stSCRGroup.iEnd - stSCRGroup.iStart) <= m_fMeanInsertSize + 3 * m_fStdDevInsertSize)
    {
        (*m_pOfs) << "vStartBound_2: " << IntToStr(vStartBound.size()) <<  " >>> ";
        for(vector<St_BorderSCR>::iterator itr = vStartBound.begin(); itr != vStartBound.end(); itr++)
        {
            (*m_pOfs) << "<" << to_string(itr->iReadsMapPos) << ", " << to_string(itr->iClipPos) << "> ";
        }
        (*m_pOfs) << endl;

        (*m_pOfs) << "vEndBound_2  : " << IntToStr(vEndBound.size()) << " >>> ";
        for(vector<St_BorderSCR>::iterator itr = vEndBound.begin(); itr != vEndBound.end(); itr++)
        {
            (*m_pOfs) << "<" << to_string(itr->iReadsMapPos) << ", " << to_string(itr->iClipPos) << "> ";
        }
        (*m_pOfs) << endl;

        (*m_pOfs) << "Yes, I am good" << endl;
        //Save boundary info into Group Info -->
        stSCRGroup.vStartSCR.insert(stSCRGroup.vStartSCR.end(), vStartBound.begin(), vStartBound.end());
        stSCRGroup.vEndSCR.insert(stSCRGroup.vEndSCR.end(), vEndBound.begin(), vEndBound.end());
        //<--
        (*m_pOfs)  << " <---- End" << endl << endl;
        return true;
    }
    else
    {
        (*m_pOfs) << "Group Range BAD" << endl;
        (*m_pOfs)  << " <---- End" << endl << endl;
        return false;
    }
}

void ClsSvDelDetect::DepthFilter(vector<St_DRGroup>& vDRGroup, vector<St_SCRGroup>& vSCRGroup,
                                 string& strBamFile, string strChrom)
{
    //for soft clip reads first
    (*m_pOfs) << "vSCRGroup Size: " << to_string(vSCRGroup.size()) << endl;
    if(!vSCRGroup.empty())
    {
        for(vector<St_SCRGroup>::iterator itr = vSCRGroup.end() - 1; itr >= vSCRGroup.begin(); itr--)
        {
            if(DepthFilterGroup(*itr, strBamFile, m_fAvgDepth, strChrom))
                vSCRGroup.erase(itr);
        }
        (*m_pOfs) << "\t " << to_string(vSCRGroup.size()) << endl;
    }

    //for discordant reads second
    (*m_pOfs) << "vDRGroup Size: " << to_string(vDRGroup.size()) << endl;
    if(!vDRGroup.empty())
    {
        for(vector<St_DRGroup>::iterator itr = vDRGroup.end() - 1; itr >= vDRGroup.begin(); itr--)
        {
            if(DepthFilterGroup(*itr, strBamFile, m_fAvgDepth, strChrom))
                vDRGroup.erase(itr);
        }
        (*m_pOfs) << "\t" << to_string(vDRGroup.size()) << endl;
    }
}

bool ClsSvDelDetect::DepthFilterGroup(St_BaseGroup& stGroup, string& strBamFile,
                                      float fThreshold, string strChrom)
{
    string strCmd = (string)"samtools depth -a -r " +
                    strChrom + ":" + to_string(stGroup.iStart) + "-" + to_string(stGroup.iEnd) + " " + strBamFile +
                    " | awk '{sum+=$3} END {if(sum==\"\") {print 0} else {print sum/NR}}'";

    string strResult = exec(strCmd.c_str());
    float fCoverage = stof(strResult.c_str());

    //if(fCoverage < fThreshold || fCoverage > 2 * fThreshold)
    if(fCoverage < fThreshold) //0.9 -> 0.5110 <> 0.8 -- 0.5217
        return false;
    else
        return true;
}

void ClsSvDelDetect::GetSpeciReads(vector<St_DRGroup>& vDRGroup, string strBamFile)
{
    //Collect the speci pair end reads: one is mapped and one is unmapped
    BamReader* pBamReader = new BamReader();

    pBamReader->Open(strBamFile);
    pBamReader->OpenIndex(strBamFile + ".bai");
    BamAlignment al;
    St_RegReads stRegReads;

    while(pBamReader->GetNextAlignment(al)) //Get each alignment result
    {
        if(!al.IsMapped() || al.QueryBases == "") //　we need both map
            continue;

        if(al.IsMapped() && !al.IsMateMapped())
        {}
        else
            continue;

        if(al.QueryBases.find('N') != string::npos ||
           al.QueryBases.find('n') != string::npos)
            continue;

        for(vector<St_DRGroup>::iterator itr = vDRGroup.begin(); itr != vDRGroup.end(); itr++)
        {
            if(al.Position >= itr->iLeftBoundary && al.Position <= itr->iRightBoundary) // Dropped into the current group
            {
                stRegReads.Clear();
                stRegReads.strName = al.Name;
                stRegReads.iReadsMapPos = al.Position;
                itr->vRegReads.push_back(stRegReads);
            }
        }
    }

    delete pBamReader;
    pBamReader = NULL;

    //output the the speci reads dropped in each range -->
    (*m_pOfs) << "********* The speci reads for each group *************" << endl;
    for(vector<St_DRGroup>::iterator itr = vDRGroup.begin(); itr != vDRGroup.end(); itr++)
    {
        (*m_pOfs) << "(" << to_string(itr->iStart) << ", " << to_string(itr->iEnd) << ")"
             << " --- " << to_string(itr->vRegReads.size()) << endl;
    }
    (*m_pOfs) << "*******************" << endl;
    //<--
}

void ClsSvDelDetect::CollectFeatures(vector<St_DRGroup>& vDRGroup, string strBamFile, string strChrom)
{
    //Now we try to collect features for each group -->
    (*m_pOfs) << "vDRGroup Sum: " << to_string(vDRGroup.size()) << endl;

    for(vector<St_DRGroup>::iterator itr = vDRGroup.begin(); itr != vDRGroup.end(); itr++)
    {
        (*m_pOfs) << endl << "---" << "(" << to_string(itr->iStart) << ", " << to_string(itr->iEnd) << ")" << endl;

        itr->stGF.iStart = itr->iStart;
        itr->stGF.iEnd = itr->iEnd;
        itr->stGF.iSpeciReadsNum = itr->vRegReads.size();
        itr->stGF.iDRReadsNum = itr->vDR.size();
        //Try to get abnormal rang of depth --> Let's do it tomorrow        
        GetDepthFeatures(*itr, strBamFile, strChrom);
        (*m_pOfs) << "\t GetDepthFeatures" << endl;

        //Get Clip features -->
        GetClipFeatures(*itr);
        (*m_pOfs) << "\t GetClipFeatures" << endl;
    }
    //<--

    //Next Step --> make some additional filtering --> Go!!
    (*m_pOfs) << endl << "***************** Temp Filter by Features *******************" << endl;
    (*m_pOfs) << "Number of Groups: " << to_string(vDRGroup.size()) << endl;
    for(vector<St_DRGroup>::iterator itr = vDRGroup.end() - 1; itr >= vDRGroup.begin(); itr--)
    {
//        float fMaxDepth = itr->stGF.GetMaxTargetDepth();
//        float fMinDepth = itr->stGF.GetMinTargetDepth();

//        if(fMinDepth < 0 || fMaxDepth < 0)
//        {
//            vDRGroup.erase(itr);
//            continue;
//        }

//        if(fMinDepth <= m_fAvgDepth * .3 ||
//           fMaxDepth >= m_fAvgDepth * 1.5)
//        {}
//        else
//            vDRGroup.erase(itr);
    }
    (*m_pOfs) << "\t After feature Filter: "  << to_string(vDRGroup.size()) << endl;
}

//Go!!  --> finish  those depths when weak up !!!
void ClsSvDelDetect::GetDepthFeatures(St_DRGroup& stDRGroup, string strBamFile, string strChrom)
{
    //Step 1: Get the avarage depth in this range --> Go First
    //  (1)--> Left_Reg_Depth
//    int iOffSet = 50;
//    int iL = stDRGroup.iLeftBoundary - iOffSet;
//    int iR = stDRGroup.iStart + stDRGroup.iStartAlignLen;
//    if(stDRGroup.bStartClipAdjusted)
//        iR = stDRGroup.iStart;

//    string strCmd = (string)"samtools depth -a -r " +
//                    "11:" + to_string(iL) + "-" + to_string(iR) +
//                    " " + strBamFile + " | awk '{sum+=$3} END {if(sum==\"\") {print 0} else {print sum/NR}}'";
//    string strResult = exec(strCmd.c_str());
//    float fCoverage = stof(strResult.c_str());
//    stDRGroup.stGF.fAvgRegDepthLeft = fCoverage;

    int iOffSet = 50;
    int iL = stDRGroup.iLeftBoundary - iOffSet;
    int iR = stDRGroup.iStart + stDRGroup.iStartAlignLen;
    if(stDRGroup.bStartClipAdjusted)
        iR = stDRGroup.iStart;
    int iLen = iR - iL + 1;

    (*m_pOfs) << "m_fAvgDepth is: " << FloatToStr(m_fAvgDepth) << endl;
    //This is for the left regular discordant range        
    GetDepthOfRangeEX(strBamFile,
                    iL, iR, iLen,
                    stDRGroup.stGF.stDCRangeLeft, strChrom);


    //  (2) --> Right_Reg_Depth
//    iL = stDRGroup.iEnd;
//    iR = stDRGroup.iRightBoundary + iOffSet;
//    strCmd = (string)"samtools depth -a -r " +
//             "11:" + to_string(iL) + "-" + to_string(iR) +
//             " " + strBamFile + " | awk '{sum+=$3} END {if(sum==\"\") {print 0} else {print sum/NR}}'";

//    strResult = exec(strCmd.c_str());
//    fCoverage = stof(strResult.c_str());
//    stDRGroup.stGF.fAvgRegDepthRight = fCoverage;

    iOffSet = 50;
    iL = stDRGroup.iEnd;
    iR = stDRGroup.iRightBoundary + iOffSet;
    iLen = iR - iL + 1;

    //This is for the right regular discordant range
    GetDepthOfRangeEX(strBamFile,
                    iL, iR, iLen,
                    stDRGroup.stGF.stDCRangeRight, strChrom);


    //Step 2: Get the average depth for the real normal depth range
    //(1) Left No SV range        
    iOffSet = 50;
    iL = stDRGroup.iLeftBoundary - iOffSet*3;
    iR = stDRGroup.iLeftBoundary - iOffSet;
    iLen = iR - iL + 1;
    //--> 1: Get the same features as other DC ranges -->
    GetDepthOfRangeEX(strBamFile,
                    iL, iR, iLen,
                    stDRGroup.stGF.stNoSvRangeLeft, strChrom);
    //--> 2: Calc average Depth
    string strCmd = (string)"samtools depth -a -r " +
                    strChrom + ":" + to_string(iL) + "-" + to_string(iR) +
                    " " + strBamFile + " | awk '{sum+=$3} END {if(sum==\"\") {print 0} else {print sum/NR}}'";
    string strResult = exec(strCmd.c_str());
    float fCoverage = stof(strResult.c_str());
    stDRGroup.stGF.stNoSvRangeLeft.fAvgDepth = fCoverage;

    //(2) Right No SV range
    iOffSet = 50;
    iL = stDRGroup.iRightBoundary + iOffSet;
    iR = stDRGroup.iRightBoundary + iOffSet*3;
    iLen = iR - iL + 1;
    //--> 1: Get the same features as other DC ranges -->
    GetDepthOfRangeEX(strBamFile,
                    iL, iR, iLen,
                    stDRGroup.stGF.stNoSvRangeRight, strChrom);
    //--> 2: Calc average Depth
    strCmd = (string)"samtools depth -a -r " +
             strChrom + ":" + to_string(iL) + "-" + to_string(iR) +
             " " + strBamFile + " | awk '{sum+=$3} END {if(sum==\"\") {print 0} else {print sum/NR}}'";
    strResult = exec(strCmd.c_str());
    fCoverage = stof(strResult.c_str());
    stDRGroup.stGF.stNoSvRangeRight.fAvgDepth = fCoverage;


    //Step 3: Find the partial depth range --> Go!! -->
//    iOffSet = 50;
//    iL = stDRGroup.iLeftBoundary - iOffSet;
//    iR = stDRGroup.iRightBoundary + iOffSet;
//    GetDepthOfRange(strBamFile,
//                    iL, iR, stDRGroup.GetAccuLen(),
//                    stDRGroup.stGF.fPartialDepth, stDRGroup.stGF.iDepthRangLen);
    iOffSet = 0;
    iL = stDRGroup.iStart + stDRGroup.iStartAlignLen - iOffSet;
    if(stDRGroup.bStartClipAdjusted)
        iL = stDRGroup.iStart - iOffSet;

    iR = stDRGroup.iEnd + iOffSet;
    iLen = iR - iL + 1;
    GetDepthOfRangeEX(strBamFile,
                    iL, iR, iLen,
                    stDRGroup.stGF.stDCRangeCore, strChrom);
    //<--
}

struct St_DepthRange
{
    int iDepth;
    int iLen;

    St_DepthRange()
    {
        Clear();
    }

    void Clear()
    {
        iDepth = -1;
        iLen = 0;
    }
};

bool sort_DepthSmallToLarge_func(St_DepthRange st1, St_DepthRange st2)
{
    if(st1.iDepth < st2.iDepth)
        return true;
    else
        return false;
}

/*其实在这里deletion存在两种情况:
//(1) 全部被delete 掉了  -->这个时候就应该都是0
//(2) 只有一个copy 被delete 掉了 --> 这个时候应该是 average depth的一半
//所以我们应该分情况去做 的去做这个事儿
//同时如果两边比较messy
  (1)存在 repeats --> 这样的话  错误的跑到了deletion的里面
     找高coverage 的区域，coverage的书程度应该是 standard的 1.5 倍

我们期待的不同 --> drop到 discordant中的和normal range的区别
(1) 突然很高的Coverage
(2) 突然很低的coverage

正常部分我就期待的应该是很正常才对

--> 这个函数我们应该重新撰写一下
*/
void ClsSvDelDetect::GetDepthOfRangeEX(string strBamFile,
                                       int iStartPos, int iEndPos, int iRangeLen,
                                       St_DepthSeg& stDepthSeg,
                                       string strChrom)
{
    if(strChrom == "")
        return;

    /*
    //Step 1: Get The coverage in every point for current range
    //Step 2: 将他们进行分组 --> 合并 merge the depth with same value
    //Step 3: 去掉一些孤岛: Filter the isolate depth value
    //Step 4: 将他们进行排序
    //Step 5: 我们得到三个相应的depth的值(lowest, highest, average)
            　目前为止，我们需要在长度上做一些限制，使得能够得到更准确的值
    */
    stDepthSeg.iSegLen = iRangeLen;

    //Step 1:
    string strCmd = (string)"samtools depth -a -r " +
             strChrom + ":" + to_string(iStartPos) + "-" + to_string(iEndPos) +
             " " + strBamFile + " | awk '{print $3}'";

    string strResult = exec(strCmd.c_str());
    vector<int> vCoverage;
    int iStart = 0;
    int iLen = 0;
    for(int i=0; i<strResult.length(); i++)
    {
        if(strResult.substr(i, 1) == "\n")
        {
            iLen = i - iStart;
            vCoverage.push_back(atoi(strResult.substr(iStart, iLen).c_str()));
            iStart = i + 1;
        }
    }

    //Step 2
    vector<St_DepthRange> vDepthRange;
    St_DepthRange stDepthRange;
    for(vector<int>::iterator itr = vCoverage.begin(); itr != vCoverage.end(); itr++)
    {
        stDepthRange.Clear();
        stDepthRange.iDepth = *itr;
        for(vector<int>::iterator subItr = itr; subItr != vCoverage.end(); subItr++)
        {
            if(stDepthRange.iDepth == *subItr)
            {
                stDepthRange.iLen++;
            }
            else
            {
                vDepthRange.push_back(stDepthRange);
                itr = subItr - 1;
                break;
            }

            if(subItr + 1 == vCoverage.end()) // the last one
            {
                vDepthRange.push_back(stDepthRange);
                itr = subItr;
                break;
            }
        }
    }

    //Step 3: Filter the isolate depth value --> 这一步我觉得不需要，因为这些并不影响大局

    //Step 4: Sort those depth
    sort(vDepthRange.begin(), vDepthRange.end(), sort_DepthSmallToLarge_func);

    for(vector<St_DepthRange>::iterator itr = vDepthRange.begin(); itr != vDepthRange.end(); itr++)
    {
        (*m_pOfs) << to_string(itr->iDepth) << endl;
        (*m_pOfs) << "\t\t" << to_string(itr->iLen) << endl;
    }
    (*m_pOfs) << to_string(vDepthRange.size()) << endl;

    //Step 5:Get the depth value: Max, Average, Min
    //1: For deletion 0/0
    int iZeroPureDepthLen = 0;    
    int iZeroWeekDepthLen = 0;
    int iZeroMaxDiff = floor(.25 * m_fAvgDepth); //Allow 1 to 1/4 depth (25% average  depth)
//    int iZeroMaxDiff = ceil(.25 * m_fAvgDepth); //Allow 1 to 1/4 depth (25% average depth)
    int iZeroWeekDepthSum = 0;
    bool bZeroFind = false;

    //2: For deletion 0/1 or 1/0
    int iZeroMixDepthMin = floor(.25 * m_fAvgDepth);
    int iZeroMixDepthMax = ceil(.5 * m_fAvgDepth);
//    int iZeroMixDepthMin = .25 * m_fAvgDepth;
//    int iZeroMixDepthMax = .5 * m_fAvgDepth + 1;
    int iZeroMixDepthLen = 0; // from depth
    int iZeroMixDepthSum = 0;
    bool bZeroMixFind = false;

    //3: For Repearts area
    //   *For repeats area --> which mean the repeats drops into the deletion range and mis mapped by mapping tool
    int iRepeatDepth = m_fAvgDepth * 1.5; // --> this may one copy (from only one strand)
    int iStrongRepeatDepth = m_fAvgDepth * 2;
    int iRepeatDepthSum = 0;
    int iStrongRepeatDepthSum = 0;
    int iRepeatDepthLen = 0;
    int iStrongRepeatDepthLen = 0;

    //Now already from small depth to large -->
    for(vector<St_DepthRange>::iterator itr = vDepthRange.begin(); itr != vDepthRange.end(); itr++)
    {
        //1: Zero / Zero Week
        if(itr->iDepth == 0)
        {
            iZeroPureDepthLen += itr->iLen;
            iZeroWeekDepthLen += itr->iLen;
        }
        else if(itr->iDepth <= iZeroMaxDiff)
        //else if(itr->iDepth < iZeroMaxDiff)
        {
            iZeroWeekDepthLen += itr->iLen;
            iZeroWeekDepthSum += itr->iDepth * itr->iLen;
        }

        //2: ZeroMix
        if(itr->iDepth > iZeroMixDepthMin && itr->iDepth <= iZeroMixDepthMax)
        //if(itr->iDepth >= iZeroMixDepthMin && itr->iDepth <= iZeroMixDepthMax)
        {
            iZeroMixDepthLen += itr->iLen;
            iZeroMixDepthSum += itr->iDepth * itr->iLen;
        }

        //3: Repeats
        if(itr->iDepth >= iRepeatDepth)
        {
            iRepeatDepthLen += itr->iLen;
            iRepeatDepthSum += itr->iDepth * itr->iLen;
        }

        if(itr->iDepth >= iStrongRepeatDepth)
        {
            iStrongRepeatDepthLen += itr->iLen;
            iStrongRepeatDepthSum += itr->iDepth * itr->iLen;
        }
    }

    //Set value to the structure -> Go!!
    stDepthSeg.iZeroPureDepthLen = iZeroPureDepthLen;

    stDepthSeg.fZeroWeekDepth = (iZeroWeekDepthLen != 0 ? (float)iZeroWeekDepthSum / iZeroWeekDepthLen : -1.);
    stDepthSeg.iZeroWeekDepthLen = iZeroWeekDepthLen;

    stDepthSeg.fZeroMixDepth = (iZeroMixDepthLen != 0 ? (float)iZeroMixDepthSum / iZeroMixDepthLen : -1.);
    stDepthSeg.iZeroMixDepthLen = iZeroMixDepthLen;

    stDepthSeg.fRepWeekDepth = (iRepeatDepthLen != 0 ? (float)iRepeatDepthSum / iRepeatDepthLen : -1.);
    stDepthSeg.iRepWeekDepthLen = iRepeatDepthLen;

    stDepthSeg.fRepStrongDepth = (iStrongRepeatDepthLen != 0 ? (float)iStrongRepeatDepthSum / iStrongRepeatDepthLen : -1.);
    stDepthSeg.iRepStrongDepthLen = iStrongRepeatDepthLen;
    //<--

    //下一步
    //For deletion 0/0
    if(iZeroPureDepthLen > iRangeLen * .15)
    {
        stDepthSeg.fStatDepth = 0;
        stDepthSeg.iStatDepthLen = iZeroPureDepthLen;
        stDepthSeg.fStatLenRatio = (float)stDepthSeg.iStatDepthLen / stDepthSeg.iSegLen;
        stDepthSeg.enDepthType = dtZeroPure;
        return;
    }

    //For repeats
    if(iStrongRepeatDepthLen > iRangeLen * .15)
    {
        stDepthSeg.fStatDepth = (float)iStrongRepeatDepthSum / iStrongRepeatDepthLen;
        stDepthSeg.iStatDepthLen = iStrongRepeatDepthLen;
        stDepthSeg.fStatLenRatio = (float)stDepthSeg.iStatDepthLen / stDepthSeg.iSegLen;
        stDepthSeg.enDepthType = dtRepStrong;
        return;
    }

    if(iZeroWeekDepthLen > iRangeLen * .30)
    {
        stDepthSeg.fStatDepth = (float)iZeroWeekDepthSum / iZeroWeekDepthLen;
        stDepthSeg.iStatDepthLen = iZeroWeekDepthLen;
        stDepthSeg.fStatLenRatio = (float)stDepthSeg.iStatDepthLen / stDepthSeg.iSegLen;
        stDepthSeg.enDepthType = dtZeroWeek;
        return;
    }

    //For deletion 0/1
    if(iZeroMixDepthLen > iRangeLen * .30)
    {
        stDepthSeg.fStatDepth = (float)iZeroMixDepthSum / iZeroMixDepthLen;
        stDepthSeg.iStatDepthLen = iZeroMixDepthLen;
        stDepthSeg.fStatLenRatio = (float)stDepthSeg.iStatDepthLen / stDepthSeg.iSegLen;
        stDepthSeg.enDepthType = dtZeroMix;
        return;
    }

    if(iRepeatDepthLen > iRangeLen * .30)
    {
        stDepthSeg.fStatDepth = (float)iRepeatDepthSum / iRepeatDepthLen;
        stDepthSeg.iStatDepthLen = iRepeatDepthLen;
        stDepthSeg.fStatLenRatio = (float)stDepthSeg.iStatDepthLen / stDepthSeg.iSegLen;
        stDepthSeg.enDepthType = dtRepWeek;
        return;
    }
    //<--

//    //Try to find the largest and the longest one ->
//    int iIndex = -1; // the index of the largest depth
//    int iDepth = -1; // largest depth
//    iLen = -1; //the length of the largest depth
//    //<--

//    int iZeroDepthLen = 0;

//    //--> we also need to find the longest 0 range -->  Go!!!
//    int iSmallestLongestIndex = -1;
//    int iSmallestDepth = 10000;
//    int iSmallestLongestLen = -1;
//    vector<int> vZeroIndex;
//    //<--
//    for(unsigned int i=0; i<vDepthRange.size(); i++)
//    {
//        if(vDepthRange[i].iDepth == 0)
//        {
//            iZeroDepthLen += vDepthRange[i].iLen;
//            vZeroIndex.push_back(i);
//        }

//        //Get the smallest and longest index -->
//        if(vDepthRange[i].iDepth < iSmallestDepth)
//        {
//            iSmallestDepth = vDepthRange[i].iDepth;
//            iSmallestLongestIndex = i;
//            iSmallestLongestLen = vDepthRange[i].iLen;
//        }
//        else if(vDepthRange[i].iDepth == iSmallestDepth)
//        {
//            if(vDepthRange[i].iLen > iSmallestLongestLen)
//            {
//                iSmallestLongestIndex = i;
//                iSmallestLongestLen = vDepthRange[i].iLen;
//            }
//        }
//        //<--

//        if(iDepth <= vDepthRange[i].iDepth) // update need to find the lagest one !! so change < to <=
//        {
//            bool bUpdate = true;
//            if(iDepth == vDepthRange[i].iDepth &&
//               vDepthRange[i].iLen < iLen)
//                bUpdate = false;

//            if(bUpdate)
//            {
//                iDepth = vDepthRange[i].iDepth;
//                iLen = vDepthRange[i].iLen;
//                iIndex = i;
//            }
//        }
//    }
//    if(iIndex < 0)
//    {
//        cout << "Bad Depth Calculation" << endl;
//        return;
//    }

//    cout << "\t\t Get target Depth" << endl;
//    cout << "iZeroDepthLen: " << to_string(iZeroDepthLen) << endl;
//    cout << "iSmallestDepth: " << to_string(iSmallestDepth) << endl;
//    cout << "iSmallestLongestIndex: " << to_string(iSmallestLongestIndex) << endl;
//    cout << "iSmallestLongestLen: " << to_string(iSmallestLongestLen) << endl;

//    float fRatio = 0.25;
//    //int iRangeLen = stDRGroup.iEnd - stDRGroup.iStart + 1;
//    //Check the range of 0 coverage: special treatment shuld be applied (may be discontinously) -->
//    bool bLowDepthExtract = false;
//    float fAvgLowDepth = -1;
//    int iLowDepthLen = 0;
//    if(iZeroDepthLen >= iRangeLen * fRatio)
//    {
//        fTargetRangDepth = 0;
//        iTargetRangLen = iZeroDepthLen;

//        cout << "Set Zero Depth Length!" << endl;
//        bLowDepthExtract = true;
//    }
//    else if(iSmallestDepth == 0 && iSmallestLongestIndex >= 0) //Get
//    {
//        int iPreIndex = 0;
//        int iTotalDepth = 0;
//        vector<int> vRangedIndex;
//        for(vector<int>::iterator itr = vZeroIndex.begin(); itr != vZeroIndex.end(); itr++)
//        {
//            int iCurIndex = *itr;

//            if(iCurIndex < iPreIndex)
//                continue;

//            vRangedIndex.clear();
//            //We allow 1 diff
//            int iMaxAllowed = 2; //change the coefficient from 1 to 2
//            int iJumpNode = 0;
//            int iMaxJumNode = 1;
//            int iAccidentLen = 2;
//            for(int i = iCurIndex - 1; i >= iPreIndex; i--)
//            {
//                if( abs(vDepthRange[i].iDepth - iSmallestDepth) <= iMaxAllowed)
//                {
//                    vRangedIndex.push_back(i);
//                    iJumpNode = 0;
//                }
//                else
//                {
//                    if(vDepthRange[i].iLen > iAccidentLen) // it is not accident
//                        break;

//                    iJumpNode++;
//                    if(iJumpNode > iMaxJumNode)
//                        break;
//                }
//            }
//            cout << "Phase 1 ranged index finished!" << endl;

//            vRangedIndex.push_back(iCurIndex);
//            iPreIndex = iCurIndex + 1;
//            cout << "Phase 2 ranged index finished!" << endl;

//            iJumpNode = 0;
//            for(unsigned int i = iCurIndex + 1; i<vDepthRange.size(); i++)
//            {
//                if( abs(vDepthRange[i].iDepth - iSmallestDepth) <= iMaxAllowed)
//                {
//                    vRangedIndex.push_back(i);
//                    iPreIndex = i;
//                    iJumpNode = 0;
//                }
//                else
//                {
//                    if(vDepthRange[i].iLen > iAccidentLen) // it is not accident
//                        break;

//                    iJumpNode++;
//                    if(iJumpNode > iMaxJumNode)
//                        break;
//                }
//            }
//            iPreIndex++;
//            cout << "Phase 3 ranged index finished!" << endl;

//            for(vector<int>::iterator itr = vRangedIndex.begin(); itr != vRangedIndex.end(); itr++)
//            {
//                iLowDepthLen += vDepthRange[*itr].iLen;
//                iTotalDepth += vDepthRange[*itr].iLen * vDepthRange[*itr].iDepth;
//            }
//        }

//        fAvgLowDepth = (float)iTotalDepth / iLowDepthLen;
//        cout << "fAvgLowDepth, TotalDepth and iLen: " << to_string(fAvgLowDepth) << ", "
//                                                      << to_string(iTotalDepth) << ", "
//                                                      << to_string(iLowDepthLen) << endl;

//        if(iLowDepthLen > iRangeLen * fRatio)
//        {
//            fTargetRangDepth = fAvgLowDepth;
//            iTargetRangLen = iLowDepthLen;
//            bLowDepthExtract = true;
//        }
//    }
//    //<--

//    if(!bLowDepthExtract)
//    {
//        //Collect maximum length
//        vector<int> vRangedIndex;
//        int iMaxAllowed = 2;

//        //This is for repetitive range
//        if(iDepth >= 2 * m_fAvgDepth)
//            iMaxAllowed = 3;

//        cout << "Do additional Depth Length Calculation" << endl;
//        cout << "\t" << to_string(vDepthRange.size()) << endl;
//        cout << "\t" << "Index: " << to_string(iIndex) << endl;

//        int iJumpNode = 0;
//        int iMaxJumNode = 1;
//        int iAccidentLen = 2;

//        for(int i = iIndex - 1; i >= 0; i--)
//        {
//            if( abs(vDepthRange[i].iDepth - iDepth) <= iMaxAllowed)
//            {
//                vRangedIndex.push_back(i);
//                iJumpNode = 0;
//            }
//            else
//            {
//                if(vDepthRange[i].iLen > iAccidentLen) // it is not accident
//                    break;

//                iJumpNode++;
//                if(iJumpNode > iMaxJumNode)
//                    break;
//            }
//        }
//        cout << "Phase 1 ranged index finished!" << endl;

//        vRangedIndex.push_back(iIndex);
//        cout << "Phase 2 ranged index finished!" << endl;

//        iJumpNode = 0;
//        for(unsigned int i = iIndex + 1; i<vDepthRange.size(); i++)
//        {
//            if( abs(vDepthRange[i].iDepth - iDepth) <= iMaxAllowed)
//            {
//                vRangedIndex.push_back(i);
//                iJumpNode = 0;
//            }
//            else
//            {
//                if(vDepthRange[i].iLen > iAccidentLen) // it is not accident
//                    break;

//                iJumpNode++;
//                if(iJumpNode > iMaxJumNode)
//                    break;
//            }
//        }
//        cout << "Phase 3 ranged index finished!" << endl;

//        int iLen = 0;
//        int iTotalDepth = 0;
//        for(vector<int>::iterator itr = vRangedIndex.begin(); itr != vRangedIndex.end(); itr++)
//        {
//            iLen += vDepthRange[*itr].iLen;
//            iTotalDepth += vDepthRange[*itr].iLen * vDepthRange[*itr].iDepth;
//            cout << "(" << to_string(vDepthRange[*itr].iDepth) << ", " << to_string(vDepthRange[*itr].iLen) << ")" << ", ";
//        }
//        cout << endl;

//        cout << "TotalDepth and iLen: " << to_string(iTotalDepth) << ", " << to_string(iLen) << endl;

//        float fAvgDepth = (float)iTotalDepth / iLen;

//        if( (iLen * .9 >= iLowDepthLen) &&
//            (fAvgDepth > m_fAvgDepth * 1.5 || iLowDepthLen < iRangeLen * .1))
//        {
//            fTargetRangDepth = fAvgDepth;
//            iTargetRangLen = iLen;
//            cout << "\t\t Get final result of target Depth (high) --- " << FloatToStr(fAvgDepth) << " --- " << to_string(iLen) << endl;
//        }
//        else
//        {
//            fTargetRangDepth = fAvgLowDepth;
//            iTargetRangLen = iLowDepthLen;
//            cout << "\t\t Get final result of target Depth (still do low) --- " << FloatToStr(fAvgLowDepth) << " --- " << to_string(iLowDepthLen) << endl;
//        }
//    }
}

void ClsSvDelDetect::GetDepthOfRange(string strBamFile,
                                     int iStartPos, int iEndPos, int iRangeLen,
                                     float& fTargetRangDepth, int& iTargetRangLen, string strChrom)
{
//    int iAvgDepth = 5;
//    int iOffSet = 50;
    if(strChrom == "")
        return;

    string strCmd = (string)"samtools depth -a -r " +
             strChrom + ":" + to_string(iStartPos) + "-" + to_string(iEndPos) +
             " " + strBamFile + " | awk '{print $3}'";

//    string strCmd = (string)"samtools depth -a -r " +
//             "11:" + to_string(stDRGroup.iStart - iOffSet) + "-" + to_string(stDRGroup.iEnd + iOffSet) +
//             " " + strBamFile + " | awk '{print $3}'";

    string strResult = exec(strCmd.c_str());
    vector<int> vCoverage;
    int iStart = 0;
    int iLen = 0;
    for(int i=0; i<strResult.length(); i++)
    {
        if(strResult.substr(i, 1) == "\n")
        {
            iLen = i - iStart;
            vCoverage.push_back(atoi(strResult.substr(iStart, iLen).c_str()));
            iStart = i + 1;
        }
    }

    //merge the depth with same value
    vector<St_DepthRange> vDepthRange;
    St_DepthRange stDepthRange;
    for(vector<int>::iterator itr = vCoverage.begin(); itr != vCoverage.end(); itr++)
    {
        stDepthRange.Clear();
        stDepthRange.iDepth = *itr;
        for(vector<int>::iterator subItr = itr; subItr != vCoverage.end(); subItr++)
        {
            if(stDepthRange.iDepth == *subItr)
            {
                stDepthRange.iLen++;
            }
            else
            {
                vDepthRange.push_back(stDepthRange);
                itr = subItr - 1;
                break;
            }

            if(subItr + 1 == vCoverage.end()) // the last one
            {
                vDepthRange.push_back(stDepthRange);
                itr = subItr;
                break;
            }
        }
    }

    for(vector<St_DepthRange>::iterator itr = vDepthRange.begin(); itr != vDepthRange.end(); itr++)
    {
        (*m_pOfs) << to_string(itr->iDepth) << endl;
        (*m_pOfs) << "\t\t" << to_string(itr->iLen) << endl;
    }
    (*m_pOfs) << to_string(vDepthRange.size()) << endl;

    //Try to find the largest and the longest one ->
    int iIndex = -1; // the index of the largest depth
    int iDepth = -1; // largest depth
    iLen = -1; //the length of the largest depth
    //<--

    int iZeroDepthLen = 0;

    //--> we also need to find the longest 0 range -->  Go!!!
    int iSmallestLongestIndex = -1;
    int iSmallestDepth = 10000;
    int iSmallestLongestLen = -1;
    vector<int> vZeroIndex;
    //<--
    for(unsigned int i=0; i<vDepthRange.size(); i++)
    {
        if(vDepthRange[i].iDepth == 0)
        {
            iZeroDepthLen += vDepthRange[i].iLen;
            vZeroIndex.push_back(i);
        }

        //Get the smallest and longest index -->
        if(vDepthRange[i].iDepth < iSmallestDepth)
        {
            iSmallestDepth = vDepthRange[i].iDepth;
            iSmallestLongestIndex = i;
            iSmallestLongestLen = vDepthRange[i].iLen;
        }
        else if(vDepthRange[i].iDepth == iSmallestDepth)
        {
            if(vDepthRange[i].iLen > iSmallestLongestLen)
            {
                iSmallestLongestIndex = i;
                iSmallestLongestLen = vDepthRange[i].iLen;
            }
        }
        //<--

        if(iDepth <= vDepthRange[i].iDepth) // update need to find the lagest one !! so change < to <=
        {
            bool bUpdate = true;
            if(iDepth == vDepthRange[i].iDepth &&
               vDepthRange[i].iLen < iLen)
                bUpdate = false;

            if(bUpdate)
            {
                iDepth = vDepthRange[i].iDepth;
                iLen = vDepthRange[i].iLen;
                iIndex = i;
            }
        }
    }
    if(iIndex < 0)
    {
        (*m_pOfs) << "Bad Depth Calculation" << endl;
        return;
    }

    (*m_pOfs) << "\t\t Get target Depth" << endl;
    (*m_pOfs) << "iZeroDepthLen: " << to_string(iZeroDepthLen) << endl;
    (*m_pOfs) << "iSmallestDepth: " << to_string(iSmallestDepth) << endl;
    (*m_pOfs) << "iSmallestLongestIndex: " << to_string(iSmallestLongestIndex) << endl;
    (*m_pOfs) << "iSmallestLongestLen: " << to_string(iSmallestLongestLen) << endl;

    float fRatio = 0.25;
    //int iRangeLen = stDRGroup.iEnd - stDRGroup.iStart + 1;
    //Check the range of 0 coverage: special treatment shuld be applied (may be discontinously) -->
    bool bLowDepthExtract = false;
    float fAvgLowDepth = -1;
    int iLowDepthLen = 0;
    if(iZeroDepthLen >= iRangeLen * fRatio)
    {
        fTargetRangDepth = 0;
        iTargetRangLen = iZeroDepthLen;

        (*m_pOfs) << "Set Zero Depth Length!" << endl;
        bLowDepthExtract = true;
    }
    else if(iSmallestDepth == 0 && iSmallestLongestIndex >= 0) //Get
    {
        int iPreIndex = 0;
        int iTotalDepth = 0;
        vector<int> vRangedIndex;
        for(vector<int>::iterator itr = vZeroIndex.begin(); itr != vZeroIndex.end(); itr++)
        {
            int iCurIndex = *itr;

            if(iCurIndex < iPreIndex)
                continue;

            vRangedIndex.clear();
            //We allow 1 diff
            int iMaxAllowed = 2; //change the coefficient from 1 to 2
            int iJumpNode = 0;
            int iMaxJumNode = 1;
            int iAccidentLen = 2;
            for(int i = iCurIndex - 1; i >= iPreIndex; i--)
            {
                if( abs(vDepthRange[i].iDepth - iSmallestDepth) <= iMaxAllowed)
                {
                    vRangedIndex.push_back(i);
                    iJumpNode = 0;
                }
                else
                {
                    if(vDepthRange[i].iLen > iAccidentLen) // it is not accident
                        break;

                    iJumpNode++;
                    if(iJumpNode > iMaxJumNode)
                        break;
                }
            }
            (*m_pOfs) << "Phase 1 ranged index finished!" << endl;

            vRangedIndex.push_back(iCurIndex);
            iPreIndex = iCurIndex + 1;
            (*m_pOfs) << "Phase 2 ranged index finished!" << endl;

            iJumpNode = 0;
            for(unsigned int i = iCurIndex + 1; i<vDepthRange.size(); i++)
            {
                if( abs(vDepthRange[i].iDepth - iSmallestDepth) <= iMaxAllowed)
                {
                    vRangedIndex.push_back(i);
                    iPreIndex = i;
                    iJumpNode = 0;
                }
                else
                {
                    if(vDepthRange[i].iLen > iAccidentLen) // it is not accident
                        break;

                    iJumpNode++;
                    if(iJumpNode > iMaxJumNode)
                        break;
                }
            }
            iPreIndex++;
            (*m_pOfs) << "Phase 3 ranged index finished!" << endl;

            for(vector<int>::iterator itr = vRangedIndex.begin(); itr != vRangedIndex.end(); itr++)
            {
                iLowDepthLen += vDepthRange[*itr].iLen;
                iTotalDepth += vDepthRange[*itr].iLen * vDepthRange[*itr].iDepth;
            }
        }

        fAvgLowDepth = (float)iTotalDepth / iLowDepthLen;
        (*m_pOfs) << "fAvgLowDepth, TotalDepth and iLen: " << to_string(fAvgLowDepth) << ", "
                                                      << to_string(iTotalDepth) << ", "
                                                      << to_string(iLowDepthLen) << endl;

        if(iLowDepthLen > iRangeLen * fRatio)
        {
            fTargetRangDepth = fAvgLowDepth;
            iTargetRangLen = iLowDepthLen;
            bLowDepthExtract = true;
        }
    }
    //<--

    if(!bLowDepthExtract)
    {
        //Collect maximum length
        vector<int> vRangedIndex;
        int iMaxAllowed = 2;

        //This is for repetitive range
        if(iDepth >= 2 * m_fAvgDepth)
            iMaxAllowed = 3;

        (*m_pOfs) << "Do additional Depth Length Calculation" << endl;
        (*m_pOfs) << "\t" << to_string(vDepthRange.size()) << endl;
        (*m_pOfs) << "\t" << "Index: " << to_string(iIndex) << endl;

        int iJumpNode = 0;
        int iMaxJumNode = 1;
        int iAccidentLen = 2;

        for(int i = iIndex - 1; i >= 0; i--)
        {
            if( abs(vDepthRange[i].iDepth - iDepth) <= iMaxAllowed)
            {
                vRangedIndex.push_back(i);
                iJumpNode = 0;
            }
            else
            {
                if(vDepthRange[i].iLen > iAccidentLen) // it is not accident
                    break;

                iJumpNode++;
                if(iJumpNode > iMaxJumNode)
                    break;
            }
        }
        (*m_pOfs) << "Phase 1 ranged index finished!" << endl;

        vRangedIndex.push_back(iIndex);
        (*m_pOfs) << "Phase 2 ranged index finished!" << endl;

        iJumpNode = 0;
        for(unsigned int i = iIndex + 1; i<vDepthRange.size(); i++)
        {
            if( abs(vDepthRange[i].iDepth - iDepth) <= iMaxAllowed)
            {
                vRangedIndex.push_back(i);
                iJumpNode = 0;
            }
            else
            {
                if(vDepthRange[i].iLen > iAccidentLen) // it is not accident
                    break;

                iJumpNode++;
                if(iJumpNode > iMaxJumNode)
                    break;
            }
        }
        (*m_pOfs) << "Phase 3 ranged index finished!" << endl;

        int iLen = 0;
        int iTotalDepth = 0;
        for(vector<int>::iterator itr = vRangedIndex.begin(); itr != vRangedIndex.end(); itr++)
        {
            iLen += vDepthRange[*itr].iLen;
            iTotalDepth += vDepthRange[*itr].iLen * vDepthRange[*itr].iDepth;
            (*m_pOfs) << "(" << to_string(vDepthRange[*itr].iDepth) << ", " << to_string(vDepthRange[*itr].iLen) << ")" << ", ";
        }
        (*m_pOfs) << endl;

        (*m_pOfs) << "TotalDepth and iLen: " << to_string(iTotalDepth) << ", " << to_string(iLen) << endl;

        float fAvgDepth = (float)iTotalDepth / iLen;

        if( (iLen * .9 >= iLowDepthLen) &&
            (fAvgDepth > m_fAvgDepth * 1.5 || iLowDepthLen < iRangeLen * .1))
        {
            fTargetRangDepth = fAvgDepth;
            iTargetRangLen = iLen;
            (*m_pOfs) << "\t\t Get final result of target Depth (high) --- " << FloatToStr(fAvgDepth) << " --- " << to_string(iLen) << endl;
        }
        else
        {
            fTargetRangDepth = fAvgLowDepth;
            iTargetRangLen = iLowDepthLen;
            (*m_pOfs) << "\t\t Get final result of target Depth (still do low) --- " << FloatToStr(fAvgLowDepth) << " --- " << to_string(iLowDepthLen) << endl;
        }
    }
}

void ClsSvDelDetect::GetClipFeatures(St_DRGroup& stDRGroup)
{
    //Get Close to Start Boundary and End Boundary respectively
    vector<St_BorderSCR> vSvStartClip;
    vector<St_BorderSCR> vSvEndClip;

    (*m_pOfs) << "stDRGroup.vDRRangedSCR Size: " << to_string(stDRGroup.vDRRangedSCR.size()) << endl;
    for(vector<St_BorderSCR>::iterator itr = stDRGroup.vDRRangedSCR.begin(); itr != stDRGroup.vDRRangedSCR.end(); itr++)
    {
        if(itr->enClipPart == cpRight)
            vSvStartClip.push_back(*itr);
        else if(itr->enClipPart == cpLeft)
            vSvEndClip.push_back(*itr);
    }

    //Sort --> Sv Start Clip
    int iOffSet = 50;
    int iNum = 1;
    int iMaxRangedNum = 0;
    if(vSvStartClip.empty() || vSvStartClip.size() == 1)
        stDRGroup.stGF.iStartClipReadsNum = vSvStartClip.size();
    else
    {
        (*m_pOfs) << "sort vSvStartClip" << endl;
        sort(vSvStartClip.begin(), vSvStartClip.end(), sort_scrClipPos_dup_func);
        (*m_pOfs) << "sort vSvStartClip <- end" << endl;

        for(vector<St_BorderSCR>::iterator itr = vSvStartClip.begin(); itr != vSvStartClip.end(); itr++)
        {
            if(itr + 1 == vSvStartClip.end())
                break;
            else
            {
                if(abs(itr->iClipPos - (itr+1)->iClipPos) <= iOffSet)
                    iNum++;
                else
                {
                    if(iNum > iMaxRangedNum)
                    {
                        iMaxRangedNum = iNum;
                        iNum = 1;
                    }
                }
            }
        }
        stDRGroup.stGF.iStartClipReadsNum = iMaxRangedNum;
    }

    //Sort --> Sv End Clip
    iNum = 1;
    iMaxRangedNum = 0;
    if(vSvEndClip.empty() || vSvEndClip.size() == 1)
        stDRGroup.stGF.iEndClipReadsNum = vSvEndClip.size();
    else
    {
        (*m_pOfs) << "sort vSvEndClip" << endl;
        sort(vSvEndClip.begin(), vSvEndClip.end(), sort_scrClipPos_dup_func);
        (*m_pOfs) << "sort vSvEndClip <- end" << endl;

        for(vector<St_BorderSCR>::iterator itr = vSvEndClip.begin(); itr != vSvEndClip.end(); itr++)
        {
            if(itr + 1 == vSvEndClip.end())
                break;
            else
            {
                if(abs(itr->iClipPos - (itr+1)->iClipPos) <= iOffSet)
                    iNum++;
                else
                {
                    if(iNum > iMaxRangedNum)
                    {
                        iMaxRangedNum = iNum;
                        iNum = 1;
                    }
                }
            }
        }
        stDRGroup.stGF.iEndClipReadsNum = iMaxRangedNum;
    }
}
