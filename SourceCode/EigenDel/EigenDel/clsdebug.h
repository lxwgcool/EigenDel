#ifndef CLSDEBUG_H
#define CLSDEBUG_H

#include "clssvdeldetect.h"

struct St_SvSupportGroup
{
    vector<St_DRGroup> vDRGroup;
    vector<St_SCRGroup> vSCRGroup;
    St_SV stSv;

    St_SvSupportGroup()
    {}

    void Clear()
    {
        vDRGroup.clear();
        vSCRGroup.clear();
        stSv.Clear();
    }
};

class ClsDebug
{
public:
    ClsDebug();
    ~ClsDebug();

public:
    void GetStdSvDelRelatedGroup(vector<St_DRGroup>& vDRGroup, vector<St_SCRGroup>& vSCRGroup,
                                 vector<St_SV>& vSvDEL, string strChrom);
    void UpdateGroupBoundary();

public:
    void SetStrBamFile(string strV1);    
    void SetMeanInserSize(int iV1);
    void SetChromName(string strChromName);
    void SetOfstream(ofstream& ofs);

private:
//    void FilterBadDR(vector<St_DiscordantReads>& vDiscdRreads);

private:
    string m_strBamFile;
    int m_iMeanInsertSize;
    vector<St_SvSupportGroup> m_vSvSupportGroup;
    string m_strChromName;
    ofstream* m_pOfs;
};

#endif // CLSDEBUG_H
