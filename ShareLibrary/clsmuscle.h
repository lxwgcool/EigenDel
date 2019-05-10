#ifndef CLSMUSCLE_H
#define CLSMUSCLE_H
#include <string>
using namespace std;

/* The target of this class is trying to do multiple alignment
 * Function
 * (1) Do multiple alignment
 * (2) Merge the aligned sequence together
 */

class ClsMuscle
{
public:
    ClsMuscle();
    ~ClsMuscle();
public:
    string Run(string strFaPath);
    string MergeAlignedSeq(string strSeqPath);

private:
    string m_strMusclePath;
};

#endif // CLSMUSCLE_H
