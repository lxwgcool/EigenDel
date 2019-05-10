#ifndef CLSCONFIG_H
#define CLSCONFIG_H

#include <string>
using namespace std;

struct St_Config
{
    //General

    //Debug
    string strVcf; //VCF
    string strBamFile; //BAMFILE
    string strRef; // Reference File
    string strPythonCode;
    string strPicard;

    string strMultiAlignSeq;
    string strPindelPath;
    string strLumpyPath;
    string strDellyPath;
    string strSpritesPath;
    string strGASVproPath;
    string strCNVnatorPath;

    string strSampleName;
    string strChrom;


    St_Config():strVcf(""), strBamFile(""), strRef(""), strPythonCode(""), strPicard(""), strMultiAlignSeq(""),
                strPindelPath(""), strLumpyPath(""), strDellyPath(""), strSpritesPath(""),
                strGASVproPath(""), strCNVnatorPath(""), strSampleName(""), strChrom("")
    {}
};

class ClsConfig
{
public:
    ClsConfig();

public:
    void ReadConfig(St_Config& stConfig, const char *cpIniPath);
    bool CheckConfig(St_Config& stConfig);
};

#endif // CLSCONFIG_H
