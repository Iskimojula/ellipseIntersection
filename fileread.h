#ifndef FILEREAD_H
#define FILEREAD_H
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stack>
#include <string>
#include <map>
#include <vector>
#include <chrono>
#include <stack>
#include <stdio.h>
#include <stdint.h>

#include "algorithm/ellipse.h"

#ifdef _MSC_VER
#define STRTOK(aptr,delm,savP) strtok_s(aptr,delm,&savP)

#define TDEBUG(...)            qDebug(__VA_ARGS__)

#else
#define STRTOK(aptr,delm,savP) strtok_r(aptr,delm,&savP);

#define TDEBUG(...)            LOGD(__VA_ARGS__)
#endif

#define READBUFFSIZE 128
#define TESTDATAPATH "E:/code/qttest/pythontest/conic/test_data.txt"
#define RESDATAPATH "E:/code/qttest/pythontest/conic/res_datac_c.csv"
class CDataCenter
{
public:
    CDataCenter(){};
    ~CDataCenter(){};

public:
    void process(const std::string filepath);

};



#endif // FILEREAD_H
