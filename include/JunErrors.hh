#ifndef JunErrors_HH
#define JunErrors_HH 1

#include <iostream>
#include <string>

#include <TString.h>

// printed color
#define BLACK "\033[30m"
#define RED "\033[31m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define BLUE "\033[34m"
#define PURPLE "\033[35m"
#define CYAN "\033[36m"
#define WHITE "\033[37m"
#define BOLD "\033[1m"
#define COLOR_END "\033[0m"
using namespace std;

void MiaoError(string errorStr="");
void MiaoWarn(string warnStr="");
void MiaoHint(string hintStr="");
#endif