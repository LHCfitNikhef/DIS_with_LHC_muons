// C++ wrapper for powheginput and powheginputstring routines by Tomas Jezo, 2017
#pragma once

#include <sstream>
#include <cstring>

// these are implemented in powheginput.f
extern "C" double powheginput_(const char*, int);
extern "C" void powheginputstring_(char*, char*, int, int);

// these are implemented in powheginput_cpp.cc
double powheginput(std::string s);
void powheginputstring(std::string option, std::string &value);
