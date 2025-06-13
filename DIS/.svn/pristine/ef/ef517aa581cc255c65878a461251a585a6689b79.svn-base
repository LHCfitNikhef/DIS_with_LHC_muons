// C++ wrapper for powheginput and powheginputstring routines by Tomas Jezo, 2017

#include "powheginput_cpp.h"

double powheginput(std::string s) {
	return powheginput_(s.c_str(), s.size());	
}

void powheginputstring(std::string option, std::string &value) {
	// c_char will be a null terminated string, which we don't want, because we are passing stuff to fortran
	// first create a char* of the same length as option (not including the null character)
	int stringalen = option.length();
	char * stringa = new char[stringalen]; 
	// copy the memory from option.c_str() [which is of length option.lenthg()+1], everything but the null character
	memcpy(stringa, option.c_str(), option.length());
	int maxlen = 100;
	char * stringout = new char[maxlen];
	powheginputstring_(stringa, stringout, option.length(), maxlen);
	int i;
	for (i = 0; i < maxlen; i++)
		if (stringout[i] == ' ') break; 
	stringout[i] = '\0';
	
	// store it in string
//	value = "hxswg_ttbjets_part_v1bNLO_decay"; 
	value = stringout;
	delete[] stringa;
	delete[] stringout;
}
