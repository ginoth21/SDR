/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>

// declaration of a function prototypes
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);
void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &, int);
void singleconvdown(std::vector<float> &, const std::vector<float> &, const std::vector<float> &, int);
void blockconvdown(std::vector<float> &, const std::vector<float> &, const std::vector<float> &, std::vector<float> &, int);
void blockconvdowni(std::vector<float> &, const std::vector<float> &, const std::vector<float> &, std::vector<float> &, int);
void blockconvdownq(std::vector<float> &, const std::vector<float> &, const std::vector<float> &, std::vector<float> &, int);
void blockconvdownsplit(std::vector<float> &, std::vector<float> &, const std::vector<float> &, const std::vector<float> &, std::vector<float> &, int);
void impulseResponseBPF(float, float, float, unsigned short int, std::vector<float> &);
void impulseResponseBPF_algorithm(float, float, float, unsigned short int, std::vector<float> &);
void delayBlock(std::vector<float> &, std::vector<float> &);

#endif // DY4_FILTER_H
