/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_PLL_H
#define DY4_PLL_H

// add headers as needed
#include <vector>

class PllState {
    public:
    bool valid;
    float integrator;
	float phaseEst;
	float feedbackI;
	float feedbackQ;
	float firstVal;
	float firstValQ;
	int trigOffset;
};

// declaration of a function prototypes
void fmPll(std::vector<float> &, const std::vector<float> &, float, float, float, float, float, PllState state);
void fmPllIQ(std::vector<float> &, std::vector<float> &, const std::vector<float> &, float, float, float, float, float, PllState state);

#endif // DY4_PLL_H