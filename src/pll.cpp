#include "dy4.h"
#include "pll.h"
#include <cmath>

void fmPll(std::vector<float> &out, const std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth, PllState s) {
    // DEFAULTS
    // ncoScale = 1.0
    // phaseAdjust = 0.0
    // normBandwidth = 0.01
    
    int size = pllIn.size();
    out.clear(); out.resize(size, 0.0);

    float Cp = 2.666;
	float Ci = 3.555;

    float Kp = (normBandwidth) * Cp;
    float Ki = (normBandwidth * normBandwidth) * Ci;

    if (!s.valid) {
        s.integrator = 0.0;
        s.phaseEst = 0.0;
        s.feedbackI = 1.0;
        s.feedbackQ = 0.0;
        s.firstVal = 1.0;
        s.trigOffset = 0;
        s.valid = true;
    }

	out[0] = s.firstVal;

    for (int k = 0; k < size; k++) {
        float errorI = pllIn[k] * (s.feedbackI);
		float errorQ = pllIn[k] * (-s.feedbackQ);

		float errorD = atan2(errorQ, errorI);

		s.integrator = s.integrator + Ki*errorD;

		s.phaseEst = s.phaseEst + Kp*errorD + s.integrator;

		s.trigOffset += 1;
		float trigArg = 2*PI*(freq/Fs)*(s.trigOffset) + s.phaseEst;
		s.feedbackI = cos(trigArg);
		s.feedbackQ = sin(trigArg);
        float outval = cos(trigArg*ncoScale + phaseAdjust);
        if (k == size - 1) {
            s.firstVal = outval;
        } else {
		    out[k+1] = outval;
        }
    }
}

void fmPllIQ(std::vector<float> &out, std::vector<float> &outQ, const std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth, PllState s) {
    // DEFAULTS
    // ncoScale = 1.0
    // phaseAdjust = 0.0
    // normBandwidth = 0.01
    
    int size = pllIn.size();
    out.clear(); out.resize(size, 0.0);
    outQ.clear(); outQ.resize(size, 0.0);

    float Cp = 2.666;
	float Ci = 3.555;

    float Kp = (normBandwidth) * Cp;
    float Ki = (normBandwidth * normBandwidth) * Ci;

    if (!s.valid) {
        s.integrator = 0.0;
        s.phaseEst = 0.0;
        s.feedbackI = 1.0;
        s.feedbackQ = 0.0;
        s.firstVal = 1.0;
        s.firstValQ = 0.0;
        s.trigOffset = 0;
        s.valid = true;
    }

	out[0] = s.firstVal;
    outQ[0] = s.firstVal;

    for (int k = 0; k < size; k++) {
        float errorI = pllIn[k] * (s.feedbackI);
		float errorQ = pllIn[k] * (-s.feedbackQ);

		float errorD = atan2(errorQ, errorI);

		s.integrator = s.integrator + Ki*errorD;

		s.phaseEst = s.phaseEst + Kp*errorD + s.integrator;

		s.trigOffset += 1;
		float trigArg = 2*PI*(freq/Fs)*(s.trigOffset) + s.phaseEst;
		s.feedbackI = cos(trigArg);
		s.feedbackQ = sin(trigArg);
        float outval = cos(trigArg*ncoScale + phaseAdjust);
        float outvalQ = sin(trigArg*ncoScale + phaseAdjust);
        if (k == size - 1) {
            s.firstVal = outval;
            s.firstValQ = outvalQ;
        } else {
		    out[k+1] = outval;
            outQ[k+1] = outvalQ;
        }
    }
}