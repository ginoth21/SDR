/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_STEREO_H
#define DY4_STEREO_H

// add headers as needed
#include <iostream>
#include <vector>

// declaration of a function prototypes
void mix(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void mix_inplace(std::vector<float> &, const std::vector<float> &);
void fast_allpass(std::vector<float> &, std::vector<float> &, std::vector<float> &);
void stereo_unmix(std::vector<float> &, std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void stereo_to_interleaved(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);

#endif // DY4_STEREO_H
