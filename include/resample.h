#ifndef DY4_RESAMPLE_H
#define DY4_RESAMPLE_H

// add headers as needed
#include <vector>

// declaration of a function prototypes
void decimate(std::vector<float> &, std::vector<float> &, int);

void upsample(std::vector<float> &, int);

void resample(std::vector<float> &, std::vector<float> &, std::vector<float> &, int, int, std::vector<float> &state);
#endif // DY4_RESAMPLE_H
