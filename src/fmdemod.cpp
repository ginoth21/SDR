#include "dy4.h"
#include "fmdemod.h"
#include <vector>

void FM_Demod(std::vector<float> &i, const std::vector<float> &q, std::vector<float> &x, float* previ, float* prevq) {
    x.clear(); x.resize(i.size(), 0.0);

    for (int k = 0; k < i.size(); k++) {
        x[k] = ((i[k] * (q[k] - *prevq)) - (q[k] * (i[k] - *previ)))/(i[k] * i[k] + q[k] * q[k]);
        *previ = i[k];
        *prevq = q[k];
    }
}
