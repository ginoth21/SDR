#include "dy4.h"
#include "stereo.h"
#include <vector>

void mix(std::vector<float> &out, const std::vector<float> &a, const std::vector<float> &b) {
    // all three arrays should be the same size
    int size = a.size();
    out.clear(); out.resize(size, 0.0);

    for (int i = 0; i < size; i++) {
        out[i] = a[i] * b[i];
    }
}

void mix_inplace(std::vector<float> &a, const std::vector<float> &b) {
    int size = a.size();
    for (int i = 0; i < size; i++) {
        a[i] = a[i] * b[i];
    }
}


void fast_allpass(std::vector<float> &out_block, std::vector<float> &in_block, std::vector<float> &state) {
    int block_size = in_block.size();
    int state_size = state.size(); // should be half standard filter size

    out_block.clear(); out_block.resize(block_size, 0.0);

    // copy state to front of out_block
    for (int i = 0; i < state_size; i++) {
        out_block[i] = state[i];
    }

    // copy in_block to out_block
    for (int i = 0; i < block_size - state_size; i++) {
        out_block[i + state_size] = in_block[i];
    }

    // copy end of in_block back into state for next block
    for (int i = 0; i < state_size; i++) {
        state[i] = in_block[block_size - state_size + i];
    }
}

void stereo_unmix(std::vector<float> &l, std::vector<float> &r, const std::vector<float> &mono, const std::vector<float> &stereo) {
    int size = mono.size();
    l.clear(); r.clear(); l.resize(size, 0.0); r.resize(size, 0.0);

    for (int i = 0; i < size; i++) {
        l[i] = 2*stereo[i] + mono[i];
        r[i] = 2*stereo[i] - mono[i];
    }
}

void stereo_to_interleaved(std::vector<float> &out, const std::vector<float> &mono, const std::vector<float> &stereo) {
    int size = mono.size();
    out.clear(); out.resize(size * 2);

    for (int i = 0; i < size; i++) {
        out[i*2]     = 2*stereo[i] + mono[i];
        out[i*2 + 1] = 2*stereo[i] - mono[i];
    }
}