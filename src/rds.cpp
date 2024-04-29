#include "dy4.h"
#include "rds.h"
#include <vector>
#include <cmath>

void rrc(float Fs, unsigned short int num_taps, std::vector<float> &h){
    float T_symbol = 1/2375.0;
    float beta = 0.90; 
    float t = 0.0;

    h.clear(); h.resize(num_taps, 0.0);

    for (int k = 0; k < num_taps; k++){
        t = float(k - num_taps/2) / Fs;
        if (t == 0){
            h[k] = 1.0 + beta*((4/PI)-1);
        } else if ((t == -T_symbol/(4*beta)) || (t == T_symbol/(4*beta))){
            h[k] = (beta/sqrt(2))*(((1+2/PI) * (sin(PI/(4*beta)))) + ((1-2/PI)*(cos(PI/(4*beta)))));
        } else{
            h[k] = (sin(PI*t*(1-beta)/T_symbol) + 4*beta*(t/T_symbol)*cos(PI*t*(1+beta)/T_symbol)) / (PI*t*(1-(4*beta*t/T_symbol)*(4*beta*t/T_symbol)) / T_symbol);
        }
    }
}

// void differential_decode(std::vector<int> &rds_bitstream, const std::vector<int> &bits, int dif_state){
//     rds_bitstream.clear(); rds_bitstream.resize(bits.size());
//     rds_bitstream[0] = bits[0] ^ dif_state;
//     for(int i = 1; i < bits.size(); i++){
//         rds_bitstream[i] = bits[i] ^ bits[i-1];
//         // std::cerr << "BIT STREAM ARRAY: " << rds_bitstream[i] << "\n";
//     }
//     dif_state = bits[bits.size() - 1];
// }

void frame_sync(const std::vector<int> &bits){
    std::vector<int> window (26,0);  //sliding window of 26 bits

    //fill window with first 26 bits
    for (int i=0; i < 26; i++){
        window[i] = bits[i];
    }

    for (int i=26; i < bits.size(); i++){
        //if valid syndrome, break loop

        //else slide window by 1 bit
        window.erase(window.begin());
        window.push_back(bits[i]);
    }
}

void mat_mult(std::vector<int> &out, std::vector<int> &mat1, std::vector<std::vector<int>> &mat2){
    int cols = mat2[0].size();
    out.clear();
    out.resize(cols);
    for(int col = 0; col < cols; col++){
        for(int el = 0; el < mat2.size(); el++){
            out[col] = out[col] ^ (mat1[el] & mat2[el][col]);
        }
    }
}

unsigned short int get_pi(const std::vector<int> &frame) {
    unsigned short int id = 0;
    for (int i = 0; i < 16; i++) {
        id = id << 1;
        id += frame[i];
    }
    return id;
}

unsigned char get_pt(const std::vector<int> &frame) {
    unsigned char pt = 0;
    for (int i = 6; i < 11; i++) {
        pt = pt << 1;
        pt += frame[i];
    }
    return pt;
}

void print_pt(unsigned char pt) {
    const std::string lookup[32] = {
        "None",
        "News",
        "Information",
        "Sports",
        "Talk",
        "Rock",
        "Classic Rock",
        "Adult Hits",
        "Soft Rock",
        "Top 40",
        "Country",
        "Oldies",
        "Soft Music",
        "Nostalgia",
        "Jazz",
        "Classical",
        "Rhythm and Blues",
        "Soft Rhythm and Blues",
        "Language",
        "Religious Music",
        "Religious Talk",
        "Personality",
        "Public",
        "College",
        "Spanish Talk",
        "Spanish Music",
        "Hip Hop",
        "Unassigned",
        "Unassigned",
        "Weather",
        "Emergency Test",
        "Emergency"
    };
    std::cerr << lookup[pt] << std::endl;
}