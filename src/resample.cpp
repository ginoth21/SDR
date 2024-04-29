#include "dy4.h"
#include "resample.h"
#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>

void decimate(std::vector<float> &x, std::vector<float> &y, int factor) {
    int new_size = x.size() / factor;
    y.clear(); y.resize(new_size);
    for (int n = 0; n < new_size; n ++) {
        y[n] = x[n*factor];
    }
}

void upsample(std::vector<float> &y, int factor) {
    int new_size=y.size()*factor;
    std::vector<float> ytemp = y;
    y.clear(); y.resize(new_size); //new vector of upscaled size with all zeros
    
    for (int n=0; n<ytemp.size(); n++){
        y[n*factor]=factor*ytemp[n]; //put input in every "factor" spots of upscaled zero vector
    }
    
    
}

void resample(std::vector<float> &y, std::vector<float> &x, std::vector<float> &h, int U, int D, std::vector<float> &state){ //add state saving!!!!!!!!!!!!
    // Clock stuff for testing
    //~ int seed = std::time(0x0);
    //~ auto start_time = std::chrono::high_resolution_clock::now();
    
    // Define these so we're not calling ~.size() on repeat
    int h_size = static_cast<int>(h.size());
    int x_size = static_cast<int>(x.size());
    int state_size = static_cast<int>(state.size());
    
    float y_long_size = (x_size * U) - 1;
    int y_size = round((y_long_size / D) + 0.5);
    y.clear(); y.resize(y_size);
   // int state_index=0;
    
    for(int n = 0; n < y_long_size; n += D){
        for(int k = n % U; k < h_size; k += U){
            // std::cerr << "n: " << n << ", k: " << k << std::endl;
            if((n-k >= 0) && ((n-k) < (x_size*U))){
                y[int(n/D)] += h[k] * x[int((n-k)/U)] * U; //U is to increase the gain back 
            } 
            else if ((n-k)/U < 0){
                y[int(n/D)] += h[k] * state[state_size + int((n-k))] * U;
            }
        }
        //~ if ( n>(y_long_size-state_size) ){ //state saving vector updates
            //~ std::cerr<<"x index is"<< state_size << "\n";
            //~ state[state_index] = x[n/U];
            //~ state_index++;
            //~ }
    }
    int state_index = 0;
	for (int b = (y_long_size - state_size); b < y_long_size; b++){ //state saving vector updates
    //  std::cerr << "y long size = " << y_long_size << ", state size = " << state_size << ", h size = " << h_size << "b="<< b/U << "\n";
      // if (b>0){
           state[state_index] = x[b/U];
       //}
		state_index++;
	}

    // for (int b = 0; b < state_size; b++){
    //     state[b] = x[x_size - state_size - b];
    // }
    // More clock stuff for testing
    //~ auto stop_time = std::chrono::high_resolution_clock::now();
    //~ std::chrono::duration<double, std::milli> Test_run_time = stop_time - start_time;
    //~ std::cerr<<"the run time was "<< Test_run_time.count() << " milliseconds \n";
}
