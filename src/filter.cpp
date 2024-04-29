/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <cmath>

// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// allocate memory for the impulse response
	h.clear(); h.resize(num_taps, 0.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	if (Fc > Fs/2) {
		std::cerr << "Cutoff must be smaller than sample rate / 2\n";
		exit(-1);
	}
	float tapsovertwo = (num_taps - 1)/2;
	float normcutoff = Fc / (Fs / 2);
	for (int i = 0; i < num_taps; i++) {
		if (i == tapsovertwo) {
			h[i] = normcutoff;
		} else {
			h[i] = normcutoff*std::sin(PI * normcutoff * (float)(i-tapsovertwo)) / (PI * normcutoff * (float)(i-tapsovertwo));
		}
		h[i] *= std::sin((float)i * PI / (float)num_taps) * std::sin((float)i * PI / (float)num_taps);
	}
}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	// allocate memory for the output (filtered) data
	y.clear(); y.resize(x.size()+h.size()-1, 0.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	for (int n = 0; n < y.size(); n++) {
		for (int k = 0; k < h.size(); k++) {
			if (n-k >= 0 && n-k < x.size()) {
				y[n] += h[k] * x[n-k];
			} else {
				y[n] += 0;
			}
		}
	}
}

void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, int block_size)
{
	// allocate memory for the output (filtered) data
	y.clear(); y.resize(x.size()+h.size()-1, 0.0);

	int num_blocks = x.size() / block_size;
	std::vector<float> xb, yb, state;
	xb.resize(block_size, 0.0);
	yb.resize(block_size, 0.0);
	state.resize(h.size() - 1, 0.0);

	for (int i = 0; i < num_blocks; i++){
		//extract block to be processed
		for (int j = 0; j < block_size; j++){
			xb[j] = x[i*block_size + j];
		}
		
		// the rest of the code in this function is to be completed by you
		// based on your understanding and the Python code from the first lab
		for (int n = 0; n < yb.size(); n++) {
			for (int k = 0; k < h.size(); k++) {
				if (n-k >= 0) {
					yb[n] += h[k] * xb[n-k];
				} else {
					// use value from the previous state 
					yb[n] += h[k] * state[(state.size() - std::abs(n-k))];
				}
			}
		}

		//store block output into final output vector
		for (int a = 0; a < block_size; a++){
			y[i*block_size + a] = yb[a];
		}

		//save state for next block
		int state_index = 0;
		for (int b = (block_size - state.size()); b < block_size; b++){
			state[state_index] = xb[b];
			state_index++;
		}
	}
}

void singleconvdown(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h , int d){
	// allocate memory for the output (filtered) data
	y.clear(); y.resize(x.size()/d, 0.0);
	
	for (int i=0 ; i<(x.size()/d) ; i++){
		for(int j=0 ; j<h.size(); j++){
			int index=i*d-j;
			if(index>=0 && index<x.size() ){
				y[i] += h[j] * x[index];
			}
		}
	}
	
}

void blockconvdown(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state, int d) {
	int x_size = x.size();
	int h_size = h.size();
	int state_size = state.size();
	
	// clear output vector
	y.clear(); y.resize(x_size/d, 0.0);
	int y_size = y.size();

	// fast convolution/downsampling in one
	for (int i = 0; i < y_size; i++) {
		for (int j = 0; j < h_size; j++) {
			int index = i*d - j;
			if (index >= 0 && index < x_size) {
				y[i] += h[j] * x[index];
			} else if (index < 0) {
				y[i] += h[j] * state[state_size + index];
			}
		}
	}

	// clear state vector (comment for improved performance?)
	//state.clear(); state.resize(h.size() - 1, 0.0);

	// copy end of input into state vector for next block
	int state_index = 0;
	for (int b = (x_size - state_size); b < x_size; b++){
		state[state_index] = x[b];
		state_index++;
	}
}

void blockconvdowni(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state, int d) {
	int x_size = x.size();
	int h_size = h.size();
	int state_size = state.size();
	
	// clear output vector
	y.clear(); y.resize(x_size/(d*2), 0.0);
	int y_size = y.size();

	// fast convolution/downsampling in one
	for (int i = 0; i < y_size; i++) {
		for (int j = 0; j < h_size; j++) {
			int index = 2*i*d - 2*j;
			if (index >= 0 && index < x_size) {
				y[i] += h[j] * x[index];
			} else if (index < 0) {
				y[i] += h[j] * state[state_size + index];
			}
		}
	}

	// clear state vector (comment for improved performance?)
	//state.clear(); state.resize(h.size() - 1, 0.0);

	// copy end of input into state vector for next block
	int state_index = 0;
	for (int b = (x_size - state_size); b < x_size; b++){
		state[state_index] = x[b];
		state_index++;
	}
}

void blockconvdownq(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state, int d) {
	int x_size = x.size();
	int h_size = h.size();
	int state_size = state.size();
	
	// clear output vector
	y.clear(); y.resize(x_size/(d*2), 0.0);
	int y_size = y.size();

	// fast convolution/downsampling in one
	for (int i = 0; i < y_size; i++) {
		for (int j = 0; j < h_size; j++) {
			int index = 2*i*d - 2*j+1;
			if (index >= 0 && index < x_size) {
				y[i] += h[j] * x[index];
			} else if (index < 0) {
				y[i] += h[j] * state[state_size + index];
			}
		}
	}

	// clear state vector (comment for improved performance?)
	//state.clear(); state.resize(h.size() - 1, 0.0);

	// copy end of input into state vector for next block
	int state_index = 0;
	for (int b = (x_size - state_size); b < x_size; b++){
		state[state_index] = x[b];
		state_index++;
	}
}

void blockconvdownsplit(std::vector<float> &yi, std::vector<float> &yq, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state, int d) {
	int x_size = x.size();
	int h_size = h.size();
	int state_size = state.size();

	// clear output vectors
	yi.clear(); yi.resize(x_size/(d*2), 0.0); yq.clear(); yq.resize(x_size/(d*2), 0.0);
	int y_size = yi.size();

	// fast convolution/downsampling in one
	for (int i = 0; i < y_size; i++) {
		for (int j = 0; j < h_size; j++) {
			int index = 2*(i*d - j);
			if (index >= 0 && index + 1 < x_size) {
				yi[i] += h[j] * x[index];
				yq[i] += h[j] * x[index + 1];
			} else if (index < 0) {
				yi[i] += h[j] * state[state_size + index];
				yq[i] += h[j] * state[state_size + index + 1];
			}
		}
	}

	// clear state vector (comment for improved performance?)
	//state.clear(); state.resize(h.size() - 1, 0.0);

	// copy end of input into state vector for next block
	int state_index = 0;
	for (int b = (x_size - state_size); b < x_size; b++){
		state[state_index] = x[b];
		state_index++;
	}
}

void impulseResponseBPF(float Fs, float Fc_upper, float Fc_lower, unsigned short int num_taps, std::vector<float> &h) {
	std::vector<float> wideLPF;
	std::vector<float> narrowLPF;

	impulseResponseLPF(Fs, Fc_upper, num_taps, wideLPF);
	impulseResponseLPF(Fs, Fc_lower, num_taps, narrowLPF);

	h.clear(); h.resize(wideLPF.size(), 0.0);
	for (int i = 0; i < h.size(); i++) {
		h[i] = wideLPF[i] - narrowLPF[i];
	}
}

void impulseResponseBPF_algorithm(float Fb, float Fe, float Fs, unsigned short int num_taps, std::vector<float> &h) {	//change name of function later
	// allocate memory for the impulse response
	h.clear(); h.resize(num_taps, 0.0);

	float normcenter = ((Fb + Fe)/2) / (Fs/2);	//can remove /2's if it helps speed up
	float normpass = (Fe - Fb) / (Fs/2);
	float tapsovertwo = (num_taps - 1)/2;

	for (int i = 0; i < num_taps; i++) {
		if (i == tapsovertwo) {
			h[i] = normpass;
		} else {
			h[i] = normpass*std::sin(PI * (float)(normpass/2) * (float)(i-tapsovertwo)) / (PI * (float)(normpass/2) * (float)(i-tapsovertwo));
		}
		h[i] *= std::cos((float)(i - tapsovertwo) * PI * (float)normcenter);
		h[i] *= std::pow(std::sin((float)(i * PI/num_taps)), 2);
	}
}

void delayBlock(std::vector<float> &block, std::vector<float> &state){
    // Iterates backwards through the state vector
    for(int n = 1; n <= state.size(); n++){
        // Inserts element from back of state vector at beginning of block vector
        block.insert(block.begin(),state[state.size() - n]);

        // Replace element of state vector that was just inserted with last element of block vector (then remove that element)
        state[state.size() - n] = block[block.size() - 1];
        block.pop_back();
    }

    block.shrink_to_fit();
}
