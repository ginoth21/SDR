/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "fmdemod.h"
#include "resample.h"
#include "pll.h"
#include "stereo.h"
#include "rds.h"
#include <chrono>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>

//#define DEBUG
#define QUEUE_ELEMS 10

void audio_function(std::queue<std::vector<float>> &queue, std::mutex& mut, std::condition_variable& cvar,
	char mode, char path, int block_size, int D, int U, int num_taps, int num_taps_final,
	float if_fs, float if_fs_final, float audio_fs, float audio_fc, float audio_fc_final, bool* run) {

	std::vector<float> mono_lpf;
	impulseResponseLPF(if_fs_final, audio_fc_final, num_taps_final, mono_lpf);

	std::vector<float> mono_down_state;
	mono_down_state.resize(mono_lpf.size() - 1, 0.0);

	std::vector<float> mono_audio;

	std::vector<float> pilot_bpf, stereo_bpf;
	impulseResponseBPF(if_fs, 19500.0, 18500.0, num_taps, pilot_bpf);
	impulseResponseBPF(if_fs, 54000.0, 22000.0, num_taps, stereo_bpf);

	std::vector<float> pilot_if, stereo_if, nco_if, mixed_if;
	std::vector<float> pilot_if_state, stereo_if_state;
	pilot_if_state.resize(pilot_bpf.size() - 1, 0.0);
	stereo_if_state.resize(stereo_bpf.size() - 1, 0.0);
	PllState ps;
	ps.valid = false;
	std::vector<float> mono_delay_if, mono_delay_if_state;
	mono_delay_if_state.resize((stereo_bpf.size() - 1)/2, 0.0);
	std::vector<float> stereo_audio, stereo_down_state;
	stereo_down_state.resize(mono_down_state.size(), 0.0);
	std::vector<float> l, r, interleaved_audio;

	// std::cerr << "RUNPOINTER aduio: " << *run << "\n";
	while (*run || !queue.empty()) {
		std::unique_lock<std::mutex> lock(mut);
		while (queue.empty()) {
			cvar.wait(lock);
		}

		std::vector<float> fmout = queue.front();
		// std::cerr << "AUDIO FMOUT SIZE: " << fmout.size() << "\n";
		queue.pop();
		cvar.notify_one();
		lock.unlock();

		switch (path) {
			case 'm':{
				if (mode=='0' | mode=='1'){
					blockconvdown(mono_audio, fmout, mono_lpf, mono_down_state, D);
				}
				else{
					//std::cerr << "start resample\n";
					resample(mono_audio,fmout,mono_lpf,U,D,mono_down_state);
					//std::cerr << "final audio filter done\n";
				}	
				// write audio to stdout
				writeStdoutBlockData(mono_audio);
				break;
			}
			case 's': case 'r': {
				auto time1 = std::chrono::high_resolution_clock::now();
				blockconvdown(pilot_if, fmout, pilot_bpf, pilot_if_state, 1);
				blockconvdown(stereo_if, fmout, stereo_bpf, stereo_if_state, 1);
				fmPll(nco_if, pilot_if, 19000.0, if_fs, 2.0, 0.0, 0.01, ps);
				mix(mixed_if, nco_if, stereo_if);
				fast_allpass(mono_delay_if, fmout, mono_delay_if_state);
				
				if (mode=='0' | mode=='1'){
					blockconvdown(mono_audio, mono_delay_if, mono_lpf, mono_down_state, D);
					blockconvdown(stereo_audio, mixed_if, mono_lpf, stereo_down_state, D);
				}
				else{
					resample(mono_audio,mono_delay_if,mono_lpf,U,D,mono_down_state);
					resample(stereo_audio,mixed_if,mono_lpf,U,D,stereo_down_state);
				}
				
				auto time2 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double, std::milli> time_dif = time2 - time1;
				//std::cerr << "case s (- interleaving and writing) took " << time_dif.count() << " ms" << std::endl;
				stereo_to_interleaved(interleaved_audio, mono_audio, stereo_audio);
				writeStdoutBlockData(interleaved_audio);
				break;
			}
		}
	}
}

void rds_function(std::queue<std::vector<float>> &queue, std::mutex& mut, std::condition_variable& cvar,
	char mode, int block_size, int num_taps, float if_fs, bool* run) {

	int rds_sps = ((mode == '2')? 16 : 35);
	float rds_fs = 2375.0 * rds_sps;
	int u = (mode == '2')? 19 : 133;
	int d = (mode == '2')? 120 : 384;
	
	// P Matrix for frame synch
	std::vector<std::vector<int>> P = {{1,0,0,0,0,0,0,0,0,0},
									   {0,1,0,0,0,0,0,0,0,0},
								       {0,0,1,0,0,0,0,0,0,0},
								       {0,0,0,1,0,0,0,0,0,0},
								       {0,0,0,0,1,0,0,0,0,0},
								       {0,0,0,0,0,1,0,0,0,0},
								       {0,0,0,0,0,0,1,0,0,0},
								       {0,0,0,0,0,0,0,1,0,0},
								       {0,0,0,0,0,0,0,0,1,0},
								       {0,0,0,0,0,0,0,0,0,1},
								       {1,0,1,1,0,1,1,1,0,0},
								       {0,1,0,1,1,0,1,1,1,0},
								       {0,0,1,0,1,1,0,1,1,1},
								       {1,0,1,0,0,0,0,1,1,1},
								       {1,1,1,0,0,1,1,1,1,1},
								       {1,1,0,0,0,1,0,0,1,1},
								       {1,1,0,1,0,1,0,1,0,1},
								       {1,1,0,1,1,1,0,1,1,0},
								       {0,1,1,0,1,1,1,0,1,1},
								       {1,0,0,0,0,0,0,0,0,1},
								       {1,1,1,1,0,1,1,1,0,0},
								       {0,1,1,1,1,0,1,1,1,0},
								       {0,0,1,1,1,1,0,1,1,1},
								       {1,0,1,0,1,0,0,1,1,1},
								       {1,1,1,0,0,0,1,1,1,1},
								       {1,1,0,0,0,1,1,0,1,1}};

	// Syndrome vectors
    std::vector<int> A = {1,1,1,1,0,1,1,0,0,0};
    std::vector<int> B = {1,1,1,1,0,1,0,1,0,0};
    std::vector<int> C = {1,0,0,1,0,1,1,1,0,0};
    std::vector<int> C_prime = {1,1,1,1,0,0,1,1,0,0};
    std::vector<int> D = {1,0,0,1,0,1,1,0,0,0};

	PllState rds_ps;
	
	std::vector<float> rds_bpf, rds_channel_state, rrcf, rds_cr_bpf, rds_cr_state, rds_cr_raw, rds_cr_filt, ncoI, ncoQ, rds_demod_mix, rds_lpf, rds_resample, rds_resample_state, rds_cos_filt, rds_cos_state, rds_channel_delay, rds_channel_delay_state;
	impulseResponseBPF(if_fs, 60000.0, 54000.0, num_taps, rds_bpf);
	rds_channel_state.resize(rds_bpf.size() - 1, 0.0);
	impulseResponseBPF(if_fs, 114500.0, 113500.0, num_taps, rds_cr_bpf);
	rds_cr_state.resize(rds_cr_bpf.size() - 1, 0.0);
	impulseResponseLPF(if_fs, 3000.0, num_taps, rds_lpf);
	rrc(rds_fs, num_taps, rrcf);
	rds_cos_state.resize(rrcf.size() - 1, 0.0);
	rds_resample_state.resize(rds_lpf.size() - 1);

	std::vector<float> rds_channel;
	
	// State for manchester decoding
	int man_state = -1;
	// State for differential decoding
	int dif_state = 0;
	// State for frame
	std::vector<int> frame_state = {};

	// PI and PT
	unsigned short int saved_PI;
	unsigned char saved_PT;

	int blk_count = 0;

	while (*run || !queue.empty()) {

		std::unique_lock<std::mutex> lock(mut);
		while (queue.empty()) {
			cvar.wait(lock);
		}

		// wait for queue, signals, etc.

		std::vector<float> fmout = queue.front();
		queue.pop();
		cvar.notify_one();
		lock.unlock();

		blockconvdown(rds_channel, fmout, rds_bpf, rds_channel_state, 1);
		fast_allpass(rds_channel_delay, rds_channel, rds_channel_delay_state);
		mix(rds_cr_raw, rds_channel, rds_channel);
		blockconvdown(rds_cr_filt, rds_cr_raw, rds_cr_bpf, rds_cr_state, 1);
		fmPllIQ(ncoI, ncoQ, rds_cr_filt, 114000.0, if_fs, 0.5, 0.0, 0.0025, rds_ps);
		mix(rds_demod_mix, rds_channel_delay, ncoI);
		resample(rds_resample, rds_demod_mix, rds_lpf, u, d, rds_resample_state);
		blockconvdown(rds_cos_filt, rds_resample, rrcf, rds_cos_state, 1);
		
		//std::cerr << "Block count: " << blk_count << "\n";
		blk_count +=  1;
		//data recovery here
		// find first zero crossing
		int zero_cross_index = 0;
		for(int i = 1; i < rds_cos_filt.size(); i++){
			if((rds_cos_filt[i-1] > 0 && rds_cos_filt[i] < 0) || (rds_cos_filt[i-1] < 0 && rds_cos_filt[i] > 0)){
				zero_cross_index = i;
				break;
			}
		}
		
		// Find first peak
		int first_peak_index = zero_cross_index;
		if(zero_cross_index >= (rds_sps + (rds_sps/2))){
			first_peak_index -= (rds_sps + (rds_sps / 2));
		}else if(zero_cross_index >= (rds_sps/2)) {
			first_peak_index -= (rds_sps / 2);
		} else{
			first_peak_index += (rds_sps/2);
		}
		
		// Find remaining peaks 
		std::vector<int> peaks;
		peaks.resize(int(rds_cos_filt.size() / rds_sps));

		for(int i = 0; i < peaks.size(); i++){
			if(rds_cos_filt[first_peak_index + i * rds_sps] > 0){
				peaks[i] = 1;
			}
			else{
				peaks[i] = 0;
			}
		}
		
		// Append state to beginning of peaks if applicable
		if(man_state != -1){
			peaks.insert(peaks.begin(),man_state);
		}
		
		// Loop through and find the allignment with fewer HH or LL
		
		std::vector<int> bits1, bits2, mybits;
		//mybits.resize(int(peaks.size() / 2), 0);
		bits1.resize(int(peaks.size() / 2), 0);
		bits2.resize(int(peaks.size() / 2), 0);
		int mistakes1 = 0;
		int mistakes2 = 0;
		for(int i = 1; i < peaks.size(); i+=2){
			if(peaks[i-1] == peaks[i]){
				mistakes1 += 1;
				bits1[i/2] = 0;
			}
			else{
				bits1[i/2] = peaks[i-1];
			}
			if((i + 1) < peaks.size()){
				if(peaks[i+1] == peaks[i]){
					mistakes2 += 1;
					bits2[i/2] = 0;
				}
				else{
					bits2[i/2] = peaks[i];
				}
			}
		}
		
		// Throw out the one with more mistakes to get Manchester decoded
		if(mistakes1 >= mistakes2){
			//std::vector<int> bits(bits2);
			//vector_copy(bits2, mybits);
			mybits = bits2;
			// Save the state if even number of peaks
			man_state = ((peaks.size() % 2) == 0) ? peaks[peaks.size() - 1] : -1;
		}
		else{
			//vector_copy(bits1, mybits);
			mybits = bits1;
			// Save the state if odd number of peaks
			man_state = ((peaks.size() % 2) == 0) ? -1 : peaks[peaks.size() - 1];
		}

		
		// Differential decoding
		std::vector<int> rds_bitstream;
		// differential_decode(rds_bitstream, bits, dif_state);
		rds_bitstream.clear();
		rds_bitstream.resize(mybits.size());
		rds_bitstream[0] = mybits[0] ^ dif_state;
		for(int i = 1; i < mybits.size(); i++){
			rds_bitstream[i] = mybits[i] ^ mybits[i-1];
			// std::cerr << "BIT STREAM ARRAY: " << rds_bitstream[i] << "\n";
		}
		dif_state = mybits[mybits.size() - 1];
		
		// Frame synchronization
		// Add state to beginning of bitstream
		if(!frame_state.empty()){
			rds_bitstream.insert(rds_bitstream.begin(), frame_state.begin(),frame_state.end());
		}
        int frame_index = 0;
        std::vector<int> frame;
        frame.clear();
        frame.resize(26);
        std::vector<int> mat_mult_result;
        bool synched = false;
        char synch_type = 'Z';
        while((synched == false) && ((frame_index + 26) < rds_bitstream.size())){
			// std::cerr << frame_index << std::endl;
            frame = std::vector<int>(rds_bitstream.begin() + frame_index, rds_bitstream.begin() + frame_index + 26);
            //printRealVector(frame);
			mat_mult(mat_mult_result,frame,P);
			// printRealVector(mat_mult_result);
            if(mat_mult_result == A){
                synch_type = 'A';
                synched = true;
            }else if(mat_mult_result == B){
                synch_type = 'B';
                synched = true;
            }else if(mat_mult_result == C){
                synch_type = 'C';
                synched = true;
            }else if(mat_mult_result == C_prime){
                synch_type = 'C';
                synched = true;
            }else if(mat_mult_result == D){
                synch_type = 'D';
                synched = true;
            }else{
                frame_index += 1;
            }
        }
		if (synched) {
			//std::cerr << synch_type << " at index " << frame_index << ", bitstream size = " << rds_bitstream.size() << std::endl;
			switch (synch_type) {
				case 'A':{
					unsigned short int new_PI = get_pi(frame);
					if(new_PI != saved_PI){
						fprintf(stderr, "Got PI: %x\n", new_PI);
						saved_PI = new_PI;
					}
					break;
				}
				case 'B':{
					unsigned char new_PT = get_pt(frame);
					if(new_PT != saved_PT){
						fprintf(stderr, "Got PT: %d ", get_pt(frame));
						print_pt(get_pt(frame));
						saved_PT = new_PT;
					}
					break;
				}
				default:
					break;
			}
		}

		// Look for next frame
		frame_index += 26;
		//synched = true;
		while(frame_index + 26 < rds_bitstream.size()){
			//std::cerr << frame_index << std::endl;
            frame = std::vector<int>(rds_bitstream.begin() + frame_index, rds_bitstream.begin() + frame_index + 26);
            //printRealVector(frame);
			mat_mult(mat_mult_result,frame,P);
			//printRealVector(mat_mult_result);rame_index + 26
            if(mat_mult_result == A){
                synch_type = 'A';
                synched = true;
            }else if(mat_mult_result == B){
                synch_type = 'B';
                synched = true;
            }else if(mat_mult_result == C){
                synch_type = 'C';
                synched = true;
            }else if(mat_mult_result == C_prime){
                synch_type = 'C';
                synched = true;
            }else if(mat_mult_result == D){
                synch_type = 'D';
                synched = true;
            }else{
                synched = false;
				//std::cerr << "Lost synchronization :( \n";
				frame_index += 1;
            }
			if(synched == true){
				//std::cerr << "Frame is type " << synch_type << " at index " << frame_index << ":)\n";
				switch (synch_type) {
					case 'A':{
						unsigned short int new_PI = get_pi(frame);
					if(new_PI != saved_PI){
						fprintf(stderr, "Got PI: %x\n", new_PI);
						saved_PI = new_PI;
					}
					break;
					}
				case 'B':{
					unsigned char new_PT = get_pt(frame);
					if(new_PT != saved_PT){
						fprintf(stderr, "Got PT: %d ", get_pt(frame));
						print_pt(get_pt(frame));
						saved_PT = new_PT;
					}
					break;
				}
				default:
					break;
				}
				frame_index += 26;
			}
		}

		// Save remaining bits to state
		int state_size = rds_bitstream.size() - frame_index;
		if(state_size > 0){
			frame_state.clear();
			frame_state.resize(state_size);
			for(int i = 0; i < state_size; i++){
				frame_state[i] = rds_bitstream[frame_index + i];
			}
		}else{
			frame_state = {};
		}
		//std::cerr << "Frame index at end of loop: " << frame_index << "\n";
		//std::cerr << "Size of frame state: " << frame_state.size() << "\n";
	}
}

int main(int argc, char* argv[])
{
//testing for runtime	
	int seed = std::time(0x0);
	std::cerr << "starting from seed" << std::hex <<seed<<std::dec<< "\n";

	// default mode 0, mono path
	char mode = '0';
	char path = 'm';
	int error = 0;

	if (argc >= 2) {
		if (argv[1][0] >= '0' && argv[1][0] <= '3') {
			mode = argv[1][0];
		} else {
			std::cerr << "Invalid mode.\n";
			error += 1;
		}
	}

	if (argc == 3) {
		if (argv[2][0] == 'm' || argv[2][0] == 's' || argv[2][0] == 'r') {
			path = argv[2][0];
		} else {
			std::cerr << "Invalid path.\n";
			error += 1;
		}
	}

	if (error != 0) {
		exit(-1);
	}

	float iq_fs = 2400000;
	float iq_fc = 100000;
	int num_taps = 101;
	float if_fs = 240000;
	float audio_fs = 48000;
	float audio_fc = 16000;
	unsigned int block_size = 100000;
	
	//declare new variables for mode selection
	int num_taps_final=num_taps;
	float if_fs_final=if_fs;
	//float audio_fs_final=audio_fs;
	float audio_fc_final=audio_fc;
	int U=1;
	int D=1;

	switch (mode) {
		case '0':{
			// mode 0 parameters, already configured as default
			D=5;
			num_taps_final=num_taps;
			iq_fs=2400000;
			if_fs=240000;
			if_fs_final=if_fs;
			audio_fs=48000;
			audio_fc=audio_fs/2;
			//std::cerr << iq_fs << " " << iq_fc << " " << iq_fs << " "
			break;
		}
		case '1':{
			// set mode 1 params here
			num_taps_final=num_taps;
			iq_fs=2304000;
			if_fs=288000;
			if_fs_final=if_fs;
			//audio_fs_final=audio_fs;
			audio_fc_final=audio_fc;
			D=8;
			break;
		}
		case '2':{
			// set mode 2 params here
			block_size = 100000;
			D=800;
			U=147;
			num_taps_final=num_taps*U;
			iq_fs=2400000;
			if_fs=240000;
			if_fs_final=if_fs*U;
			audio_fs = 44100;
			audio_fc=audio_fs/2;
			break;
		}
		case '3':{
			// set mode 3 params here
			//block_size = 60000;
			D=1600;
			U=441;
			num_taps_final=num_taps*U;
			iq_fs=1440000;
			if_fs=160000;
			if_fs_final=if_fs*U;
			audio_fs = 44100;
			audio_fc=audio_fs/2;
			block_size = 60000;
			break;
		}
		default:
			// catchall, do nothing
			break;
	}
	//std::cerr << "starting\n";
	std::vector<float> front_end_lpf;
	impulseResponseLPF(iq_fs, iq_fc, num_taps, front_end_lpf);
	std::vector<float> i_down_state, q_down_state;
	i_down_state.resize(front_end_lpf.size() - 1, 0.0);
	q_down_state.resize(front_end_lpf.size() - 1, 0.0);

	float i_ = 0.0;
	float q_ = 0.0;

	std::vector<float> block_data(block_size * 2);
	std::vector<float> i_data, q_data;
	std::vector<float> i_down, q_down;
	std::vector<float> fmout;

	// For audio thread
	std::queue<std::vector<float>> if_queue;
	std::mutex if_mutex;
	std::condition_variable if_cvar;
	
	// For RDS thread
	std::queue<std::vector<float>> rds_queue;
	std::mutex rds_mutex;
	std::condition_variable rds_cvar;
	
	// Shared between audio and RDS threads
	bool run = true;
	bool* rp = &run;

	std::thread audio_thread = std::thread(audio_function, std::ref(if_queue), \
		std::ref(if_mutex), std::ref(if_cvar), mode, path,
		block_size, D, U, num_taps,
		num_taps_final, if_fs, if_fs_final,
		audio_fs, audio_fc, audio_fc_final, rp);
	
	std::thread rds_thread;
	if (path == 'r'){
		rds_thread = std::thread(rds_function, std::ref(rds_queue), \
			std::ref(rds_mutex), std::ref(rds_cvar), mode, block_size, num_taps, 
			if_fs, rp);
	}

	//std::cerr << "start block calcs\n";
	for (unsigned int block_id = 0; std::cin.rdstate() == 0; block_id++) {
		// time measurement start
		auto start_time = std::chrono::high_resolution_clock::now();

		// read block from stdin
		
		readStdinBlockData(block_size * 2, block_id, block_data);
		//std::cerr << "starting conv and downsample\n";
		
		//std::cerr << "Read block " << block_id << std::endl;

		auto read_done = std::chrono::high_resolution_clock::now();

		blockconvdownsplit(i_down, q_down, block_data, front_end_lpf, q_down_state, int(iq_fs/if_fs));
		
		//std::cerr << "front end done\n" << i_down.size() << "size \n" ;
		auto frontend_filt_down_done = std::chrono::high_resolution_clock::now();

		// fm demod
		FM_Demod(i_down, q_down, fmout, &i_, &q_);
		//std::cerr << "dm demod done\n";
		auto demod_done = std::chrono::high_resolution_clock::now();

		// Audio thread
		std::unique_lock<std::mutex> lock_audio(if_mutex);
		while (if_queue.size() >= QUEUE_ELEMS) {
			if_cvar.wait(lock_audio);
		}

		if_queue.push(fmout);
		if_cvar.notify_one();
		lock_audio.unlock();
		
		// RDS thread
		if (path == 'r'){
			std::unique_lock<std::mutex> lock_rds(rds_mutex);
			while (rds_queue.size() >= QUEUE_ELEMS) {
				rds_cvar.wait(lock_rds);
			}

			// std::cerr << "FMOUT SIZE GOING INTO RDS QUEUE: " << fmout.size() << "\n";
			rds_queue.push(fmout);
			rds_cvar.notify_one();
			lock_rds.unlock();
		}

		// TIMING STUFF
		auto stop_time = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double, std::milli> Test_run_time = stop_time - start_time;
		//std::chrono::duration<double, std::milli> read_time = read_done - start_time;
		//std::chrono::duration<double, std::milli> deint_time = deinterleave_done - read_done;
		//std::chrono::duration<double, std::milli> front_time = frontend_filt_down_done - deinterleave_done;
		//std::chrono::duration<double, std::milli> demod_time = demod_done - frontend_filt_down_done;
		//std::chrono::duration<double, std::milli> back_time = backend_filt_down_done - demod_done;
		
		//std::cerr << "read time " << read_time.count() << " deint time " << deint_time.count() << " front time " << front_time.count() << " demod time " << demod_time.count() << " back time " << back_time.count() << std::endl;
		//std::cerr << "block processing took " << Test_run_time.count() << " ms" << std::endl;
	}

	// Stop running and join threads
	*rp = false;

	audio_thread.join();
	if (path == 'r'){
		rds_thread.join();
	}
	return 0;
}
