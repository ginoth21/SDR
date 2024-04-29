/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_RDS_H
#define DY4_RDS_H

// add headers as needed
#include <iostream>
#include <vector>

// declaration of a function prototypes
void rrc(float, unsigned short int, std::vector<float> &);
void frame_sync(const std::vector<int> &);
void mat_mult(std::vector<int> &, std::vector<int> &, std::vector<std::vector<int>> &);
// void differential_decode(std::vector<int> &, const std::vector<int> &, int);
unsigned short int get_pi(const std::vector<int> &);
unsigned char get_pt(const std::vector<int> &);
void print_pt(unsigned char);

#endif // DY4_RDS_H
