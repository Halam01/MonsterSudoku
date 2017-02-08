#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>
#include <string.h>
#include <time.h>

#include<cstdlib>
#include<ctime>

int generator(std::string input_file_name, std::string output_file_name);

//Creates a blank board (initial constraints included)
std::pair<int, std::set<int> >** make_board(int N);

//Updates the surrounding cells, given a  particular cell and the value that it was just updated with
void update_cells(std::pair<int, std::set<int> > **grid, int row, int col, int N, int val, int p, int q);

//Checks whether or not a particular value will violate a particular cell's constraints
bool constraint_violation(std::pair<int, std::set<int> > **grid, int row, int col, int N, int val, int p, int q);

//Generates a random, vacant cell, given a board
std::pair<int, int> generate_cell(std::pair<int, std::set<int> >** grid, int N);

bool fill_cell(std::pair<int, std::set<int> > **grid, std::pair<int, int> cell, int p, int q, int N);

char num_to_letter(int num);

int str_to_int(std::string inputs);

int letter_to_num(char letter);

std::pair<int, std::set<int> >** problem_reader(std::string input_file_name);

int getBacktracks();

int getAssignments();

time_t getEndTime();

time_t getStartTime();

std::string getStatus();

std::pair<int, int> getFirstUnassignedVar();

bool isComplete(std::pair<int, std::set<int> >** assignment);

std::pair<int, std::set<int> >** RecursiveBackTrack(std::pair<int, std::set<int> >** assignment);

std::pair<int, std::set<int> >** BackTrackSearch();

//new additions as of 3/2/16
std::pair<int, std::set<int> >** LCVRecursiveBackTrack(std::pair<int, std::set<int> >** assignment);
int update_cells_preview(std::pair<int, std::set<int> > **grid, int row, int col, int N, int val, int p, int q);
int find_LCV(std::pair<int, int> var, std::set<int> domain);
std::pair<int, int> getDH_MRV(std::pair<int, int>, std::pair<int, int>);
int get_DH_constraint(std::pair<int, std::set<int> > **grid, int row, int col, int N, int p, int q);
std::pair<int, int> getDH();
std::pair<int, int> getMRV();
std::vector<std::string> split_string(std::string str, char delim);

#endif MAIN_H