#include <stdio.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>
#include <string.h>
#include <time.h>
#include <math.h>


#include<cstdlib>
#include<ctime>

#include "main.h"

//Given an input file with M, N, P, Q, will generate a corresponding Sudoku board (not guaranteed to be solvable)
int generator(std::string input_file_name, std::string output_file_name){
	std::ifstream myfile;
	std::string line, item;
	std::vector<std::string> inputs;
	int M, N, P, Q;
	myfile.open(input_file_name.c_str());
	getline(myfile, line);
	int index = 0;
	while ((index = line.find(" ")) != std::string::npos){
		item = line.substr(0, index);
		inputs.push_back(item);
		line.erase(0, index + 1);
	}
	inputs.push_back(line);
	M = str_to_int(inputs[0]);
	N = str_to_int(inputs[1]);
	P = str_to_int(inputs[2]);
	Q = str_to_int(inputs[3]);


	if (N != P*Q){
		std::cout << "Invalid parameters. N != P*Q." << std::endl;
		std::exit(1);
	}
	if (M > N*N){
		std::cout << "Invalid parameters. M > N*N." << std::endl;
		std::exit(1);
	}

	myfile.close();

	std::vector<int> dim;
	dim.push_back(P);
	dim.push_back(Q);
	dim.push_back(N);
	dim.push_back(M);

	//return result;
	std::srand(std::time(0));
	std::pair<int, std::set<int> >** game_board;

	while (true){
		game_board = make_board(dim[2]);
		bool success = true;
		for (int i = 0; i<dim[3]; i++){		// for each M
			std::pair<int, int> cell = generate_cell(game_board, dim[2]);//generate a random, vacant cell
			if (!fill_cell(game_board, cell, dim[0], dim[1], dim[2])){	//fill cell with random token=
				success = false;										//if there are no more possible
				break;													//solutions for that cell, start over
			}
		}
		if (success)
			break;
	}

	std::ofstream output(output_file_name.c_str());
	output << dim[2] << " " << dim[0] << " " << dim[1] << std::endl; //print header info  << dim[3] << " " 
	for (int p = 0; p<dim[2]; p++){			//print resulting board
		for (int q = 0; q<dim[2]; q++){
			if (game_board[p][q].first > 9){
				output << num_to_letter(game_board[p][q].first) << " ";
			}
			else{
				output << game_board[p][q].first << " ";
			}
		}
		output << "\n";
	}
	return 0;
}

//Creates a blank board (initial constraints included)
std::pair<int, std::set<int> >** make_board(int N){
	std::set<int> setup;
	for (int i = 1; i<N + 1; i++){
		setup.insert(i);
	}
	std::pair<int, std::set<int> > constraint = std::make_pair(0, setup);

	std::pair<int, std::set<int> >** grid = new std::pair<int, std::set<int> >*[N];
	for (int i = 0; i<N; i++){
		grid[i] = new std::pair<int, std::set<int> >[N];
		for (int t = 0; t<N; t++){
			grid[i][t] = constraint;
		}
	}
	return grid;
};

//Updates the surrounding cells, given a  particular cell and the value that it was just updated with
void update_cells(std::pair<int, std::set<int> > **grid, int row, int col, int N, int val, int p, int q){
	for (int i = 0; i<N; i++){
		grid[row][i].second.erase(val);
		grid[i][col].second.erase(val);
	}

	int row_start = row - row%p;
	int col_start = col - col%q;

	for (int i = row_start; i<(p + row_start); i++){
		for (int t = col_start; t<(q + col_start); t++){
			if (i == row && t == col){}
			//do nothing -- for now, to avoid uneccesary iterating 
			else
				grid[i][t].second.erase(val);
		}
	}

};

//Checks whether or not a particular value will violate a particular cell's constraints
bool constraint_violation(std::pair<int, std::set<int> > **grid, int row, int col, int N, int val, int p, int q){

	for (int i = 0; i<N; i++){
		if (grid[row][i].first == val){ //check associated row
			//std::cout << "{row} conflicts with row: " << row << ", col: " << i << std::endl;
			return false;
		}
		if (grid[i][col].first == val){ //check associated column
			//std::cout << "{col} conflicts with row: " << i << ", col: " << col << std::endl;
			return false;
		}
	}

	int row_start = row - row%p;
	int col_start = col - col%q;

	for (int i = row_start; i<(p + row_start); i++){ //check associated block
		for (int t = col_start; t<(q + col_start); t++){
			if (grid[i][t].first == val){
				//std::cout << "{block} conflicts with row: " << i << ", col: " << t << std::endl;
				return false;
			}
		}
	}
	return true;
};

//Generates a random, vacant cell, given a board
std::pair<int, int> generate_cell(std::pair<int, std::set<int> >** grid, int N){
	int row, col;
	while (true){
		row = std::rand() % (N)+0;
		col = std::rand() % (N)+0;
		if (grid[row][col].first == 0)
			break;
	}
	return std::make_pair(row, col);
};


//Fills a cell with a random number (within its constraints), provided a random cell to begin with
//Returns a bool depending on whether the token was able to be filled or not (i.e. whether it's list
//of possible tokens is empty or not)
bool fill_cell(std::pair<int, std::set<int> > **grid, std::pair<int, int> cell, int p, int q, int N){
	int row = cell.first;
	int col = cell.second;

	//while(true){ 
	if (grid[row][col].second.size() == 0)
		return false; //if cell never satisfies constraints, fail globally
	int index = rand() % grid[row][col].second.size();
	std::set<int>::const_iterator num(grid[row][col].second.begin());
	advance(num, index);

	//if(constraint_violation(grid, row, col, N, *num, p, q)){
	grid[row][col].first = (*num);
	update_cells(grid, row, col, N, *num, p, q);
	//	break; 
	//}
	//}
	return true;
};

//Convert type int to a chacter representing a letter
char num_to_letter(int num){
	return static_cast<char>(num + 55);
};

//Converts a character representing a letter to a type int
int letter_to_num(char letter){
	if (letter - '0' < 10)
		return letter - '0';
	else
		return static_cast<int>(letter - 55);
};

//Converts type string to type int
int str_to_int(std::string inputs){
	std::stringstream ss;
	ss << inputs;
	int to_return;
	ss >> to_return;
	return to_return;
};

//Given an initial sudoku board via a text file, will read and store the represented board into a 2D array, as 
//well as update it's corresponding domain values according to the values given
std::pair<int, std::set<int> >** problem_reader(std::string input_file_name){
	std::ifstream myfile;
	std::string line, item;
	std::vector<std::string> inputs;
	int N, P, Q;
	myfile.open(input_file_name.c_str());
	getline(myfile, line);
	int index = 0;
	while ((index = line.find(" ")) != std::string::npos){
		item = line.substr(0, index);
		inputs.push_back(item);
		line.erase(0, index + 1);
	}
	inputs.push_back(line);


	N = str_to_int(inputs[0]);
	P = str_to_int(inputs[1]);
	Q = str_to_int(inputs[2]);
	if (N != P*Q){
		std::cout << "Invalid parameters. N != P*Q." << std::endl;
		std::exit(1);
	}

	std::pair<int, std::set<int> >** game_board;
	//getline(myfile, line);
	int row = 0;
	game_board = make_board(N);

	while (!myfile.eof()){
		getline(myfile, line);
		if (line == "" || line == "\n")
			break;
		line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
		for (int i = 0; i < N; i++){
			game_board[row][i].first = letter_to_num(line[i]);
			update_cells(game_board, row, i, N, letter_to_num(line[i]), P, Q);
		}
		row += 1;
	}
	myfile.close();

	return game_board;
};


class BTSSolver{

	//Variables
	int Assignments;
	int Backtracks;
	time_t StartTime;
	time_t EndTime;
	int TimeLimit;
	std::pair<int, std::set<int> >**  board;
	int N;
	int p;
	int q;
	std::pair<int, std::set<int> >** failure;
	std::string status;
	bool FC;
	bool LCV_constraint;
	bool MRV;
	bool DH;

	//Constructors
public: BTSSolver(std::pair<int, std::set<int> >**  board, int N, int p, int q, int TimeLimit){
			this->board = board;
			Assignments = 0;
			Backtracks = 0;
			this->N = N;
			this->p = p;
			this->q = q;
			std::set<int> empt;
			StartTime = time(0);
			this->TimeLimit = TimeLimit;
			status = "error";
			FC = false;
			LCV_constraint = false;
			MRV = false;
			DH = false;

}

		int getBacktracks(){
			return Backtracks;
		}

		int getAssignments(){
			return Assignments;
		}

		time_t getEndTime(){
			return EndTime;
		}

		time_t getStartTime(){
			return StartTime;
		}

		std::string getStatus(){
			return status;
		}

		void setFC(){
			FC = true;
		}

		void setLCV_constraint(){
			LCV_constraint = true;
		}

		void setMRV(){
			MRV = true;
		}

		void setDH(){
			DH = true;
		}

		//Loops over the board, and finds the first unassigned variable
		std::pair<int, int> getFirstUnassignedVar(){
			for (int i = 0; i < N; i++){
				for (int j = 0; j < N; j++){
					if (board[i][j].first == 0){
						//std::cout << i << " " << j << std::endl;
						return std::make_pair(i, j);
					}
				}
			}
		}

		//Loops over the board, and finds the cell with highest DH value
		std::pair<int, int> getDH(){
			int track = 0;
			int max = 0;
			std::pair<int, int> DH_cell;

			for (int i = 0; i < N; i++){
				for (int j = 0; j < N; j++){
					if (track == 0 && board[i][j].first == 0){
						DH_cell = std::make_pair(i, j);
						max = get_DH_constraint(board, i, j, N, p, q);
						track = 1;
					}
					else{
						if (board[i][j].first == 0){
							int temp = get_DH_constraint(board, i, j, N, p, q);
							if (temp > max){
								max = temp;
								DH_cell = std::make_pair(i, j);
							}
						}
					}
				}
			}
			return DH_cell;
		}

		//Loops over a particular domain, and finds the minimum remaning value. 
        //This function may also have the option of implementing DH as a tie-breaker, if DH is activated as a token.
		std::pair<int, int> getMRV(){
			std::pair <int, int> MRV;
			int lowestRemaining = N;
			int track = 0;
			for (int i = 0; i < N; i++){
				for (int j = 0; j < N; j++){
					if (track == 0 && board[i][j].first == 0){
						MRV = std::make_pair(i, j);
						track = 1;
					}
					else{
						if (board[i][j].first == 0){
							if (DH == true){
								if (board[i][j].second.size() == lowestRemaining){
									MRV = getDH_MRV(std::make_pair(i, j), MRV);
								}
								else if (board[i][j].second.size() < lowestRemaining){
									MRV = std::make_pair(i, j);
									lowestRemaining = board[i][j].second.size();
								}
							}
							else{
								if (board[i][j].second.size() <= lowestRemaining){
									MRV = std::make_pair(i, j);
									lowestRemaining = board[i][j].second.size();
								}
							}
						}
					}
				}
			}
			return MRV;
		}

		//Acts as a tie-breaker for the MRV heuristic, choosing instead the cell that is associated with the highest number of unassigned cells around it
		std::pair<int, int> getDH_MRV(std::pair<int, int> p1, std::pair<int, int> p2){
			int a = get_DH_constraint(board, p1.first, p1.second, N, p, q);
			int b = get_DH_constraint(board, p2.first, p2.second, N, p, q);
			if (a > b)
				return p1;
			else
				return p2;
		}

		//Returns the number of unassigned cells surrounding a particular cell (row-wise, col-wise, and sub-grid wise)
		int get_DH_constraint(std::pair<int, std::set<int> > **grid, int row, int col, int N, int p, int q){
			std::set<std::pair<int, int> > counted;
			int count = 0;
			for (int i = 0; i<N; i++){
				if (i != col && grid[row][i].first == 0){
					counted.insert(std::make_pair(row, i));
					count += 1;
				}

				if (i != row && grid[i][col].first == 0){
					counted.insert(std::make_pair(i, col));
					count += 1;
				}
			}
			int row_start = row - row%p;
			int col_start = col - col%q;

			for (int i = row_start; i<(p + row_start); i++){
				for (int t = col_start; t< (q + col_start); t++){
					if (i == row && t == col){} //do nothing for now
					else{
						bool is_in = counted.find(std::make_pair(i, t)) != counted.end();
						if (grid[i][t].first == 0 && !is_in)
							count += 1;
					}
				}
			}
			return count;
		}


		//Finds LCV, given a particular cell
		int find_LCV(std::pair<int, int> var, std::set<int> domain){
			int LCV;
			int num_of_constraints;
			int temp;
			int track = 0;

			std::set<int>::iterator it;
			for (it = domain.begin(); it != domain.end(); it++){
				if (track == 0){
					num_of_constraints = update_cells_preview(board, var.first, var.second, N, (*it), p, q);
					LCV = (*it);
					track = 1;
				}
				temp = update_cells_preview(board, var.first, var.second, N, (*it), p, q);
				if (temp < num_of_constraints){
					num_of_constraints = temp;
					LCV = (*it);
				}
			}
			return LCV;
		}

		bool isComplete(std::pair<int, std::set<int> >** assignment){
			for (int i = 0; i < N; i++){
				for (int j = 0; j < N; j++){
					if (assignment[i][j].first == 0){
						return false;
					}
				}
			}
			return true;
		}

		std::set<int> getDomain(std::pair<int, int> var){
			return board[var.first][var.second].second;
		}



		//Returns the number of cells that would have their value crossed out
		int update_cells_preview(std::pair<int, std::set<int> > **grid, int row, int col, int N, int val, int p, int q){
			std::set<std::pair<int, int> > counted;
			int count = 0;
			for (int i = 0; i<N; i++){
				bool is_in = grid[row][i].second.find(val) != grid[row][i].second.end();
				if (i != col && grid[row][i].first == 0 && is_in == true){
					counted.insert(std::make_pair(row, i));
					count += 1;
				}

				is_in = grid[i][col].second.find(val) != grid[i][col].second.end();
				if (i != row && grid[i][col].first == 0 && is_in == true){
					counted.insert(std::make_pair(i, col));
					count += 1;
				}
			}
			int row_start = row - row%p;
			int col_start = col - col%q;

			for (int i = row_start; i<(p + row_start); i++){
				for (int t = col_start; t< (q + col_start); t++){
					if (i == row && t == col){} //do nothing -- for now, to avoid uneccesary iterating 
					else{
						bool in_counted = (counted.find(std::make_pair(i, t)) != counted.end());
						bool is_in_subgrid = (grid[i][t].second.find(val) != grid[i][t].second.end());
						if (grid[i][t].first == 0 && is_in_subgrid == true && !in_counted)
							count += 1;
					}
				}
			}
			return count;
		}

		//Keeps track of the cells and their value that was crossed out, in order to correctly backtrack
		std::set<std::pair<std::pair<int, int>, int> > FCUpdateCells(std::pair<int, std::set<int> > **grid, int row, int col, int N, int val, int p, int q){
			std::set<std::pair<std::pair<int, int>, int> > deleted;

			for (int i = 0; i<N; i++){
				deleted.insert(std::make_pair(std::make_pair(row, i), val));
				deleted.insert(std::make_pair(std::make_pair(i, col), val));
				grid[row][i].second.erase(val);
				grid[i][col].second.erase(val);
			}

			int row_start = row - row%p;
			int col_start = col - col%q;

			for (int i = row_start; i<(p + row_start); i++){
				for (int t = col_start; t<(q + col_start); t++){
					if (i == row && t == col){}
					//do nothing -- for now, to avoid uneccesary iterating 
					else{
						deleted.insert(std::make_pair(std::make_pair(i, t), val));
						grid[i][t].second.erase(val);
					}
				}
			}

			return deleted;
		}

		//Recursively backtracks through the board, iteratively filling in each cell until a solution or "no solution" is found.
		//This function may also implement the usage of other tokens, such as FC, MRV, and LCV, depending on what was inputted by the user
		std::pair<int, std::set<int> >** RecursiveBackTrack(std::pair<int, std::set<int> >** assignment){
			std::set<std::pair<std::pair<int, int>, int> > deleted;
			std::pair<int, int> var;
			if (isComplete(assignment)){
				return assignment;
			}
			if (MRV == true)
				var = getMRV();
			else{
				if (DH == true){
					var = getDH();
				}
				else
					var = getFirstUnassignedVar();
			}
			std::set<int> domain = getDomain(var);
			for (std::set<int>::iterator it = domain.begin(); it != domain.end(); it++){
				if (time(0) - StartTime >= TimeLimit){
					status = "timeout";
					break;
				}
				if (constraint_violation(assignment, var.first, var.second, N, *it, p, q)){
					assignment[var.first][var.second].first = (*it);
					Assignments++;
					if (FC == true)
						deleted = FCUpdateCells(assignment, var.first, var.second, N, *it, p, q);
					std::pair<int, std::set<int> >** result = RecursiveBackTrack(assignment);
					if (result != failure){
						status = "success";
						return result;
					}
					assignment[var.first][var.second].first = 0;
					for (std::set<std::pair<std::pair<int, int>, int> >::iterator it2 = deleted.begin(); it2 != deleted.end(); it2++){
						assignment[it2->first.first][it2->first.second].second.insert(it2->second);
					}
					Backtracks++;
				}
			}
			return failure;
		}

		std::pair<int, std::set<int> >** LCVRecursiveBackTrack(std::pair<int, std::set<int> >** assignment){ //using LCV
			std::set<std::pair<std::pair<int, int>, int> > deleted;
			std::pair<int, int> var;
			if (isComplete(assignment)){
				return assignment;
			}
			if (MRV == true)
				var = getMRV();
			else{
				if (DH == true)
					var = getDH();
				else
					var = getFirstUnassignedVar();
			}
			std::set<int> domain = getDomain(var);
			//for (std::set<int>::iterator it = domain.begin(); it != domain.end(); it++){
			while (status != "success"){							// ---Added this while loop to pick new LCV's until 1 works	
				if (domain.size() == 0){							// ---Added this If to check when no more possible LCV's to try
					return failure;
				}
				int LCV = find_LCV(var, domain);
				domain.erase(LCV);										// ---Added this line to erase an LCV we already checked from local variable domain
				if (time(0) - StartTime >= TimeLimit){
					status = "timeout";
					return failure;
				}
				if (constraint_violation(assignment, var.first, var.second, N, LCV, p, q)){
					assignment[var.first][var.second].first = LCV;
					Assignments++;
					if (FC == true)
						deleted = FCUpdateCells(assignment, var.first, var.second, N, LCV, p, q);
					std::pair<int, std::set<int> >** result = LCVRecursiveBackTrack(assignment);
					if (result != failure){
						status = "success";
						return result;
					}
					assignment[var.first][var.second].first = 0;
					for (std::set<std::pair<std::pair<int, int>, int> >::iterator it2 = deleted.begin(); it2 != deleted.end(); it2++){
						assignment[it2->first.first][it2->first.second].second.insert(it2->second);
					}
					Backtracks++;
				}
			}
			return failure;
		}

		std::pair<int, std::set<int> >** BackTrackSearch(){
			if (LCV_constraint == true)
				return LCVRecursiveBackTrack(board);
			return RecursiveBackTrack(board);
		}


};

//Returns a set of substrings from s, delimited by the chosen "delimiter" input
std::vector<std::string> split_string(std::string str, char delim){
	std::vector<std::string> tokens;
	std::stringstream ss(str);
	std::string token;
	while (std::getline(ss, token, delim)){
		tokens.push_back(token);
	}
	return tokens;
}

//Given an input file containing M N P Q, generate boards and write them to file names listed in output.
void GenerateFiles(std::string input, std::string output){
	std::ifstream myfile;
	std::string line;
	std::string fileoutline;
	std::ifstream outfile;

	myfile.open(input.c_str());
	outfile.open(output.c_str());
	getline(myfile, line);
	getline(outfile, fileoutline);
	while (fileoutline != ""){
		generator(line, fileoutline);
		getline(outfile, fileoutline);
	}
}

//Used to test the average results (total nodes, time taken, and standard deviation) of 10 randomly generated sudoku boards
void RunTests(std::string input, std::string output, int time_limit, std::string token1, std::string token2, std::string token3){
	double nodes_totals = 0.0;
	double time_totals = 0.0;
	int completed = 0;
	std::vector<double> avgs;

	std::string boards = "board.txt";

	generator(input, boards);

	for (int k = 0; k < 10; k++){
		std::cout << "Started New Test" << std::endl;
		time_t total_start = time(0);

		std::ifstream myfile;
		std::string line, item;
		std::vector<std::string> inputs;
		int N, P, Q;
		myfile.open(boards.c_str());
		getline(myfile, line);
		int index = 0;
		while ((index = line.find(" ")) != std::string::npos){
			item = line.substr(0, index);
			inputs.push_back(item);
			line.erase(0, index + 1);
		}
		inputs.push_back(line);
		N = str_to_int(inputs[0]);
		P = str_to_int(inputs[1]);
		Q = str_to_int(inputs[2]);


		std::pair<int, std::set<int> >**  board = problem_reader(boards);

		BTSSolver BTS = BTSSolver::BTSSolver(board, N, P, Q, time_limit);

		if (token1 == "FC" || token2 == "FC" || token3 == "FC")
			BTS.setFC();
		if (token1 == "LCV" || token2 == "LCV" || token3 == "LCV")
			BTS.setLCV_constraint();
		if (token1 == "MRV" || token2 == "MRV" || token3 == "MRV")
			BTS.setMRV();
		if (token1 == "DH" || token2 == "DH" || token3 == "DH")
			BTS.setDH();



		time_t preproc_start = time(0);
		time_t preproc_end = time(0);

		std::pair<int, std::set<int> >** solved = BTS.BackTrackSearch();

		time_t backtrack_end = time(0);

		int assignments = BTS.getAssignments();
		int backtracks = BTS.getBacktracks();

		std::ofstream file_out;
		file_out.open(output.c_str());
		file_out << "TOTAL_START=" << total_start << std::endl;
		file_out << "PREPROCESSING_START=" << preproc_start << std::endl;
		file_out << "PREPROCESSING_DONE=" << preproc_end << std::endl;
		file_out << "SEARCH_START=" << BTS.getStartTime() << std::endl;
		file_out << "SEARCH_DONE=" << backtrack_end << std::endl;
		file_out << "SOLUTION_TIME=" << backtrack_end - BTS.getStartTime() << std::endl;
		file_out << "STATUS=" << BTS.getStatus() << std::endl;



		if (BTS.getStatus() != "success"){
			file_out << "SOLUTION=(";
			for (int i = 0; i < N*N; i++){
				if (i == (N*N) - 1){
					file_out << 0 << ")";
				}
				else{
					file_out << 0 << ",";
				}
			}
		}
		else{
			file_out << "SOLUTION=(";
			for (int i = 0; i < N; i++){
				for (int j = 0; j < N; j++){
					if (solved[i][j].first > 9){
						if (i == N - 1 && j == N - 1){
							file_out << num_to_letter(solved[i][j].first) << ")";
						}
						else{
							file_out << num_to_letter(solved[i][j].first) << ",";
						}
					}
					else{
						if (i == N - 1 && j == N - 1){
							file_out << solved[i][j].first << ")";
						}
						else{
							file_out << solved[i][j].first << ",";
						}
					}
				}
			}
		}

		file_out << std::endl;
		file_out << "COUNT_NODES=" << assignments << std::endl;
		file_out << "COUNT_DEADENDS=" << backtracks << std::endl;

		if (BTS.getStatus() == "success" || BTS.getStatus() == "error"){
			 completed += 1;
		}
		//else{
			//generator(input, boards);
			//continue;
		//}

		time_totals += difftime(backtrack_end, BTS.getStartTime());
		nodes_totals += assignments;
		avgs.push_back(difftime(backtrack_end, BTS.getStartTime()));

		generator(input, boards);
	}
		double avg_time;
		double avg_nodes;
		double solvable;
		double std_dev;


		avg_time = time_totals / 10.0;
		avg_nodes = nodes_totals / 10.0;
		solvable = completed / 10.0;
		solvable = solvable * 100.0;

		double sum=0;
		double meansub;
		for (std::vector<double>::iterator it = avgs.begin(); it != avgs.end(); ++it){
			meansub = (*it - avg_time);
			meansub = meansub*meansub;
			sum += meansub;
		}
		std_dev = sum / 10.0;
		std_dev = sqrt(std_dev);

		std::cout << "Avg time= " << avg_time << std::endl;
		std::cout << "Avg nodes= " << avg_nodes << std::endl;
		std::cout << "% Solved= " << solvable << std::endl;
		std::cout << "std_dev= " << std_dev << std::endl;
}

int main(int argc, char* argv[] ){
	//--Below, commented out code was used for compiling/testing on Visual Studios--
	//std::string input, output, tokens;
	//int time_limit;

	//std::cout << "Enter input file name: ";
	//std::cin >> input;

	//std::cout << "Enter output file name: ";
	//std::cin >> output;

	//std::cout << "Enter a time limit in seconds: ";
	//std::cin >> time_limit;

	//std::cout << "Enter tokens separated by spaces: ";
	//std::cin.ignore(256, '\n');
	//std::getline(std::cin, tokens);

	//std::vector<std::string> token_list = split_string(tokens, ' ');

	std::string input = argv[1];
	std::string output = argv[2];
	int time_limit = std::atoi(argv[3]);
	

	if (argc > 4){
		for (int i = 4; i < argc; i++){
			if (std::string(argv[i]) == "GEN"){
				generator(input, input);
			}
			else{
				continue;
			}
		}
	}


	//if ((std::find(token_list.begin(), token_list.end(), "GEN") != token_list.end()))
	//	generator(input,input);

	time_t total_start = time(0);

	std::ifstream myfile;
	std::string line, item;
	std::vector<std::string> inputs;
	int N, P, Q;
	myfile.open(input.c_str());
	getline(myfile, line);
	int index = 0;
	while ((index = line.find(" ")) != std::string::npos){
		item = line.substr(0, index);
		inputs.push_back(item);
		line.erase(0, index + 1);
	}
	inputs.push_back(line);
	N = str_to_int(inputs[0]);
	P = str_to_int(inputs[1]);
	Q = str_to_int(inputs[2]);


	std::pair<int, std::set<int> >**  board = problem_reader(input);

	BTSSolver BTS = BTSSolver::BTSSolver(board, N, P, Q, time_limit);

	if (argc > 4){
		for (int i = 4; i < argc; i++){
			if (std::string(argv[i]) == "FC")
				BTS.setFC();
			if (std::string(argv[i]) == "LCV")
				BTS.setLCV_constraint();
			if (std::string(argv[i]) == "MRV")
				BTS.setMRV();
			if (std::string(argv[i]) == "DH")
				BTS.setDH();
		}
	}

	//if ((std::find(token_list.begin(), token_list.end(), "FC") != token_list.end()))
	//	BTS.setFC();
	//if ((std::find(token_list.begin(), token_list.end(), "LCV") != token_list.end()))
	//	BTS.setLCV_constraint();
	//if ((std::find(token_list.begin(), token_list.end(), "MRV") != token_list.end()))
	//	BTS.setMRV();
	//if ((std::find(token_list.begin(), token_list.end(), "DH") != token_list.end()))
	//	BTS.setDH();


	time_t preproc_start = time(0);
	time_t preproc_end = time(0);

	std::pair<int, std::set<int> >** solved = BTS.BackTrackSearch();

	time_t backtrack_end = time(0);

	int assignments = BTS.getAssignments();
	int backtracks = BTS.getBacktracks();

	std::ofstream file_out;
	file_out.open(output.c_str());
	file_out << "TOTAL_START=" << total_start << std::endl;
	file_out << "PREPROCESSING_START=" << preproc_start << std::endl;
	file_out << "PREPROCESSING_DONE=" << preproc_end << std::endl;
	file_out << "SEARCH_START=" << BTS.getStartTime() << std::endl;
	file_out << "SEARCH_DONE=" << backtrack_end << std::endl;
	file_out << "SOLUTION_TIME=" << difftime(backtrack_end, BTS.getStartTime())<< std::endl;
	file_out << "STATUS=" << BTS.getStatus() << std::endl;

	if (BTS.getStatus() != "success"){
		file_out << "SOLUTION=(";
		for (int i = 0; i < N*N; i++){
			if (i == (N*N) - 1){
				file_out << 0 << ")";
			}
			else{
				file_out << 0 << ",";
			}
		}
	}
	else{
		file_out << "SOLUTION=(";
		for (int i = 0; i < N; i++){
			for (int j = 0; j < N; j++){
				if (solved[i][j].first > 9){
					if (i == N - 1 && j == N - 1){
						file_out << num_to_letter(solved[i][j].first) << ")";
					}
					else{
						file_out << num_to_letter(solved[i][j].first) << ",";
					}
				}
				else{
					if (i == N - 1 && j == N - 1){
						file_out << solved[i][j].first << ")";
					}
					else{
						file_out << solved[i][j].first << ",";
					}
				}
			}
		}
	}

	file_out << std::endl;
	file_out << "COUNT_NODES=" << assignments << std::endl;
	file_out << "COUNT_DEADENDS=" << backtracks << std::endl;


	////for debug, prints out board in grid form
	//file_out << std::endl;
	//for (int i = 0; i < N; i++){
	//	for (int j = 0; j < N; j++){
	//		file_out << solved[i][j].first << " ";
	//	}
	//	file_out << std::endl;
	//}
	//std::cout << "Val of N was: " << N << std::endl;
	return 0;
	return 0;
}