#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char** argv) {
    std::vector<std::vector<float> > grid;
    std::ifstream file;
    std::string line;
    file.open(argv[1], std::ios::in);
    int row = 0, col = 0;
    while(std::getline(file, line)) {
        std::vector<std::string> tokens;
        std::string token;
        std::istringstream iss(line);
        while (std::getline(iss, token, ',')) {
            float val = strtof(token.c_str(), NULL);
            grid[row][col] = val;
            col++;
        }
        col = 0;
        row++;
    }

    for (int i = 0; i < grid.size(); i++) {
        for (int j = 0; j < grid[i].size(); j++) {
            std::cout << grid[i][j] << " ";
        }
        std::cout << "\n";
    }

    return 0;
}
