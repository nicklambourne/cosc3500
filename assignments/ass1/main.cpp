#include<iostream>
#include<fstream>
#include<string>
using namespace std;


// Finds the longest common substring between the two provided strings
string lcs(string x, string y) {
    int len_x = x.length();
    int len_y = y.length();
    vector<vector<int>> table (
        len_x,
        vector<int> (len_y, 0)        
    );
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= len_y; j++) {
            if (i == 0 || j == 0) {
                table[i][j] = 0;
            } else if (x[i - 1] == y[j - 1]) {
                table[i][j] = table[i - 1][j - 1] + 1;
            } else {
                table[i][j] = max(table[i - 1][j], table[i, j-1]);
            }
        }
    }
    return "value";
}


int main(int argc, char** argv) {
    // Check args and capture input/output file names
    if (argc != 3) {
        cerr << "Invalid call:\n Usage: ass1 <input_file> <output_file>" \
             << endl;
        exit(1);
    }
    string file_name = argv[1];
    string output_file_name = argv[2];

    ifstream input_file;
    input_file.open(file_name);

    string string_x;
    string string_y;

    if (input_file.is_open()) {
        string buffer;
        if (!getline(input_file, string_x) ||
            !getline(input_file, string_y)) {
            cerr << "Error reading from provided file." << endl;
            exit(3);
        }
        // cout << string_x << endl;
        // cout << string_y << endl;
    } else {
        cerr << "Could not open provided file." << endl;
        exit(2);
    }

    string longest = lcs(string_x, string_y);
    cout << longest << endl;
    
    ofstream output_file;
    output_file.open(output_file_name);

    if (output_file.is_open()) {
        output_file << longest;
    } else {
        cerr << "Could not write to given output file." << endl;
        exit(4);
    }
}
