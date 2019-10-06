#include<algorithm>
#include<fstream>
#include<iostream>
#include<string>
#include<vector>

using namespace std;


#define CELL_LENGTH 3


// Left pads a string with spaces to reach CELL_LENGTH
string normalise_length(string number) {
    string normalised_string = number;
    while (normalised_string.length() < CELL_LENGTH) {
        normalised_string = " " + normalised_string;
    }
    return normalised_string;
}


// Produces a printable string representation of the given table
string show_table(vector<vector<int>> table, string x, string y) {
    x = " " + x;
    y = " " + y;
    string str_table = "      ";
    for (int i = 0; i < (int) table.size(); i++) {
        str_table += "  ";
        str_table += y[i];
    }
    str_table.append("\n      ");
    for (int i = 0; i < (int) table[1].size(); i++) {
        str_table.append(normalise_length(to_string(i)));
    }
    str_table += "\n";
    for (int i = 0; i < (int) table.size(); i++) {
        str_table += "  ";
        str_table +=  x[i];
        str_table.append(normalise_length(to_string(i)));
        for (int j = 0; j < (int) table[i].size(); j++) {
            string value = to_string(table[i][j]);
            str_table.append(normalise_length(value));
        }
        str_table.append("\n");
    }
    return str_table;
}


// Reverses a string in place
void reverse_string(string& x) {
    reverse(x.begin(), x.end());
}


// Takes two strings and an LCS table and reconstructs the actual LCS string
string reconstruct_lcs(vector<vector<int>> table, string x, string y) {
    string lcs = "";
    x = " " + x;
    y = " " + y;
    int i = x.length();
    int j = y.length();
    while (i != 0 && j != 0) {
        if (x[i] == y[j]) {  // Characters match
            lcs.append(x, i, 1);
            i = max(i - 1, 0);
            j = max(j - 1, 0);
        } else if (i == 0) {  // Cannot go any further left
            j--;
        } else if (j == 0) {  // Cannot go any further up
            i--;
        } else if (table[i][j-1] > table[i-1][j]) {  // Left is larger  
            j--;
        } else { // Up is larger or they are the same
            i--;
        }
    }
    reverse_string(lcs);
    return lcs;
}


// Finds the longest common subsequence between the two provided strings
string lcs(string x, string y) {
    int m = x.length();
    int n = y.length();

    if (m == 0 || n == 0) {
        return "";
    }

    vector<vector<int>> table (
        m,
        vector<int> (n, 0)        
    );

    int i, j, diag;
    //#pragma omp parallel default(none) shared(x, y, table, m, n) private(i, j, diag) 
    //#pragma omp for
    for (diag = 2; diag < m + n; diag++) {
        //#pragma omp parallel for
        cout << "?" << endl; 
        for (i = min(m, diag - 1); i < max(1, diag - 1); i++) {
            j = diag - 1;
            
            if (j > n - 1) {
                continue;
            }
            cout << show_table(table, x, y) << endl;
            if (x[i - 1] == y[j - 1]) {
                table[i][j] = table[i - 1][j - 1] + 1;
            } else {
                table[i][j] = max(table[i - 1][j], table[i][j - 1]);
            }
        }
    }

    // cout << show_table(table, x, y);
    return reconstruct_lcs(table, x, y);
}


int main(int argc, char** argv) {
    // Check args and capture input/output file names
    if (argc != 3) {
        cerr << "Invalid call:\n Usage: ass1 <input_file> <output_file>" \
             << endl;
        exit(1);
    }

    // Get arguments
    string file_name = argv[1];
    string output_file_name = argv[2];

    ifstream input_file;
    input_file.open(file_name);

    string string_x;
    string string_y;

    // Read in sequences to compare
    if (input_file.is_open()) {
        string buffer;
        if (!getline(input_file, string_x) ||
            !getline(input_file, string_y)) {
            cerr << "Error reading from provided file." << endl;
            exit(3);
        }
    } else {
        cerr << "Could not open provided file." << endl;
        exit(2);
    }

    // Get LCS and write to stdout and file
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
