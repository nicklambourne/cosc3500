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
    for (int i = 0; i < (int) table.size(); i++) {
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
    int i = x.length();
    int j = y.length();
    while (!(i == 0 && j == 0)) {
        printf("(i: %d, j: %d): %d\n", i, j, table[i][j]);
        if (x[i] == y[j]) {  // Characters match
            lcs.append(x, i, 1);
            printf("%c\n", x[i]);
        }  
        if (i == 0 && j == 0) {  // No more to parse
            break;
        } else if (i == 0) {  // Cannot go any further left
            j--;
        } else if (j == 0) {  // Cannot go any futher up
            i--;
        } else if (table[i-1][j] > table[i][j-1]) {  // Up is larger  
            i--;
        } else if (table[i-1][j] < table[i][j-1]) {  // Left is larger
            j--;
        } else if (table[i-1][j] == table[i][j-1]) {  // Equal, move diagonal
            i--;
            j--;            
        } else {  // Decrement the larger index
            if (i > j) {
                i--;
            } else {
                j--;
            }
        }
    }
    reverse_string(lcs);
    return lcs;
}

// Finds the longest common subsequence between the two provided strings
string lcs(string x, string y) {
    int len_x = x.length();
    int len_y = y.length();

    vector<vector<int>> table (
        len_x + 1,
        vector<int> (len_y + 1, 0)        
    );
    for (int i = 0; i <= len_x; i++) {
        for (int j = 0; j <= len_y; j++) {
            if (i == 0 || j == 0) {
                table[i][j] = 0;
            } else if (x[i - 1] == y[j - 1]) {
                table[i][j] = table[i - 1][j - 1] + 1;
            } else {
                table[i][j] = max(table[i - 1][j], table[i][j - 1]);
            }
        }
    }
    cout << show_table(table, x, y);
    return reconstruct_lcs(table, x, y);
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
