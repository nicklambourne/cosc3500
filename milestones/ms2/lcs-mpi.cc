#include <mpi.h>
#include <stdio.h>

using namespace std;

void gather_strings(string input_file, string* x, string* y) {
    // Read in sequences to compare
    ifstream input_file;
    input_file.open(file_name);

    if (input_file.is_open()) {
        string buffer;
        if (!getline(input_file, *x) ||
            !getline(input_file, *y)) {
            cerr << "Error reading from provided file." << endl;
            exit(3);
        }
    } else {
        cerr << "Could not open provided file." << endl;
        exit(2);
    }
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

    string string_x;
    string string_y;

    gather_strings(file_name, string_x, string_y);

    cout << string_x << endl;
    cout << string_y << endl;
}
