#include <mpi.h>
#include <stdio.h>
#include <fstream>
#include <vector>


using namespace std;


// Left pads a string with spaces to reach CELL_LENGTH
string pad(string contents, int length) {
    string normalised_string = contents;
    while ((int) normalised_string.size() < length) {
        normalised_string = " " + normalised_string;
    }
    return normalised_string;
}

void print_table(vector<vector<int>> table, string a, string b) {
    cout << " ";
    for (int i = 0; i < (int) b.length(); i++) {
        cout << pad(b.substr(i, 1), 3);
    }
    cout << endl;
    for (int y = 0; y < (int) table.size(); y++) {
        cout << a.substr(y, 1);
        for (int x = 0; x < (int) table[y].size(); x++) {
            cout << pad(to_string(table[y][x]), 3);
        }
        cout << endl;
    }
}

int end_x(int x, int y, int width) {
    return min(x + y, width);
}   

int get_x(int index, int height) {
    if (index < height) {
        return 0;
    } else {
        return index - (height - 1);
    }
}

int get_y(int index, int height) {
    if (index >= height) {
        return height - 1;
    } else {
        return index;
    }
}

int calculate_cell(vector<vector<int>> table, int x, int y,
                   string a, string b) {
    int left = x > 0 && y >= 0 ? table[y][x - 1] : 0;
    int top = x >= 0 && y > 0 ? table[y - 1][x] : 0;
    int top_left = x > 0 && y > 0 ? table[y - 1][x - 1] : 0;
    if (x >= 0 && y >= 0 && a[y] == b[x]) {
        return top_left + 1;
    } else {
        return max(top, left);
    }
}

void diagonal_lcs(vector<vector<int>> table, string a, string b) {
    int height = a.length();
    int width = b.length();
    for (int i = 0; i < height + width - 1; i++) {
        int current_x = get_x(i, height);
        int current_y = get_y(i, height);
        int end = end_x(current_x, current_y, width - 1);
        cout << "start position - x: "  << current_x << ", y: " << current_y << endl;
        cout << "end_x: " << end << endl;
        print_table(table, a, b);
        
        for (; current_x <= end;) {
            int result = calculate_cell(table, current_x, current_y, a, b);
            cout << "(x: "  << current_x << ", y: " << current_y << ", result: " << result << ")"  << endl;
            table[current_y][current_x] = result;
            current_x++;
            current_y--;
        }
    }
    print_table(table, a, b);
}

vector<vector<int>> construct_table(string a, string b) {
    int height = a.length();
    int width = b.length();
    vector<vector<int>> table (
        height,
        vector<int> (width, 0)        
    );
    return table;
}

void gather_strings(string file_name, string* x, string* y) {
    // Read in sequences to compare
    fstream input_file;
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

    string string_a;
    string string_b;

    gather_strings(file_name, &string_a, &string_b);

    cout << string_a << endl;
    cout << string_b << endl;

    vector<vector<int>> table = construct_table(string_a, string_b);

    diagonal_lcs(table, string_a, string_b);

}
