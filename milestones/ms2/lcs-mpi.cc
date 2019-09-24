#include <mpi.h>
#include <stdio.h>
#include <fstream>
#include <vector>


using namespace std;

vector<vector<int>> construct_table(int height, int width) {
    vector<vector<int>> table (
        height,
        vector<int> (width, 0)        
    );
    return table;
}


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

int get_top_left(vector<vector<int>> table, int x, int y, vector<int> top, vector<int> left) {
    if (y == 0) {
        return top[x];
    }
    if (x == 0) {
        return left[y-1];
    }
    return table[y - 1][x - 1];
}

int calculate_cell(vector<vector<int>> table, int x, int y,
                   string a, string b, vector<int> top, vector<int> left) {
    int left_cell = x > 0 && y >= 0 ? table[y][x - 1] : left[y];
    int top_cell = x >= 0 && y > 0 ? table[y - 1][x] : top[x+1];
    int top_left = get_top_left(table, x, y, top, left);
    if (x >= 0 && y >= 0 && a[y] == b[x]) {
        return top_left + 1;
    } else {
        return max(top_cell, left_cell);
    }
}

void diagonal_lcs(vector<vector<int>> table, string a, string b, vector<int> top, vector<int> left) {
    int height = a.length();
    int width = b.length();
    vector<vector<vector<int>>> indices;
    for (int i = 0; i < height + width - 1; i++) {
        int start_x = get_x(i, height);
        int start_y = get_y(i, height);
        int limit_x = end_x(start_x, start_y, width - 1);
        vector<vector<int>> diagonal;
        for (int x = start_x, y = start_y; x <= limit_x; x++, y--) {
            diagonal.push_back(vector<int>{x, y});
        }
        indices.push_back(diagonal);
    }
    
    for (int i = 0; i < (int) indices.size(); i++) {
        vector<vector<int>> diagonal = indices[i];
        #pragma omp parallel for shared(table, diagonal, a, b, top, left)
        for (int j = 0; j < (int) diagonal.size(); j++) {
            vector<int> cell = diagonal[j];
            int x = cell[0];
            int y = cell[1];
            int result = calculate_cell(table, x, y, a, b, top, left);
            table[y][x] = result;
        }
    }

    print_table(table, a, b);
}

vector<int> extract_solution(vector<vector<int>> table) {
    vector<int> single_dim;
    for (int y = 0; y < (int) table.size(); y++) {
        for (int x = 0; x < (int) table[y].size(); x++) {
            single_dim.push_back(table[y][x]);
        }
    }
    return single_dim;
}

void process_block(vector<vector<int>> table, vector<int> top, vector<int> left, string a, string b) {
    cout << "x" << endl; 
}

vector<string> get_substrings(string str, int size) {
    vector<string> substrings;
    for (int i = 0; i < (int) str.length(); i += size) {
        int limit = min(str.length() - i, size);
        substrings.push_back(str.substr(i, limit));
    }
    return substrings;
}

void lcs_parallel(string a, string b) {
    // vector<vector<int>> table = construct_table(a.length(), b.length());
    // vector<int> top (a.length() + 1, 0);
    // vector<int> left (b.length(), 0);

    //diagonal_lcs(table, a, b, top, left);
    vector<string> substrings = get_substrings(a, 3);
    for (int i = 0; i < (int) substrings.size(); i++) {
        cout << substrings[i] << endl;
    }
}

vector<int> get_mpi_dimensions(int num_procs) {
    if (num_procs < 3) {
        return vector<int>{1, 1};
    }
    
    int root = 1;

    int num_workers = num_procs - root;

    return vector<int>{num_workers, num_workers};

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

    lcs_parallel(string_a, string_b);

}
