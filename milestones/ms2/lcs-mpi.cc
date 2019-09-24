#include <mpi.h>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <math.h>


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
        int limit = min((int) str.length() - i, size);
        substrings.push_back(str.substr(i, limit));
    }
    return substrings;
}

vector<int> get_mpi_dimensions(int num_procs, string a, string b) {
     // Height, Width
    vector<int> dims = {1, 1};
    if (num_procs < 3) {
        return dims;
    }
    
    int root = 1;
    int num_workers = num_procs - root;

    if (!(a.length() == 1)) {
        dims[0] = num_workers;
    }

    if (!(b.length() == 1)) {
        dims[1] = num_workers;
    }

    return dims;
}

typedef struct {
    int id;
    int start_x;
    int start_y;
    int end_x;
    int end_y;
} SectionInfo;

vector<vector<SectionInfo>> produce_sections(string a, string b, vector<int> mpi_dims) {
    int cell_width = mpi_dims[1];
    int cell_height = mpi_dims[0];
    int total_width = a.length();
    int total_height = b.length();
    int norm_width = ceil(total_width / (float) cell_width);
    int norm_height = ceil(total_height / (float) cell_height);
    cout << norm_height << endl;
    vector<vector<SectionInfo>> sections;
    for (int i = 0; i < cell_width + cell_height - 1; i++) {
        int start_x = get_x(i, cell_height);
        int start_y = get_y(i, cell_height);
        int limit_x = end_x(start_x, start_y, cell_width - 1);
        vector<SectionInfo> diagonal;
        for (int x = start_x, y = start_y; x <= limit_x; x++, y--) {
            SectionInfo info;
            info.id = i;
            info.start_x = max(0, x * norm_width);
            info.start_y = max(0, y * norm_height);
            info.end_x = min((x + 1) * norm_width - 1, total_width - 1);
            info.end_y = min((y + 1) * norm_height - 1, total_height - 1);
            cout << "[x:" << x << ", y:" << y << "]" << endl;
            cout << "(sx: " << info.start_x << ", ex: " << info.end_x << 
                        " - sy: " << info.start_y << ", ey: " << info.end_y <<
                ")" << endl;
            diagonal.push_back(info);
        }
        sections.push_back(diagonal);
    }
    return sections;
} 

void lcs_parallel(string a, string b) {
    // vector<vector<int>> table = construct_table(a.length(), b.length());
    // vector<int> top (a.length() + 1, 0);
    // vector<int> left (b.length(), 0);

    //diagonal_lcs(table, a, b, top, left);
    MPI_Init(NULL, NULL);

    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    num_procs = 3;
    vector<int> dims = get_mpi_dimensions(num_procs, a, b);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_rank);
    
    // Root Process
    if (my_mpi_rank == 0) {
        vector<vector<SectionInfo>> sections = produce_sections(a, b, dims);
    } else {
        // Worker processes
    }

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
