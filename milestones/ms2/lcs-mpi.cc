#include <mpi.h>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;


typedef struct {
    int rank;
    int start_x;
    int start_y;
    int end_x;
    int end_y;
} SectionInfo;

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

int get_top_left(vector<vector<int>> table, int x, int y, int* top, int* left) {
    if (y == 0) {
        return top[x];
    }
    if (x == 0) {
        return left[y-1];
    }
    return table[y - 1][x - 1];
}

int calculate_cell(vector<vector<int>> table, int x, int y,
                   string a, string b, int* top, int* left) {
    int left_cell = x > 0 && y >= 0 ? table[y][x - 1] : left[y];
    int top_cell = x >= 0 && y > 0 ? table[y - 1][x] : top[x+1];
    int top_left = get_top_left(table, x, y, top, left);
    if (x >= 0 && y >= 0 && a[y] == b[x]) {
        return top_left + 1;
    } else {
        return max(top_cell, left_cell);
    }
}

void diagonal_lcs(vector<vector<int>> table, string a, string b, int* top, int* left) {
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

int* get_top(vector<vector<int>> table, SectionInfo info) {
    int top_size = info.end_x - info.start_x + 2;
    int* top = (int*) calloc(top_size, sizeof(int));

    if (info.start_x == 0 || info.start_y == 0) {
        top[0] = 0;
    } else {
        top[0] = table[info.start_y - 1][info.start_x - 1];
    }

    if (!(info.start_y == 0)) {
        for (int i = 1; i < top_size; i++) {
            top[i] = table[info.start_y - 1][info.start_x + i];
        }
    }

    return top;
}

int* get_left(vector<vector<int>> table, SectionInfo info) {
    int left_size = info.end_y - info.start_y + 1;
    int* left = (int*) calloc(left_size, sizeof(int));
    if (!(info.start_x == 0)) {
        for (int i = 0; i < left_size; i++) {
            left[i] = table[info.start_y + i][info.start_x - 1];
        }
    }
    return left;
}

int* extract_solution(vector<vector<int>> table) {
    vector<int> single_dim = malloc(sizeof(int) * (table.size() * table[0].size()));
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
        for (int x = start_x, y = start_y, rank = 1; x <= limit_x; x++, y--, rank++) {
            SectionInfo info;
            info.rank = i == 0 ? 0 : rank;
            info.start_x = max(0, x * norm_width);
            info.start_y = max(0, y * norm_height);
            info.end_x = min((x + 1) * norm_width - 1, total_width - 1);
            info.end_y = min((y + 1) * norm_height - 1, total_height - 1);
            cout << "[x:" << x << ", y:" << y << ", r:" << info.rank << "]" << endl;
            cout << "(sx: " << info.start_x << ", ex: " << info.end_x << 
                        " - sy: " << info.start_y << ", ey: " << info.end_y <<
                ")" << endl;
            diagonal.push_back(info);
        }
        sections.push_back(diagonal);
    }
    return sections;
} 

void repopulate(vector<vector<int>> table, SectionInfo info, int* contents) {
    int index = 0;
    for (int y = info.start_y; y <= info.end_y; y++) {
        for (int x = info.start_x; x < info.end_x; x++, index++) {
            table[y][x] = contents[index];
        }
    }
}

// Reverses a string in place
void reverse_string(string& x) {
    reverse(x.begin(), x.end());
}

// Takes two strings and an LCS table and reconstructs the actual LCS string
string reconstruct_lcs(vector<vector<int>> table, string x, string y) {
    string lcs = "";
    int i = x.length() - 1;
    int j = y.length() - 1;
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

void lcs_parallel(string a, string b) {
    MPI_Init(NULL, NULL);

    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    num_procs = 3;
    vector<int> dims = get_mpi_dimensions(num_procs, a, b);
    int my_mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_rank);
    vector<vector<int>> table = construct_table(a.length(), b.length());
    
    // Root Process
    if (my_mpi_rank == 0) {
        vector<vector<SectionInfo>> sections = produce_sections(a, b, dims);
        // Process first section

        for (int diagonal = 1; diagonal < (int) sections.size() - 1; diagonal++) {
            for (int index = 0, rank = 1; index < (int) sections[diagonal].size() - 1; index++, rank++) {
                // Send top
                SectionInfo section = sections[diagonal][index];
                int* top = get_top(table, section);
                int section_data[4] = {section.start_x, section.end_x, section.start_y, section.end_y};
                    // Send info
                MPI_Send(
                    *section_data,
                    4,
                    MPI_INT,
                    rank,
                    0, // Tag
                    MPI_COMM_WORLD);
                
                    // Send top data
                MPI_Send(
                    top,
                    size,
                    MPI_INT,
                    rank,
                    0, // Tag
                    MPI_COMM_WORLD);

                // Send left
                int* left = get_left(table, SectionInfo);

                    // Send left data
                MPI_Send(
                    left,
                    size,
                    MPI_INT,
                    rank,
                    0, // Tag
                    MPI_COMM_WORLD);
            }

            for (int index = 0, rank = 1; index < (int) sections[diagonal].size() - 1; index++, rank++) {
                // Receive 
                SectionInfo section = sections[diagonal][index];
                int* section_content = malloc(sizeof(int) * size)
                MPI_Recv(
                    section_content,
                    size,
                    MPI_INT,
                    rank,
                    0,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE);
                // Populate table
                repopulate(table, section, section_content);
            }
        }
        // Process last section
        SectionInfo last = sections[sections.size() - 1][0];

        // Backtrace

        // Write out

        // Finalise
        MPI_Finalize();
        exit(0);
    } else {
    // Worker processes
        // Receive Section Info
        int section[4];
        MPI_Recv(
            *section,
            4,
            MPI_INT,
            0,
            0,
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);

        int start_x = section[0];
        int end_x = section[1];
        int start_y = section[2];
        int end_y = section[3];

        int top_size = (end_x - start_x + 2);
        int left_size = (end_y - start_y + 1);

        // Receive Top Data
        int* top = malloc(sizeof(int) * top_size);
        MPI_Recv(
            top,
            top_size,
            MPI_INT,
            0,
            0,
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);

        // Receive Left Data
        int* left = malloc(sizeof(int) * left_size);
        MPI_Recv(
            left,
            left_size,
            MPI_INT,
            0,
            0,
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);

        string a_sub = a.substr(start_y, left_size);
        string b_sub = b.substr(start_x, top_size - 1);
        vector<vector<int>> table_sub = construct_table(a_sub.length(), b_sub.lenth());

        int* result = diagonal_lcs(table, a_sub, b_sub, top, left);

       // Pass back to root
       MPI_Send(result,
                (int) (a.length() * b.length()),
                MPI_INT,
                rank,
                0, // Tag
                MPI_COMM_WORLD);
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
