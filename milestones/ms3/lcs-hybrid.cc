#include <mpi.h>
#include <stdio.h>
#include <fstream>
#include <vector>
#include<string>
#include<iostream> 
#include<algorithm> 
#include <math.h>

using namespace std;

// Information pertaining to the dimensions and location of a subsection.
typedef struct {
    int rank;
    int index;
    int start_x;
    int start_y;
    int end_x;
    int end_y;
} SectionInfo;

/**
 * Creates a two dimentional vector of specified size and populates it with
 * zeros.
 * @param height - target height of table
 * @param width - target width of table
 * @return - two-dimensional vector (table)
*/
vector<vector<int>> construct_table(int height, int width) {
    vector<vector<int>> table (
        height,
        vector<int> (width, 0)        
    );
    return table;
}

/** 
 * Left pads a string with spaces to reach a total specified length
 * @param contents - the string to pad
 * @param length - the target length of the final string
 * @return - a new string of the target length
 */
string pad(string contents, int length) {
    string normalised_string = contents;
    while ((int) normalised_string.size() < length) {
        normalised_string = " " + normalised_string;
    }
    return normalised_string;
}

/**
 * Prints the given table to stdout with row/col indexes and relevant
 * letters of each string.
 * @param table - the two-dimensional vector to print
 * @param a - the string represented on the vertical axis
 * @param b - the string represented on the horizontal axis
 */
void print_table(vector<vector<int>> &table, string a, string b) {
    cout << "     ";
    for (int i = 0; i < (int) b.length(); i++) {
        cout << pad(to_string(i), 3);
    }
    cout << endl << "     ";
    for (int i = 0; i < (int) b.length(); i++) {
        cout << pad(b.substr(i, 1), 3);
    }
    cout << endl;
    for (int y = 0; y < (int) table.size(); y++) {
        cout << pad(to_string(y), 3) << " " << a.substr(y, 1);
        for (int x = 0; x < (int) table[y].size(); x++) {
            cout << pad(to_string(table[y][x]), 3);
        }
        cout << endl;
    }
}

/**
 * Finds the x-coordinate for the last cell in an (anti-) diagonal.
 * @param x - the x-coordinate of the starting cell of the diagonal
 * @param y - the y-coordinate of the starting cell of the diagonal
 * @param width - the total width of the table
 * @return - the x-coordinate of the last cell in the diagonal
 */
int end_x(int start_x, int start_y, int width) {
    return min(start_x + start_y, width);
}   

/**
 * Provides the x-coordinate for a cell in a diagonal based on index.
 * @param index - the index number of the cell
 * @param height - the total height of the table
 */
int get_x(int index, int height) {
    if (index < height) {
        return 0;
    } else {
        return index - (height - 1);
    }
}

/**
 * Provides the y-coordinate for a cell in a diagonal based on index.
 * @param index - the index number of the cell
 * @param height - the total height of the table
 */
int get_y(int index, int height) {
    if (index >= height) {
        return height - 1;
    } else {
        return index;
    }
}

/**
 * Gets the value in the cell to top left of the given coordinates.
 * @param table - table in which the cells reside
 * @param x - x-coordinate of original cell
 * @param y - y-coordinate of original cell
 * @param top - the row of values to the top of the specified coordinates
 * @param left - the column of values to the left of the specified coordinates
 * @return - the value of the cell to the top left of the given coordinates
 *           or 0, if unreachable.
 */
int get_top_left(vector<vector<int>> &table, int x, int y, int* top, int* left) {
    if (y == 0) {
        return top[x];
    }
    if (x == 0) {
        return left[y-1];
    }
    return table[y - 1][x - 1];
}

/**
 * Calculates the value for a cell in the table, based on max of the top and 
 * left adjacent cells. If the corresponding letters of a and b match at those 
 * indices, instead return the value of the top-left diagonally adjacent cell
 * plus one.
 * @param x - x-coordinate of target cell
 * @param y - y-coordinate of target cell
 * @param top - the row of values to the top of the specified coordinates
 * @param left - the column of values to the left of the specified coordinates
 * @return - the calculated value to populate the target cell with
 */
int calculate_cell(vector<vector<int>> &table, int x, int y,
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

void diagonal_lcs(vector<vector<int>> &table, string a, string b, 
                  int* top, int* left) {
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
}

int* get_top(vector<vector<int>> &table, SectionInfo info) {
    int top_size = info.end_x - info.start_x + 2;
    int* top = (int*) calloc(top_size, sizeof(int));

    if (info.start_x == 0 || info.start_y == 0) {
        top[0] = 0;
    } else {
        top[0] = table[info.start_y - 1][info.start_x - 1];
    }

    if (!(info.start_y == 0)) {
        for (int i = 1, x_index = 0; i < top_size; i++, x_index++) {
            top[i] = table[info.start_y - 1][info.start_x + x_index];
        }
    }

    return top;
}

int* get_left(vector<vector<int>> &table, SectionInfo info) {
    int left_size = info.end_y - info.start_y + 1;
    int* left = (int*) calloc(left_size, sizeof(int));
    if (!(info.start_x == 0)) {
        for (int i = 0; i < left_size; i++) {
            left[i] = table[info.start_y + i][info.start_x - 1];
        }
    }
    return left;
}

int* extract_solution(vector<vector<int>> &table) {
    int size = (table.size() * table[0].size());
    int* single_dim = (int*) malloc(sizeof(int) * size);
    int index = 0;
    for (int y = 0; y < (int) table.size(); y++) {
        for (int x = 0; x < (int) table[y].size(); x++, index++) {
            single_dim[index] = table[y][x];
        }
    }
    return single_dim;
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
    int total_width = b.length();
    int total_height = a.length();
    int norm_width = ceil(total_width / (float) cell_width);
    int norm_height = ceil(total_height / (float) cell_height);
    vector<vector<SectionInfo>> sections;
    for (int i = 0, index = 0; i < cell_width + cell_height - 1; i++) {
        int start_x = get_x(i, cell_height);
        int start_y = get_y(i, cell_height);
        int limit_x = end_x(start_x, start_y, cell_width - 1);
        vector<SectionInfo> diagonal;
        for (int x = start_x, y = start_y, rank = 1; x <= limit_x; x++, y--, rank++, index++) {
            SectionInfo info;
            info.rank = i == 0 ? 0 : rank;
            info.index = index;
            info.start_x = max(0, x * norm_width);
            info.start_y = max(0, y * norm_height);
            info.end_x = min((x + 1) * norm_width - 1, total_width - 1);
            info.end_y = min((y + 1) * norm_height - 1, total_height - 1);
            // cout << "[i: " << index << ", x:" << x << ", y:" << y << ", r:" << info.rank << "]" << endl;
            // cout << "(sx: " << info.start_x << ", ex: " << info.end_x << 
            //             " - sy: " << info.start_y << ", ey: " << info.end_y <<
            //     ")" << endl;
            diagonal.push_back(info);
        }
        sections.push_back(diagonal);
    }
    return sections;
} 

// Reverses a string in place
void reverse_string(string& x) {
    reverse(x.begin(), x.end());
}

// Takes two strings and an LCS table and reconstructs the actual LCS string
string reconstruct_lcs(vector<vector<int>> &table, string a, string b) {
    string lcs = "";
    int y = a.length() - 1;
    int x = b.length() - 1;
    // print_table(table, a, b);
    while (x >= 0 && y >= 0) {
        if (a[y] == b[x]) {  // Characters match
            lcs.append(a, y, 1);
            y = max(y - 1, -1);
            x = max(x - 1, -1);
        } else if (y == 0) {  // Cannot go any further left
            x--;
        } else if (x == 0) {  // Cannot go any further up
            y--;
        } else if (table[y][x - 1] > table[y - 1][x]) {  // Left is larger 
            x--;
        } else { // Up is larger or they are the same
            y--;
        }
    }
    reverse_string(lcs);
    return lcs;
}

void repopulate(vector<vector<int>> &table, SectionInfo info, int* contents, string a, string b) {
    for (int y = info.start_y, index = 0; y <= info.end_y; y++) {
        for (int x = info.start_x; x <= info.end_x; x++, index++) {
            table[y][x] = contents[index];
        }
    }
}

void send_section_info(SectionInfo section, int rank) {
    int section_data[5] = {section.start_x, section.end_x, section.start_y, section.end_y, section.index};
    // Send section info
    MPI_Send(
        &section_data,
        5,
        MPI_INT,
        rank,
        0, // Tag
        MPI_COMM_WORLD);
}

void send_top(SectionInfo section, vector<vector<int>> &table, int rank) {
    // Send top data
    int* top = get_top(table, section);
    int top_size = (section.end_x - section.start_x + 2);
    MPI_Send(
        top,
        top_size,
        MPI_INT,
        rank,
        0, // Tag
        MPI_COMM_WORLD);
}

void send_left(SectionInfo section, vector<vector<int>> &table, int rank) {
    // Send left
    int* left = get_left(table, section);
    int left_size = (section.end_y - section.start_y + 1);
    // Send left data
    MPI_Send(
        left,
        left_size,
        MPI_INT,
        rank,
        0, // Tag
        MPI_COMM_WORLD);
}

void send_close(int rank) {
    int close_data[4] = {-1, -1, -1, -1};
    // Send close
    MPI_Send(
        &close_data,
        4,
        MPI_INT,
        rank,
        0, // Tag
        MPI_COMM_WORLD);
}

void send_section(int* section_data, int size) {
    // Pass back to root
    MPI_Send(
        section_data,
        size,
        MPI_INT,
        0, // Destination
        0, // Tag
        MPI_COMM_WORLD);
}

int* receive_section(SectionInfo section, int rank) {
    int section_size = (section.end_x - section.start_x + 1) * (section.end_y - section.start_y + 1);
    int* section_content  = (int*) malloc(sizeof(int) * section_size);
    MPI_Recv(
        section_content,
        section_size,
        MPI_INT,
        rank,
        0,
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE);
    return section_content;
}

int* receive_border(int size) {
    int* border = (int*) malloc(sizeof(int) * size);
    MPI_Recv(
        border,
        size,
        MPI_INT,
        0,
        0,
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE);
    return border;
}

SectionInfo receive_info(int rank) {
    int section[5];
    MPI_Recv(
        &section,
        5,
        MPI_INT,
        0,
        0,
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE);
    

    // cout << "rank: " << rank << " index: " << section[4] << endl;
    // cout << "sx: " << section[0] << ", ex: " << section[1] << endl;
    // cout << "sy: " << section[2] << ", ey: " << section[3] << endl << endl;

    SectionInfo info;

    info.start_x = section[0];
    info.end_x = section[1];
    info.start_y = section[2];
    info.end_y = section[3];
    info.index = section[4];
    info.rank = rank;
    
    return info;
}

void process_section(vector<vector<int>> &table, SectionInfo info, string a, string b) {
    int top_size = (info.end_x - info.start_x + 2);
    int left_size = (info.end_y - info.start_y + 1);
    int end_a = left_size;
    int end_b = top_size - 1;
    string a_sub = a.substr(info.start_y, end_a);
    string b_sub = b.substr(info.start_x, end_b);
    int* top = get_top(table, info);
    int* left = get_left(table, info);
    vector<vector<int>> section_table = construct_table(a_sub.length(), b_sub.length());

    diagonal_lcs(section_table, a_sub, b_sub, top, left);

    // print_table(section_table, a_sub, b_sub);

    int* result = extract_solution(section_table);

    repopulate(table, info, result, a, b);
}

int* process_worker_section(SectionInfo info, string a, string b, int* top, int* left, int section_height, int section_width) {
    string a_sub = a.substr(info.start_y, section_height);
    string b_sub = b.substr(info.start_x, section_width);

    vector<vector<int>> section_table = construct_table(section_height, section_width);

    diagonal_lcs(section_table, a_sub, b_sub, top, left);

    // print_table(section_table, a_sub, b_sub);

    return extract_solution(section_table);
}

string lcs_parallel(string a, string b, string outfile) {
    if (a.length() == 0 || b.length() == 0) {
        return "";
    }
    
    MPI_Init(NULL, NULL);

    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    vector<int> dims = get_mpi_dimensions(num_procs, a, b);
    int my_mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_rank);
    vector<vector<int>> table = construct_table(a.length(), b.length());
    
    // Root Process
    if (my_mpi_rank == 0) {
        vector<vector<SectionInfo>> sections = produce_sections(a, b, dims);
        // Process first section
        process_section(table, sections[0][0], a, b);

        for (int diagonal = 1; diagonal < (int) sections.size() - 1; diagonal++) {

            for (int index = 0, rank = 1; index < (int) sections[diagonal].size(); index++, rank++) {
                // Send top
                SectionInfo section = sections[diagonal][index];
                
                send_section_info(section, rank);
                send_top(section, table, rank);
                send_left(section, table, rank);
            }

            for (int index = 0, rank = 1; index < (int) sections[diagonal].size(); index++, rank++) {

                // Receive 
                SectionInfo section = sections[diagonal][index];
                int* section_content = receive_section(section, rank);
        
                repopulate(table, section, section_content, a, b);

            }
        }

        for (int rank = 1; rank < num_procs; rank++) {
            send_close(rank);
        }

        // Process last section
        if (sections.size() > 1) {
            SectionInfo last = sections[sections.size() - 1][0];
            process_section(table, last, a, b);
        }

        // print_table(table, a, b);
        
        // Backtrace
        string lcs = reconstruct_lcs(table, a, b);

        MPI_Finalize();
        return lcs;

    } else {
    // Worker processes
        while (1) {
            SectionInfo section = receive_info(my_mpi_rank);

            if (section.start_x == -1) {
                MPI_Finalize();
                exit(0);
            }

            int section_height = section.end_y - section.start_y + 1;
            int section_width = section.end_x - section.start_x + 1;

            int* top = receive_border(section_width + 1);
            int* left = receive_border(section_height);

            int* result = process_worker_section(section, a, b, top, left, section_height, section_width);

            int size = section_height * section_width;

            send_section(result, size);
        }
        
    }

}

/**
 * Gets the input strings from the specified file and populates the two given
 * strings.
 * @param input_file_path - file to read strings from
 * @param a - pointer to first string to populate
 * @param b - pointer to second string to populate
 */ 
void gather_strings(string input_file_path, string* a, string* b) {
    // Read in sequences to compare
    fstream input_file;
    input_file.open(input_file_path);

    if (input_file.is_open()) {
        string buffer;
        if (!getline(input_file, *a) ||
            !getline(input_file, *b)) {
            cerr << "Error reading from provided file." << endl;
            exit(3);
        }
    } else {
        cerr << "Could not open provided file." << endl;
        exit(2);
    }
}

/**
 * Writes the given LCS out to the file specified and prints it to stdout
 * @param lcs - the string to write
 * @param output_file_path - the path of the file to write to
 */ 
void write_lcs(string lcs, string output_file_path) {
    ofstream output_file;
    output_file.open(output_file_path);

    if (output_file.is_open()) {
        output_file << lcs;
    } else {
        cerr << "Could not write to given output file." << endl;
        exit(4);
    }

    cout << lcs << endl;
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

    string string_a, string_b;
    gather_strings(file_name, &string_a, &string_b);

    string lcs = lcs_parallel(string_a, string_b, output_file_name);

    write_lcs(lcs, output_file_name);
}
