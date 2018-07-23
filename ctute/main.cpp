// Created by Nicholas Lambourne on 22/7/18.

#include "stdio.h"
#include "stdlib.h"


int parse_arguments(int* height, int* width, FILE* grid_file) {
    if (sscanf(argv[1], "%d", &height) == EOF ||
        sscanf(argv[2], "%d", &width) == EOF) {
        printf("Invalid grid dimensions\n");
        return 2;
    }
    grid_file = fopen(argv[3], "r");
    if (grid_file == NULL) {
        printf("Invalid grid file\n");
        return 3;
    }
    return 0;
}

void free_grid(char** grid, int height) {
    for (int i = 0; i < height; i++) {
        free(grid[i]);
    }
}

int populate_grid(FILE* file, char** grid, int height, int width) {
    for (int i = 0; i < height; i++) {
        char* row = (char*) malloc(sizeof(char) * width);
        for (int j = 0; j < width; j++) {
            row[j] = fgetc(file);
        }
        if (i != height - 1 && fgetc(file) != '\n') {
            return 1;
        }
        grid[i] = row;
    }
    return 0;
}

int write_grid_to_file(FILE* file, char** grid) {

}

int main (int argc, char** argv) {
    if (argc != 4) {
        printf("Usage: match height width filename\n");
        return 1;
    }
    int height, width;
    FILE* grid_file;
    args = parse_arguments(&height, &width, grid_file);
    if (args != 0) {
        return args;
    }
    char** grid = (char**) malloc(sizeof(char*) * height);
    if (populate_grid(grid_file, grid, height, width) != 0) {
        printf("Error reading grid contents\n");
        return 4;
    }
    printf("Hello, world!\n");
    free_grid(grid, height);
    return 0;
}
