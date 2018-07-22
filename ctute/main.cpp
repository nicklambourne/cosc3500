// Created by Nicholas Lambourne on 22/7/18.

#include "stdio.h"
#include "stdlib.h"

int populate_grid(FILE* file, char** grid, int height, int width) {
    for (int i = 0; i < height; i++) {
        char* row = (char*) malloc(sizeof(char) * width);
        for (int j = 0; j < width; j++) {
            row[j] = fgetc(file);
        }
        if (fgetc(file) != '\n') {
            return 1;
        }
        grid[i] = row;
    }
    return 0;

}

int main (int argc, char** argv) {
    if (argc != 4) {
        printf("Usage: match height width filename\n");
        return 1;
    }
    int height;
    int width;
    if (sscanf(argv[1], "%d", &height) == EOF ||
        sscanf(argv[2], "%d", &width) == EOF) {
        printf("Invalid grid dimensions\n");
        return 2;
    }
    FILE* grid_file = fopen(argv[3], "r");
    if (grid_file == NULL) {
        printf("Invalid grid file");
        return 3;
    }
    char** grid = (char**) malloc(sizeof(char*) * height);
    if (populate_grid(grid_file, grid, height, width) != 0) {
        printf("Error reading grid contents");
        return 4;
    }
    printf("Hello, world!\n");
    return 0;
}
