#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<string>

using namespace std;


int main(int argc, char** argv) {
    if (argc != 2) {
        cout << "Invalid number of arguments!\n";
        return 1;
    }
    
    ifstream inFile;

    inFile.open(argv[1]);

    string buffer;

    if (!inFile) {
        cout << "Error opening file.\n";
        return 2;
    }

    int sum = 0;
    int count = 0;
    float average;

    for (int i = 0; i < 2; i++) {
        int dummy, good;
        inFile >> dummy;
        inFile.get();
        inFile >> good;
        sum += good;
        count++;
        getline(inFile, buffer, '\n');
    }
    average = sum / count;

    cout << average << "\n";
}
