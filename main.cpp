#include "Matrix.h"
#include <cstdlib>
#include <string>

void div();
void div(string);
void showAllDimensions(Matrix&, Matrix&, Matrix&, char);

int main() {

    srand(time(0));
    string file = "matrix.txt";

    // populate matrix.txt
    const int ROWS = 10;
    const int COLS = ROWS;
    const int RAND = 1000;
    ofstream outFile;
    outFile.open(file);

    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            outFile << rand() % RAND << " ";
        }
        outFile << endl;
    }

    outFile.close();


    Matrix matrix(file);

    div("Before RREF");
    matrix.display();

    try {
        vector<double> sols = matrix.rref();
        int numElements = sols.size();
        div("After RREF");
        matrix.display();
        cout << "Solutions work: " << matrix.solve(sols, numElements) << endl;
        matrix.parametricVectorForm();
    } catch (const exception& e) {
        cerr << e.what() << endl;
    }

    div();

    return 0;
}

void div() {
    cout << endl;
    for (int i = 0; i < 30; i++) cout << "-";
    cout << endl << endl;
}

void div(string text) {
    cout << endl;
    cout << text << "  ";
    for (int i = 0; i < 30; i++) cout << "-";
    cout << endl << endl;
}

void showAllDimensions(Matrix& m1, Matrix& m2, Matrix& m3, char character) {
    cout << "Matrix 1: ";
    m1.dimensions();
    cout << "\t\t" << character << endl;
    cout << "Matrix 2: ";
    m2.dimensions();
    cout << "\t\t=" << endl;
    cout << "Matrix 3: ";
    m3.dimensions();
}
