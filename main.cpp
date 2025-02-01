#include "Matrix.h"
#include <cstdlib>
#include <string>

void div();
void div(string);
void showAllDimensions(Matrix&, Matrix&, Matrix&, char);

int main() {

    srand(time(0));
    int maxRand = 10;

    div("Multpilying 2 matricies");
    int r = 5;
    int c = 2;
    Matrix matrix1(r, c);
    matrix1.setRandom(maxRand);
    matrix1.display();
    div();

    int r2 = c;
    int c2 = 3;
    Matrix matrix2(r2, c2);
    matrix2.setRandom(maxRand);
    matrix2.display();
    div();

    Matrix matrix3;
    matrix3 = matrix1 * matrix2;
    matrix3.display();
    div();

    showAllDimensions(matrix1, matrix2, matrix3, 'x');
    div("Adding two matricies");

    Matrix matrix4(3, 3);
    matrix4.setRandom(maxRand);
    matrix4.display();
    div();

    Matrix matrix5(3, 3);
    matrix5.setRandom(maxRand);
    matrix5.display();
    div();

    Matrix matrix6 = matrix4 + matrix5;
    matrix6.display();
    div();

    showAllDimensions(matrix4, matrix5, matrix6, '+');
    div("Subtracting two matricies");

    matrix4.display();
    div();
    matrix5.display();
    div();
    matrix6 = matrix4-matrix5;
    matrix6.display();
    div();

    showAllDimensions(matrix4, matrix5, matrix6, '-');
    div("Pulling matrix from file & ensuring that solving works");

    string file = "matrix.txt";
    Matrix matrix7(file);
    matrix7.display();

    int sols[3] {2, 3, 4};
    if (matrix7.solve(sols, 3)) {
        cout << "This solution set works" << endl;
    } else {
        cout << "This solution set does not work" << endl;
    }

    div("Row swap on last matrix");
    matrix7.swapRows(1, 2);
    matrix7.display();

    div("Row multiplication on last row");
    matrix7.multiplyRow(2, -1);
    matrix7.display();

    div("Returning to normal...");
    matrix7.multiplyRow(2, -1);
    matrix7.swapRows(1, 2);
    matrix7.display();

    div("Row replacement on last row with first row, scalar = -5");
    matrix7.replaceRow(2, 0, -5);
    matrix7.display();
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
