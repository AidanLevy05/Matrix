#include "Matrix.h"
using namespace std;

int main()
{

  srand(time(0));

  string file;
  cout << "Please enter a file name: " << endl;
  cin >> file;

  Matrix matrix1(file);
  vector<double> solutions;

  try
  {
    cout << "Before RREF:" << endl
         << endl;
    matrix1.display();
    cout << endl;

    solutions = matrix1.rref();

    cout << "After RREF:" << endl
         << endl;
    matrix1.display();
    cout << endl;
  }
  catch (const exception &e)
  {
    cout << e.what() << endl;
  }
}
