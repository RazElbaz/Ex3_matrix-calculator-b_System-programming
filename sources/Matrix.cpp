#include <vector>
#include <iostream>
#include "Matrix.hpp"
using namespace std;

namespace zich {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adding operators
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Matrix Matrix::operator+(const Matrix &other) {
        //This operator can only accept matrices of this size
        if (this->Rows != other.Rows || this->Columns != other.Columns) {
            throw invalid_argument("The two matrices should have the same amount of columns and rows");
        }
        vector<double> ans;
        for (int i = 0; i < Rows; i++) {
            for (int j = 0; j < Columns; j++) {
                //Inserting into the vector ans the final answer
                ans.push_back((unsigned int) (this->Values[(double) ((i * Columns) + j)] +other.Values[(double) ((i * Columns) + j)]));
            }
        }
        //Return of a new matrix type object with the updated data
        Matrix matrix(ans, other.Rows, other.Columns);
        return matrix;
    }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Matrix &Matrix::operator+=(const Matrix &other) {
        //This operator can only accept matrices of this size
        if (other.Rows != this->Rows || other.Columns != this->Columns) {
            throw invalid_argument("The two matrices should have the same amount of columns and rows");
        }

        for (int i = 0; i < other.Rows; i++) {
            for (int j = 0; j < other.Columns; j++) {
                //Inserting into this matrix the final answer
                this->Values[(double) ((i * Columns) + j)] += other.Values[(double) ((i * Columns) + j)];
            }
        }
        //Return of pointer to the matrix on which we performed the operation
        return *this;
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Matrix operator+(Matrix &other) {
        vector<double> ans;
        for (int i = 0; i < other.Rows; i++) {
            for (int j = 0; j < other.Columns; j++) {
                //Inserting into the vector ans the final answer
                ans.push_back((double) other.Values[(double) ((i * other.Columns) + j)]);
            }
        }
        //Return of a new matrix type object with the updated data
        Matrix matrix(ans, other.Rows, other.Columns);
        return matrix;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Subtraction operator
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Matrix Matrix::operator-(const Matrix &other) {
        //This operator can only accept matrices of this size
        if (other.Rows != this->Rows || other.Columns != this->Columns) {
            throw invalid_argument("The two matrices should have the same amount of columns and rows");
        }
        vector<double> ans;
        for (int i = 0; i < other.Rows; i++) {
            for (int j = 0; j < other.Columns; j++) {
                //Inserting into the vector ans the final answer
                ans.push_back(this->Values[(double) ((i * Columns) + j)] - other.Values[(double) ((i * Columns) + j)]);
            }
        }
        //Return of a new matrix type object with the updated data
        Matrix matrix(ans, other.Rows, other.Columns);
        return matrix;
    }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Matrix &Matrix::operator-=(const Matrix &other) {
        //This operator can only accept matrices of this size
        if (other.Rows != this->Rows || other.Columns != this->Columns) {
            throw invalid_argument("The two matrices should have the same amount of columns and rows");
        }
        for (int i = 0; i < other.Rows; i++) {
            for (int j = 0; j < other.Columns; j++) {
                //Inserting into this matrix the final answer
                this->Values[(double) ((i * Columns) + j)] -= other.Values[(double) ((i * Columns) + j)];
            }
        }
        //Return of pointer to the matrix on which we performed the operation
        return *this;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Matrix operator-(Matrix &other) {
        vector<double> ans;
        for (int i = 0; i < other.Rows; i++) {
            for (int j = 0; j < other.Columns; j++) {
                //Inserting into the vector ans the final answer
                ans.push_back((-1) * other.Values[(double) ((i * other.Columns) + j)]);
            }
        }
        //Return of a new matrix type object with the updated data
        Matrix matrix(ans, other.Rows, other.Columns);
        return matrix;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Comparison operator
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool Matrix::operator>(Matrix &other) {
        //This operator can only accept matrices of this size
        if (other.Rows != this->Rows || other.Columns != this->Columns) {
            throw invalid_argument("The two matrices should have the same amount of columns and rows");
        }

        double ans = 0;
        double ans_other = 0;
        for (int i = 0; i < other.Rows; i++) {
            for (int j = 0; j < other.Columns; j++) {
                //Scheme of the variables in the current matrix and save
                ans = ans + this->Values[(double) ((i * Columns) + j)];
                //Scheme of the variables in the other matrix and save
                ans_other += other.Values[(double) ((i * Columns) + j)];
            }
        }
        return ans > ans_other;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool Matrix::operator>=(Matrix &other) {
        //This operator can only accept matrices of this size
        if (other.Rows != this->Rows || other.Columns != this->Columns) {
            throw invalid_argument("The two matrices should have the same amount of columns and rows");
        }
        double ans = 0;
        double ans_other = 0;
        for (int i = 0; i < other.Rows; i++) {
            for (int j = 0; j < other.Columns; j++) {
                //Scheme of the variables in the current matrix and save
                ans = ans + this->Values[(double) ((i * Columns) + j)];
                //Scheme of the variables in the other matrix and save
                ans_other += other.Values[(double) ((i * Columns) + j)];
            }
        }
        return ans >= ans_other;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool Matrix::operator<(Matrix &other) {
        //This operator can only accept matrices of this size
        if (other.Rows != this->Rows || other.Columns != this->Columns) {
            throw invalid_argument("The two matrices should have the same amount of columns and rows");
        }
        double ans = 0;
        double ans_other = 0;
        for (int i = 0; i < other.Rows; i++) {
            for (int j = 0; j < other.Columns; j++) {
                //Scheme of the variables in the current matrix and save
                ans = ans + this->Values[(double) ((i * Columns) + j)];
                //Scheme of the variables in the other matrix and save
                ans_other += other.Values[(double) ((i * Columns) + j)];
            }
        }
        return ans < ans_other;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool Matrix::operator<=(Matrix &other) {
        //This operator can only accept matrices of this size
        if (other.Rows != this->Rows || other.Columns != this->Columns) {
            throw invalid_argument("The two matrices should have the same amount of columns and rows");
        }
        double ans = 0;
        double ans_other = 0;
        for (int i = 0; i < other.Rows; i++) {
            for (int j = 0; j < other.Columns; j++) {
                //Scheme of the variables in the current matrix and save
                ans = ans + this->Values[(double) ((i * Columns) + j)];
                //Scheme of the variables in the other matrix and save
                ans_other += other.Values[(double) ((i * Columns) + j)];
            }
        }
        return ans <= ans_other;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool Matrix::operator==(const Matrix &other) const {
        //This operator can only accept matrices of this size
        if (other.Rows != Rows || other.Columns != Columns) {
            throw invalid_argument("The two matrices should have the same amount of columns and rows");
        }

        for (int i = 0; i < Rows; i++) {
            for (int j = 0; j < Columns; j++) {
                if (this->Values[(double) ((i * Columns) + j)] != other.Values[(double) ((i * Columns) + j)]) {
                    return false;
                }
            }
        }
        return true;
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    bool Matrix::operator!=(const Matrix &other) {
        //This operator can only accept matrices of this size
        if (other.Rows != Rows || other.Columns != Columns) {
            throw invalid_argument("The two matrices should have the same amount of columns and rows");
        }
        for (int i = 0; i < Rows; i++) {
            for (int j = 0; j < Columns; j++) {
                if (this->Values[(double) ((i * Columns) + j)] != other.Values[(double) ((i * Columns) + j)]) {
                    return true;
                }
            }
        }
        return false;
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Increasement operator
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ++i
    Matrix &Matrix::operator++() {
        for (int i = 0; i < this->Rows; i++) {
            for (int j = 0; j < this->Columns; j++) {
                this->Values[(double) ((i * Columns) + j)] += 1;
            }
        }
        //Return of pointer to the matrix on which we performed the operation
        return *this;
    }

    //i++
    Matrix Matrix::operator++(int) {
        //Create a new matrix object with data like that of the current matrix
        Matrix same = *this;
        for (int i = 0; i < this->Rows; i++) {
            for (int j = 0; j < this->Columns; j++) {
                this->Values[(double) ((i * Columns) + j)] += 1;
            }
        }
        return same;
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // --i
    Matrix &Matrix::operator--() {
        for (int i = 0; i < this->Rows; i++) {
            for (int j = 0; j < this->Columns; j++) {
                //Progress of each value in the matrix in 1
                this->Values[(double) ((i * Columns) + j)] -= 1;
            }
        }
        //Return of pointer to the matrix on which we performed the operation
        return *this;
    }

    //i--
    Matrix Matrix::operator--(int) {
        //Create a new matrix object with data like that of the current matrix
        Matrix same = *this;
        for (int i = 0; i < this->Rows; i++) {
            for (int j = 0; j < this->Columns; j++) {
                //Subtraction of any value in the matrix in 1
                this->Values[(double) ((i * Columns) + j)] -= 1;
            }
        }
        return same;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Scalar doubling
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Matrix operator*(const double scalar, Matrix &other) {
        //create a vector with the appropriate values
        vector<double> ans((unsigned long) (other.Rows * other.Columns), 0);
        //Go over each part of the matrix and multiply the variables
        for (int i = 0; i < other.Rows; i++) {
            for (int j = 0; j < other.Columns; j++) {
                ans[(unsigned int) ((i * other.Columns) + j)] = (other.Values[(double) ((i * other.Columns) + j)] *
                                                                 scalar);
            }
        }
        //Return of a new matrix type object with the updated data
        Matrix matrix(ans, other.Rows, other.Columns);
        return matrix;
    }

    Matrix operator*(Matrix &other, const double scalar) {
        //create a vector with the appropriate values
        vector<double> ans((unsigned long) (other.Rows * other.Columns), 0);
        //Go over each part of the matrix and multiply the variables
        for (int i = 0; i < other.Rows; i++) {
            for (int j = 0; j < other.Columns; j++) {
                ans[(unsigned int) ((i * other.Columns) + j)] = (other.Values[(double) ((i * other.Columns) + j)] *
                                                                 scalar);
            }
        }
        //Return of a new matrix type object with the updated data
        Matrix matrix(ans, other.Rows, other.Columns);
        return matrix;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Matrix &Matrix::operator*=(const Matrix &other) {
        //Since we ran over the * operator we can use it
        *this = *this * other;
        //Return of pointer to the matrix on which we performed the operation
        return *this;
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Multiplication of two matrices
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
    Matrix Matrix::operator*(const Matrix &other) {
        //In multiplication matrices the number of columns of the first matrix should be equal to the number of rows of the second matrix
        if (this->Columns != other.Rows) {
            throw invalid_argument("The two matrices should have the same amount of columns and rows");
        }
        // Multiplication matrix's dimension by definition.
        vector<double> add_mat((unsigned long) (Rows * other.Columns), 0);
        //Go over each part of the matrix
        for (unsigned long i = 0, m = 0; i < Rows; ++i) {
            for (unsigned long j = 0; j < other.Columns; ++j) {
                //saving current sum
                double sum = 0;
                for (unsigned long k = 0; k < Columns; ++k) {
                    sum += Values[i * (unsigned long) (Columns) + k] *
                           other.Values[k * (unsigned long) (other.Columns) + j];
                }
                add_mat.at(m++) = sum;
            }
        }
        //Return of a new matrix type object with the updated data
        Matrix add{add_mat, Rows, other.Columns};
        return add;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Matrix &Matrix::operator*=(const double scalar) {
        for (int i = 0; i < this->Rows; i++) {
            for (int j = 0; j < this->Columns; j++) {
                this->Values[(double) ((i * Columns) + j)] *= scalar;
            }
        }
        //Return of pointer to the matrix on which we performed the operation
        return *this;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// friend global IO operators
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ostream &operator<<(ostream &output, const Matrix &matrix) {
        //Go over each part of the matrix
        for (int i = 0; i < matrix.Rows; ++i) {
            output << "[";
            for (int j = 0; j < matrix.Columns; ++j) {
                output << matrix.Values[(double) ((i * matrix.Columns) + j)];
                //If we have reached the last column
                if (j < matrix.Columns - 1) {
                    output << " ";
                }
            }
            output << "]";
            //If we not have reached the last row
            if (i != (matrix.Rows - 1)) { output << "\n"; }
        }
        return output;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    vector <string> Split_by_character(string will_split, char character) {
        int length = will_split.length();
        vector <string> ans;
        string current;
        for (int i = 0; i < length; i++) {
            //we don't want ' ' in the string
            if (will_split[(unsigned long) i] == ' ') { break; }
            //we want to take numbers only so we don't want to take the characters of the structure of the matrix
            if (will_split[(unsigned long) i] != ']' && will_split[(unsigned long) i] != '[') {
                if (will_split[(unsigned long) i] == character) {
                    ans.push_back(current);
                    //Reset the sting to an empty sting
                    current = "";
                } else {
                    //Taking the numbers only
                    current += will_split[(unsigned long) i];
                }
            }
        }
        ans.push_back(current);
        return ans;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    istream &operator>>(istream &input, Matrix &matrix) {
        string input_matrix;
        char Input = '0';
        //Vector variables that will keep us the numbers of the matrix
        vector <string> Columns;
        vector <string> Rows;
        //The absorbed matrix will be retained in the answer variable
        vector<double> ans;

        //Absorb until the character arrives \ n
        while (Input != '\n') {
            Input = input.get();
            //write every char to the string
            input_matrix += Input;
        }
        unsigned long size = input_matrix.length();
        input_matrix.resize(size - 1); //Remove last character from end of a string

        for (size_t j = 0; j < size - 2; j++) {
            if (input_matrix.at(j) == ']') {
                if (input_matrix.at(j + 1) != ',' || input_matrix.at(j + 2) != ' ') {throw invalid_argument("invalid argument");}}
        }
        Rows = Split_by_character(input_matrix, ',');
        for (size_t i = 0; i < Rows.size(); i++) {
            Columns = Split_by_character(Rows[i], ' ');
            for (size_t j = 0; j < Columns.size(); j++) {
                //Parses str interpreting its content as a floating-point number, which is returned as a value of type double.
                ans.push_back(stod(Columns[j]));
            }}
        //Update the matrix with the data received at the input
        matrix.Rows = Rows.size();
        matrix.Columns = Columns.size();
        matrix.Values = ans;
        return input;
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
