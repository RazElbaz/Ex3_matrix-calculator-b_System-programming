#include <map>
#include <vector>
#include <string>
#include <iostream>
using namespace std;

namespace zich{
    class Matrix{
    private:
        int Rows;
        int Columns;
        std::vector<double> Values;

    public:
        //****constructor****//
        Matrix(std::vector<double> MatrixVector, int Rows, int Columns) {
            // check if the matrix values are not negative
            if (Rows <= 0 || Columns <= 0) {
                throw invalid_argument("The matrix values cannot be negative");}
            if (MatrixVector.size() != (Rows * Columns)) {
                throw invalid_argument("The number of members of the matrix is not the same as the number of vector numbers");}
            this->Rows = Rows;
            this->Columns = Columns;
            this->Values=MatrixVector;
        }

        // Adding operators
        Matrix operator+(Matrix const &other);
        Matrix& operator+=(const Matrix &other);
        friend Matrix operator+(Matrix &other);

        // Subtraction operator
        friend Matrix operator-(const Matrix &matrix1, const Matrix &matrix2);
        Matrix& operator-=(const Matrix &other);
        friend Matrix operator-(Matrix &other);
        Matrix operator-(const Matrix &other);

        // Comparison operator
        bool operator>( Matrix &other);
        bool operator>=( Matrix &other);
        bool operator<( Matrix &other);
        bool operator<=( Matrix &other);
        bool operator==(const Matrix &other) const;
        bool operator!=(const Matrix &other) ;

        // Increasement operator
        //++i
        Matrix& operator++();
        //i++
        Matrix operator++(int);
        //i--
        Matrix operator--(int);
        //--i
        Matrix& operator--();

        // Scalar doubling
        friend Matrix operator*(const double scalar, Matrix &other);
        friend Matrix operator*(Matrix &other, const double scalar);
        Matrix& operator*=(const double scalar);

        // Multiplication of two matrices
        Matrix operator*(const Matrix &other);
        Matrix& operator*=(const Matrix &other);

        // friend global IO operators
        vector<string> Split_by_character(string will_split, char character);
        friend ostream& operator<< (ostream& output, const Matrix &matrix);
        friend istream& operator>> (istream& input, Matrix &matrix);
    };
}