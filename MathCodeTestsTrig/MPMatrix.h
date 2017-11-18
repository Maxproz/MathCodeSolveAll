#pragma once


//#ifndef MPMatrix
//#define MPMatrix

#include <vector>

class MPMatrix
{
public:

	//MPMatrix(size_t rows, size_t cols);
	MPMatrix() = default;
	MPMatrix(size_t columns, size_t rows) :
		mCols(columns), mRows(rows),
		Data(mCols, std::vector<double>(mRows))
	{}

	double& operator()(size_t i, size_t j);
	double operator()(size_t i, size_t j) const;

	MPMatrix(const MPMatrix&) = default; // copy constructor
	MPMatrix(MPMatrix&&) = default; // move constructor
	MPMatrix& operator=(const MPMatrix&) = default; // copy assign
	MPMatrix& operator=(MPMatrix&&) = default; // move assign
	virtual ~MPMatrix() = default;
	
	//MPMatrix& operator =(MPMatrix& o);
	MPMatrix& operator +(MPMatrix& o);
	MPMatrix& operator -(MPMatrix& o);
	MPMatrix& operator *(MPMatrix& o);
	//MPMatrix& operator ^(int power);

	friend bool operator== (const MPMatrix &c1, const MPMatrix &c2);
	friend bool operator!= (const MPMatrix &c1, const MPMatrix &c2);

	MPMatrix& MultiplyByScaler(const double& ScalerValue);

	void print();

	friend void printGuass(const MPMatrix& InMatrix);
	friend std::vector<double> gauss(const MPMatrix& InMatrix);

	const size_t GetRows() const { return mRows; }

	friend const bool AreInverses(const MPMatrix& first, const MPMatrix& second);

	bool HasInverse2by2Matrix();

	MPMatrix GetInverseMatrix2By2();

private:
	size_t mRows;
	size_t mCols;
	std::vector<std::vector<double>> Data;
	
};

void printGuass(const MPMatrix& InMatrix);

std::vector<double> gauss(const MPMatrix& InMatrix);

void Create3By3MatrixAndUseGuassMethod();

std::vector<std::vector<double> > make_idty_matrix(int n);

const bool AreInverses(const MPMatrix& first, const MPMatrix& second);

bool operator== (const MPMatrix &c1, const MPMatrix &c2);
bool operator!= (const MPMatrix &c1, const MPMatrix &c2);


// PROPERTIES OF DETERMINANTS
// #1 If the matrix is in upper triangular form,
// the determinant equals the product of entries down the main diagonal.

// #2 When two rows are interchanged, the determinant changes sign.

// #3 If either two rows or two columns are identical, the determinant equals zero.

// #4 If a matrix contains either a row of zeros or a column of zeros, 
// the determinant equals zero.

// #5 The determinant of an inverse matrix A^-1
// is the reciprocal of the determinant of the matrix A.

// #6 If any row or column is multiplied by a constant, 
// the determinant is multiplied by the same factor.


double GetDeterminent2By2Matrix(const MPMatrix &Matrix);
//double GetDeterminent3By3Matrix(const MPMatrix &Matrix);

void Solve2By2MatrixCramersRule(
	const double& a1,
	const double& b1,
	const double& c1,
	const double& a2,
	const double& b2,
	const double& c2);

double Solve2By2MatrixForYCramersRule(
	const double& a1,
	const double& b1,
	const double& c1,
	const double& a2,
	const double& b2,
	const double& c2);

double GetDeterminent3By3Matrix(
	const double& a1,
	const double& b1,
	const double& c1,
	const double& a2,
	const double& b2,
	const double& c2,
	const double& a3,
	const double& b3,
	const double& c3);


void Solve3By3MatrixCramersRule(
	const double& a1,
	const double& b1,
	const double& c1,
	const double& d1,
	const double& a2,
	const double& b2,
	const double& c2,
	const double& d2,
	const double& a3,
	const double& b3,
	const double& c3,
	const double& d3);

//#endif // !MPMatrix.h