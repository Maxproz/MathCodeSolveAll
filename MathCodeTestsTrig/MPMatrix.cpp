#include "MPMatrix.h"

#include <iostream>

double& MPMatrix::operator()(size_t i, size_t j)
{
	return Data[i][j];
}

double MPMatrix::operator()(size_t i, size_t j) const
{
	return Data[i][j];
}

/*
MPMatrix& MPMatrix::operator =(MPMatrix& o)
{
	Data.resize(o.Data.size());

	for (int i = 0; i < Data.size(); i++)
	{
		Data[i].resize(o.Data[i].size());
	}
	for (int i = 0; i < Data.size(); i++)
	{
		for (int j = 0; j<Data[i].size(); j++)
		{
			Data[i][j] = o.Data[i][j];
		}
	}
	return *this;
}*/


MPMatrix& MPMatrix::operator +(MPMatrix& o)
{
	for (int i = 0; i < Data.size(); i++)
	{
		for (int j = 0; j<Data[i].size(); j++)
		{
			Data[i][j] = (Data[i][j] + o.Data[i][j]);
		}
	}

	return *this;
}

MPMatrix& MPMatrix::operator -(MPMatrix& o)
{
	for (int i = 0; i < Data.size(); i++)
	{
		for (int j = 0; j < Data[i].size(); j++)
		{
			Data[i][j] = Data[i][j] - o.Data[i][j];
		}
	}

	return *this;
}


MPMatrix& MPMatrix::operator *(MPMatrix& o)
{
	if (Data[0].size() != o.Data.size()) return *this;

	MPMatrix tm;
	tm.Data.resize(Data.size());
	for (int i = 0; i < tm.Data.size(); i++)
	{
		tm.Data[i].resize(o.Data[0].size());
	}

	for (int i = 0; i < tm.Data.size(); i++)
	{
		for (int j = 0; j < tm.Data[i].size(); j++)
		{
			tm.Data[i][j] = 0;

			for (int c = 0; c < Data[i].size(); c++)
			{
				tm.Data[i][j] += Data[i][c] * o.Data[c][j];
			}

		}
	}

	*this = tm;
	return *this;
}

/*
MPMatrix& MPMatrix::operator ^(int power)
{
	MPMatrix  tM2;
	tM2 = *this;

	//   not <= below \/ because first time counts as 2
	for (int i = 1; i<power; ++i)
		*this = (*this) * (tM2);

	return *this;
}
*/

MPMatrix & MPMatrix::MultiplyByScaler(const double & ScalerValue)
{
	for (int i = 0; i < Data.size(); i++)
	{
		for (int j = 0; j < Data[i].size(); j++)
		{
			Data[i][j] = (Data[i][j] * ScalerValue);
		}
	}

	return *this;

}

void MPMatrix::print()
{
	for (int i = 0; i < Data.size(); i++)
	{
		for (int j = 0; j < Data[i].size(); j++)
		{
			std::cout << Data[i][j] << ' ';
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

bool MPMatrix::HasInverse2by2Matrix()
{

	double a, b, c, d;
	a = this->operator()(0, 0);
	b = this->operator()(0, 1);
	c = this->operator()(1, 0);
	d = this->operator()(1, 1);

	if ((a * d) - (b * c) == 0)
	{
		return false;
	}

	return true;
}

MPMatrix MPMatrix::GetInverseMatrix2By2()
{
	MPMatrix TempMatrix = *this;

	if (!(HasInverse2by2Matrix()))
	{
		throw std::exception("No inverse");
		//return TempMatrix;
	}
	
	double a, b, c, d;
	a = TempMatrix(0, 0);
	b = TempMatrix(0, 1);
	c = TempMatrix(1, 0);
	d = TempMatrix(1, 1);

	double scalar = 1.0 / ((a * d) - (b * c));

	std::swap(TempMatrix(0, 0), TempMatrix(1, 1));
	TempMatrix(0, 1) = TempMatrix(0, 1) * -1;
	TempMatrix(1, 0) = TempMatrix(1, 0) * -1;

	TempMatrix.MultiplyByScaler(scalar);

	return TempMatrix;
}



void printGuass(const MPMatrix& InMatrix)
{
	std::vector<std::vector<double>> A = InMatrix.Data;

	int n = A.size();
	for (int i = 0; i<n; i++)
	{
		for (int j = 0; j<n + 1; j++)
		{
			std::cout << A[i][j] << "\t";
			if (j == n - 1)
			{
				std::cout << "| ";
			}
		}
		std::cout << "\n";
	}
	std::cout << std::endl;
}


std::vector<double> gauss(const MPMatrix& InMatrix)
{
	std::vector< std::vector<double> > A = InMatrix.Data;

	int n = A.size();

	for (int i = 0; i<n; i++)
	{
		// Search for maximum in this column
		double maxEl = abs(A[i][i]);
		int maxRow = i;
		for (int k = i + 1; k<n; k++)
		{
			if (abs(A[k][i]) > maxEl)
			{
				maxEl = abs(A[k][i]);
				maxRow = k;
			}
		}

		// Swap maximum row with current row (column by column)
		for (int k = i; k<n + 1; k++)
		{
			double tmp = A[maxRow][k];
			A[maxRow][k] = A[i][k];
			A[i][k] = tmp;
		}

		// Make all rows below this one 0 in current column
		for (int k = i + 1; k<n; k++)
		{
			double c = -A[k][i] / A[i][i];
			for (int j = i; j<n + 1; j++)
			{
				if (i == j)
				{
					A[k][j] = 0;
				}
				else
				{
					A[k][j] += c * A[i][j];
				}
			}
		}
	}

	// Solve equation Ax=b for an upper triangular matrix A
	std::vector<double> x(n);
	for (int i = n - 1; i >= 0; i--)
	{
		x[i] = A[i][n] / A[i][i];
		for (int k = i - 1; k >= 0; k--)
		{
			A[k][n] -= A[k][i] * x[i];
		}
	}
	return x;
}

void Create3By3MatrixAndUseGuassMethod()
{
	int32_t MatrixTestSize = 3;
	uint8_t ExtraAugmentedRowSlow = MatrixTestSize + 1;

	MPMatrix MatrixTestGuass(MatrixTestSize, ExtraAugmentedRowSlow);

	// Read input data
	std::cout << "Enter Input Data for the Matrix\n";
	std::cout << "Enter One Number at a time and press enter\n";
	std::cout << "The first 3 numbers are row 1, the second 3 row 2 etc...\n";

	for (int i = 0; i < MatrixTestSize; i++)
	{
		for (int j = 0; j < MatrixTestSize; j++)
		{
			std::cin >> MatrixTestGuass(i, j);
		}
	}

	std::cout << "Enter Augmented row Data for rows 1-3\n";
	std::cout << "Enter One number at a time and press enter.\n";

	for (int i = 0; i < MatrixTestSize; i++)
	{
		std::cin >> MatrixTestGuass(i, MatrixTestSize);
	}

	std::cout << "Preparing Result.. Please Wait\n";

	// Print input
	printGuass(MatrixTestGuass);

	// Calculate solution

	std::vector<double> x(MatrixTestSize);
	x = gauss(MatrixTestGuass);

	// Print result
	std::cout << "Result:\t (";

	for (int i = 0; i < MatrixTestSize; i++)
	{
		std::cout << x[i];

		const uint8_t CommaSizeCheck = i + 1;
		if ((CommaSizeCheck) < MatrixTestSize)
		{
			std::cout << ",";
		}
	}
	std::cout << ")\n";
	std::cout << std::endl;

}


std::vector<std::vector<double> > make_idty_matrix(int n)
{
	std::vector<std::vector<double> > idty(n, std::vector<double>(n, 0));
	for (int i = 0; i < n; ++i)
	{
		idty[i][i] = 1;
	}
	
	return idty;
}


const bool AreInverses(const MPMatrix & first, const MPMatrix & second)
{
	size_t Rows = first.GetRows();
	
	MPMatrix temp;
	temp.Data = make_idty_matrix(Rows);


	MPMatrix templeft = first;
	MPMatrix tempright = second;

	MPMatrix temp2 = templeft * tempright;

	if (temp == temp2)
	{
		return true;
	}

	return false;
}

bool operator==(const MPMatrix & c1, const MPMatrix & c2)
{
	for (int i = 0; i < c1.Data.size(); i++)
	{
		for (int j = 0; j < c1.Data[i].size(); j++)
		{
			if (c1.Data[i][j] != c2.Data[i][j])
			{
				return false;
			}
		}
	}

	return true;
}

bool operator!=(const MPMatrix & c1, const MPMatrix & c2)
{
	return !(operator==(c1,c2));
}

double GetDeterminent2By2Matrix(const MPMatrix& Matrix)
{
	MPMatrix TempMatrix = Matrix;

	double a = TempMatrix(0, 0);
	double b = TempMatrix(0, 1);
	double c = TempMatrix(1, 0);
	double d = TempMatrix(1, 1);

	double ad = a*d;
	double cb = c*b;

	return ad - cb;
}

void Solve2By2MatrixCramersRule(const double & a1, const double & b1, const double & c1, const double & a2, const double & b2, const double & c2)
{
	MPMatrix Numerator(2,2);
	MPMatrix Denominator(2, 2);

	Numerator(0, 0) = c1;
	Numerator(0, 1) = b1;
	Numerator(1, 0) = c2;
	Numerator(1, 1) = b2;

	Denominator(0, 0) = a1;
	Denominator(0, 1) = b1;
	Denominator(1, 0) = a2;
	Denominator(1, 1) = b2;

	double NumeratorDet = GetDeterminent2By2Matrix(Numerator);
	double DenominatorDet = GetDeterminent2By2Matrix(Denominator);

	if (DenominatorDet == 0.0)
	{
		throw std::domain_error("Denomator Cannot Equal 0");
	}

	double XResult = NumeratorDet / DenominatorDet;
	double YResult = Solve2By2MatrixForYCramersRule(a1, b1, c1, a2, b2, c2);


	std::cout << "The solution is: (" << XResult << ", " << YResult << ")" << std::endl;

}

double Solve2By2MatrixForYCramersRule(const double & a1, const double & b1, const double & c1, const double & a2, const double & b2, const double & c2)
{
	MPMatrix Numerator(2, 2);
	MPMatrix Denominator(2, 2);

	Numerator(0, 0) = a1;
	Numerator(0, 1) = c1;
	Numerator(1, 0) = a2;
	Numerator(1, 1) = c2;

	Denominator(0, 0) = a1;
	Denominator(0, 1) = b1;
	Denominator(1, 0) = a2;
	Denominator(1, 1) = b2;

	double NumeratorDet = GetDeterminent2By2Matrix(Numerator);
	double DenominatorDet = GetDeterminent2By2Matrix(Denominator);

	if (DenominatorDet == 0.0)
	{
		throw std::domain_error("Denomator Cannot Equal 0");
	}

	return NumeratorDet / DenominatorDet;
}

double GetDeterminent3By3Matrix(const double & a1, const double & b1, const double & c1, const double & a2, const double & b2, const double & c2, const double & a3, const double & b3, const double & c3)
{
	double a1b2c3 = a1*b2*c3;
	double b1c2a3 = b1*c2*a3;
	double c1a2b3 = c1*a2*b3;

	double a3b2c1 = a3*b2*c1;
	double b3c2a1 = b3*c2*a1;
	double c3a2b1 = c3*a2*b1;

	double Determinent{ 0.0 };
	Determinent = a1b2c3 + b1c2a3 + c1a2b3 - (a3b2c1)-(b3c2a1)-(c3a2b1);
	
	return Determinent;
}

void Solve3By3MatrixCramersRule(const double & a1, const double & b1, const double & c1, const double & d1, const double & a2, const double & b2, const double & c2, const double & d2, const double & a3, const double & b3, const double & c3, const double & d3)
{
	MPMatrix Denominator(3, 3);
	Denominator(0, 0) = a1;
	Denominator(0, 1) = b1;
	Denominator(0, 2) = c1;
	Denominator(1, 0) = a2;
	Denominator(1, 1) = b2;
	Denominator(1, 2) = c2;
	Denominator(2, 0) = a3;
	Denominator(2, 1) = b3;
	Denominator(2, 2) = c3;

	double DenominatorDet = GetDeterminent3By3Matrix(a1,b1,c1,a2,b2,c2,a3,b3,c3);

	MPMatrix SolveXNumerator(3, 3);
	SolveXNumerator(0, 0) = d1;
	SolveXNumerator(0, 1) = b1;
	SolveXNumerator(0, 2) = c1;
	SolveXNumerator(1, 0) = d2;
	SolveXNumerator(1, 1) = b2;
	SolveXNumerator(1, 2) = c2;
	SolveXNumerator(2, 0) = d3;
	SolveXNumerator(2, 1) = b3;
	SolveXNumerator(2, 2) = c3;

	double XNumeratorDet = GetDeterminent3By3Matrix(d1, b1, c1, d2, b2, c2, d3, b3, c3);

	MPMatrix SolveYNumerator(3, 3);
	SolveYNumerator(0, 0) = a1;
	SolveYNumerator(0, 1) = d1;
	SolveYNumerator(0, 2) = c1;
	SolveYNumerator(1, 0) = a2;
	SolveYNumerator(1, 1) = d2;
	SolveYNumerator(1, 2) = c2;
	SolveYNumerator(2, 0) = a3;
	SolveYNumerator(2, 1) = d3;
	SolveYNumerator(2, 2) = c3;

	double YNumeratorDet = GetDeterminent3By3Matrix(a1, d1, c1, a2, d2, c2, a3, d3, c3);

	MPMatrix SolveZNumerator(3, 3);
	SolveZNumerator(0, 0) = a1;
	SolveZNumerator(0, 1) = b1;
	SolveZNumerator(0, 2) = d1;
	SolveZNumerator(1, 0) = a2;
	SolveZNumerator(1, 1) = b2;
	SolveZNumerator(1, 2) = d2;
	SolveZNumerator(2, 0) = a3;
	SolveZNumerator(2, 1) = b3;
	SolveZNumerator(2, 2) = d3;

	double ZNumeratorDet = GetDeterminent3By3Matrix(a1, b1, d1, a2, b2, d2, a3, b3, d3);

	double XResult = XNumeratorDet / DenominatorDet;
	double YResult = YNumeratorDet / DenominatorDet;
	double ZResult = ZNumeratorDet / DenominatorDet;

	std::cout << "Solution: (" << XResult << ", " << YResult << ", " << ZResult << ")\n";

}


//double GetDeterminent3By3Matrix(const MPMatrix & Matrix)
//{
//
//
//	return 0.0;
//}

