#ifndef __MYMATRIX_H__
#define __MYMATRIX_H__

#include <iostream>
#include <algorithm>

using namespace std;

// 定义矩阵类
template<typename __T>
class TestMat
{
private:
	int row;
	int col;
	int type;
	__T **data;
public:
	TestMat(const int, const int, const int);
	TestMat(const TestMat<__T>&); // 拷贝构造函数
	~TestMat();
	// 运算符重载
	TestMat  operator+ (const TestMat<__T>&);
	TestMat  operator- (const TestMat<__T>&);
	TestMat  operator* (const TestMat<__T>&);
	TestMat& operator= (const TestMat<__T>&);

	TestMat  Transpose(); // 转置
	TestMat  num_multiply(__T); // 数乘
	TestMat  block_multi(const TestMat<__T>&, int); // 分块矩阵乘法
	TestMat  Hadamard(const TestMat<__T>&);
	TestMat<float>  Inverse(); // 求逆
	void  mat_show();
};


template<typename __T>
TestMat<__T>::TestMat(const int __row, const int __col, const int __type)
{
	row = __row;
	col = __col;
	type = __type;

	if (row > 0 && col > 0)
	{
		data = new __T*[row];
		for (int i = 0; i < row; ++i)
			data[i] = new __T[col];
		if (type == 0)
		{
			for (int i = 0; i < row; i++)
			{
				for (int j = 0; j < col; j++)
				{
					data[i][j] = (__T)0;
				}
			}
		}
		else if (type == 1)
		{
			int size = row < col ? row : col;
			for (int i = 0; i < size; i++)
			{
				data[i][i] = (__T)1;
			}
		}
		else if (type == 2)
		{
			for (int i = 0; i < row; i++)
			{
				for (int j = 0; j < col; j++)
				{
					data[i][j] = round((-1 + static_cast <__T> (rand()) / (static_cast <double> (RAND_MAX / 2))) * 1000) / 1000;
				}
			}
		}
		else
		{
			for (int i = 0; i < row; i++)
			{
				for (int j = 0; j < col; j++)
				{
					if (i == j)
						data[i][j] = j + 1;
					else
						data[i][j] = 1;
				}
			}
		}
	}
	else
	{
		row = 0;
		col = 0;
		type = -1;
		data = NULL;
	}
	return;
}

template<typename __T>
TestMat<__T>::TestMat(const TestMat<__T>& __temp)
{
	row = __temp.row;
	col = __temp.col;
	type = __temp.type;
	if (row > 0 && col > 0)
	{
		data = new __T*[row];
		for (int i = 0; i < row; ++i)
			data[i] = new __T[col];
		for (int i = 0; i < row; ++i)
			for (int j = 0; j < col; ++j)
				data[i][j] = __temp.data[i][j];
	}
	else
	{
		row = 0;
		col = 0;
		type = -1;
		data = NULL;
	}
	return;
}

template<typename __T>
TestMat<__T>::~TestMat()
{
	if (data)
	{
		for (int i = 0; i < row; ++i)
			delete[] data[i];
		delete[]data;
	}
	return;
}

template<typename __T>
TestMat<__T> TestMat<__T>::operator+(const TestMat<__T>& B)
{
	if (this->row == B.row&&this->col == B.col)
	{
		for (int i = 0; i < row; ++i)
			for (int j = 0; j < col; ++j)
				this->data[i][j] += B.data[i][j];
		return *this;
	}
	else
	{
		TestMat<__T> NullMatrix(0, 0);
		string Warning = "No matching matrix";
		throw Warning;
		return NullMatrix;
	}
}

template<typename __T>
TestMat<__T> TestMat<__T>::operator-(const TestMat<__T>& B)
{
	if (this->row == B.row&&this->col == B.col)
	{
		for (int i = 0; i < row; ++i)
			for (int j = 0; j < col; ++j)
				this->data[i][j] -= B.data[i][j];
		return *this;
	}
	else
	{
		TestMat<__T> NullMatrix(0, 0);
		string Warning = "No matching matrix";
		throw Warning;
		return NullMatrix;
	}
}

template<typename __T>
TestMat<__T>& TestMat<__T>::operator=(const TestMat<__T>& B)
{
	if (data)
	{
		for (int i = 0; i < row; ++i)
			delete[] data[i];
		delete data;
	}
	row = B.row;
	col = B.col;
	if (row > 0 && col > 0)
	{
		data = new __T*[row];
		for (int i = 0; i < row; ++i)
			data[i] = new __T[col];
		for (int i = 0; i < row; ++i)
			for (int j = 0; j < col; ++j)
				data[i][j] = B.data[i][j];
	}
	else
	{
		row = 0;
		col = 0;
		data = NULL;
	}
	return *this;
}

template<typename __T>
TestMat<__T> TestMat<__T>::operator*(const TestMat<__T>& B)
{
	TestMat<__T> NullMatrix(0, 0, 0);
	if (!this->row || !this->col || !B.row || !B.col)
	{
		string Warning = "No matching matrix";
		throw Warning;
	}
	else if (this->col != B.row)
	{
		string Warning = "No matching matrix";
		throw Warning;
	}
	else
	{
		TestMat<__T> Temp(this->row, B.col, 0);
		__T trans;
		for (int i = 0; i < Temp.row; ++i)
			for (int j = 0; j < Temp.col; ++j)
			{
				trans = 0;
				for (int k = 0; k < this->col; ++k)
					trans += this->data[i][k] * B.data[k][j];
				Temp.data[i][j] = trans;
			}
		return Temp;
	}
	return NullMatrix;
}

template<typename __T>
TestMat<__T> TestMat<__T>::Transpose()
{
	TestMat<__T> temp(this->col, this->row, 0);
	for (int i = 0; i < this->row; ++i)
		for (int j = 0; j < this->col; ++j)
			temp.data[j][i] = this->data[i][j];
	return temp;
}

template<typename __T>
TestMat<__T> TestMat<__T>::num_multiply(__T num)
{
	TestMat<__T> tempMat(this->col, this->row, 0);
	for (int i = 0; i < this->row; i++)
	{
		for (int j = 0; j < this->col; j++)
		{
			tempMat.data[i][j] = this->data[i][j] * num;
		}
	}
	return tempMat;
}

template<typename __T>
TestMat<__T> TestMat<__T>::block_multi(const TestMat<__T>& B, int block_size)
{
	TestMat<__T> NullMatrix(0, 0, 0);
	if (!this->row || !this->col || !B.row || !B.col)
	{
		string Warning = "No matching matrix";
		throw Warning;
	}
	else if (this->col != B.row)
	{
		string Warning = "No matching matrix";
		throw Warning;
	}
	else
	{
		TestMat<__T> Temp(this->row, B.col, 0);
		int chunk = 1;
#pragma omp for schedule (static, chunk)
		for (int bi = 0; bi < this->row; bi += block_size) {
			for (int bj = 0; bj < this->col; bj += block_size) {
				for (int bk = 0; bk < B.col; bk += block_size) {
					for (int i = bi; i < min(bi + block_size, this->row); ++i) {
						for (int j = bj; j < min(bj + block_size, B.col); ++j) {
							for (int k = bk; k < min(bk + block_size, B.col); ++k) {
								Temp.data[i][j] += this->data[i][k] * B.data[k][j];
							}
						}
					}
				}
			}
		}
		return Temp;
	}
	return NullMatrix;
}

template<typename __T>
TestMat<__T> TestMat<__T>::Hadamard(const TestMat<__T>& B)
{
	TestMat<__T> NullMatrix(0, 0, 0);
	if (!this->row || !this->col || !B.row || !B.col)
	{
		std::string Warning = "No matching matrix";
		throw Warning;
	}
	else if (this->row != B.row || this->col != B.col)
	{
		std::string Warning = "No matching matrix";
		throw Warning;
	}
	else
	{
		TestMat<__T> tempMat(this->row, this->col, 0);
		for (int i = 0; i < this->row; ++i)
			for (int j = 0; j < this->col; ++j)
				tempMat.data[i][j] = this->data[i][j] * B.data[i][j];
		return tempMat;
	}
	return NullMatrix;
}

template<typename __T>
TestMat<float> TestMat<__T>::Inverse()
{
	TestMat<float> tempMat(this->row, 2 * this->col, 0);
	for (int i = 0; i < this->row; i++)
	{
		for (int j = 0; j < this->col; j++)
		{
			tempMat.data[i][j] = this->data[i][j];
			if (i == j)
				tempMat.data[i][this->col + j] = (float)1;
		}
	}
	TestMat<float> resMat(this->row, this->col, 0);
	if (this->row != this->col)
	{
		string Warning = "No matching matrix";
		throw Warning;
	}
	__T temp = 0.0, r = 0.0;
	float minvalue = 1e-6;
	int flag;

	for (int j = 0; j < this->row; j++)
	{
		flag = j;
		for (int i = j + 1; i < this->row; i++)
			if (tempMat.data[i][j] > tempMat.data[flag][j])
				flag = i;

		if (fabs(tempMat.data[flag][j]) < minvalue) {
			cout << "Elements are too small to deal with !" << endl;
			exit(0);
		}
		if (flag != j) {
			for (int k = 0; k < 2 * this->row; k++)
			{
				temp = tempMat.data[j][k];
				tempMat.data[j][k] = tempMat.data[flag][k];
				tempMat.data[flag][k] = temp;
			}
		}
		for (int i = 0; i < this->row; i++)
		{
			if (i != j)
			{
				r = tempMat.data[i][j];
				for (int k = 0; k < 2 * this->row; k++)
					tempMat.data[i][k] -= tempMat.data[j][k] *
					r / tempMat.data[j][j];
			}
			else
			{
				r = tempMat.data[i][j];
				for (int k = 0; k < 2 * this->row; k++)
					tempMat.data[i][k] /= r;
			}
		}
	}
	for (int i = 0; i < this->row; i++)
	{
		for (int j = 0; j < this->col; j++)
		{
			resMat.data[i][j] = tempMat.data[i][this->col + j];
		}
	}
	return resMat;
}

template<typename __T>
void TestMat<__T>::mat_show()
{
	int matRow = this->row, matCol = this->col;
	if (matRow < 12 || matCol < 12)
	{
		int size = matRow <= matCol ? matRow : matCol;
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				cout << this->data[i][j] << "   ";
			}
			cout << endl;
		}
		cout << endl;
	}
	else
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << this->data[i][j] << "   ";
			}
			cout << endl;
		}
		cout << endl;
		int midRow = matRow / 2, midCol = matCol / 2;
		for (int i = midRow - 2; i < midRow + 2; i++)
		{
			for (int j = matCol - 2; j < matCol + 2; j++)
			{
				cout << this->data[i][j] << "   ";
			}
			cout << endl;
		}
		cout << endl;
		for (int i = matRow - 4; i < matRow; i++)
		{
			for (int j = matCol - 4; j < matCol; j++)
			{
				cout << this->data[i][j] << "   ";
			}
			cout << endl;
		}
		cout << endl;
	}
}

#endif