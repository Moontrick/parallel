#include <iomanip>
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <time.h>
#include <omp.h>
using namespace std;
#define PI (3.14159265358979323846)
//Function for simple initialization of input signal elements
void DummyDataInitialization(complex<double>* mas, int size) {
	for (int i = 0; i < size; i++) {
		mas[i] = (PI / 2.0) - ((PI * (i + 1)) / 1024.0);
	}
}
//Function for memory allocation and data initialization
void ProcessInitialization(complex<double>*& inputSignal,
	complex<double>*& outputSignal, int& size) {
	cout << "Input signal length = " << size << endl;
	inputSignal = new complex<double>[size];
	outputSignal = new complex<double>[size];
	DummyDataInitialization(inputSignal, size);
	//Initialization of input signal elements - tests
	//RandomDataInitialization(inputSignal, size);
	//Computational experiments
	//RandomDataInitialization(inputSignal, size);
}
//Function for computational process temination
void ProcessTermination(complex<double>*& inputSignal,
	complex<double>*& outputSignal) {
	delete[] inputSignal;
	inputSignal = NULL;
	delete[] outputSignal;
	outputSignal = NULL;
}
void BitReversing(complex<double>* inputSignal,
	complex<double>* outputSignal, int size) {
	int j = 0, i = 0;
//#pragma omp parallel for
	for (i = 0; i < size; ++i) {
		if (j > i) {
			outputSignal[i] = inputSignal[j];
			outputSignal[j] = inputSignal[i];
		}
		else {
			if (j == i) {
				outputSignal[i] = inputSignal[i];
			}
		}
		int m = size >> 1;

		while ((m >= 1) && (j >= m))
		{
			j -= m;
			m = m >> 1;
		}
		j += m;
	}
}
__inline void Butterfly(complex<double>* signal,
	complex<double> u, int offset, int butterflySize) {
	complex<double> tem = signal[offset + butterflySize] * u;
	signal[offset + butterflySize] = signal[offset] - tem;
	signal[offset] += tem;
}
void ParallelFFTCalculation(complex<double>* signal, int size) {
	int m = 0;
	for (int tmp_size = size; tmp_size > 1; tmp_size /= 2, m++);
	for (int p = 0; p < m; p++)
	{
		int butterflyOffset = 1 << (p + 1);
		int butterflySize = butterflyOffset >> 1;
		double coeff = PI / butterflySize;
//#pragma omp parallel for
		for (int i = 0; i < size / butterflyOffset; i++)
			for (int j = 0; j < butterflySize; j++)
				Butterfly(signal, complex<double>(cos(-j * coeff),
					sin(-j * coeff)), j + i * butterflyOffset, butterflySize);
	}
}
// FFT computation
void ParallelFFT(complex<double>* inputSignal,
	complex<double>* outputSignal, int size) {
	BitReversing(inputSignal, outputSignal, size);
	ParallelFFTCalculation(outputSignal, size);
}
void PrintSignal(complex<double>* signal, int size) {
	cout << "Result signal" << endl;
	for (int i = 0; i < size; i++)
		cout << signal[i] << endl;
}

void TestResult(complex<double>* inputSignal, complex<double>* outputSignal, int size)
{
	complex<double>* testFunc = new complex<double>[size];
	for (int i = 0; i < size; i++) {
		testFunc[i] = outputSignal[0].real() / 4;
		for (int k = 1; k < size; k++) {
			testFunc[i] += (outputSignal[k].real() / 2 * cos((k * 2 * PI
				 * (i + 1))/1024.0) - outputSignal[k].imag() / 2 * sin((k * 2 * PI * (i + 1)) / 1024.0));
		}
		testFunc[i] /= size / 2;
	}
	vector<double> testFunctoch(size);
	for (int i = 0; i < size; i++) {
		testFunctoch[i] = 0;
		for (int k = 1; k < size * 16; k++)
			testFunctoch[i] += sin((k * 2 * PI * (i + 1))/1024.0) / k;
	}


	for (int i = 0; i < size; i ++) {
		cout  << testFunc[i].real() << setw(10) << testFunctoch[i] << setw(10) << inputSignal[i].real() << endl;
	}
	
}

int main()
{
	setlocale(LC_ALL, "Rus");
	complex<double>* inputSignal = NULL;

	complex<double>* outputSignal = NULL;
	int size = 1024;
	double startTime;
	double duration;
	double minDuration = DBL_MAX;
	cout << "Fast Fourier Transform" << endl;
	// Memory allocation and data initialization
	ProcessInitialization(inputSignal, outputSignal, size);
	ParallelFFT(inputSignal, outputSignal, size);
	int half_size = size / 2;
	for (int i = 0; i < size; i++) {
		if (abs(outputSignal[i].real()) > 0.000001 || abs(outputSignal[i].imag()) > 0.000001) {
			outputSignal[i] = sqrt(pow(outputSignal[i].real(), 2) + pow(outputSignal[i].imag(), 2)) / half_size;
			//cout << outputSignal[i].real() << endl;
		}
		else
			outputSignal[i] = 0;
	}
	ProcessTermination(inputSignal, outputSignal);

	size = 1024;
	ProcessInitialization(inputSignal, outputSignal, size);
	ParallelFFT(inputSignal, outputSignal, size);
	// Result signal output
	//PrintSignal(outputSignal, size);
	// Computational process termination
	TestResult(inputSignal, outputSignal, size);
	
	ProcessTermination(inputSignal, outputSignal);
	return 0;
}
