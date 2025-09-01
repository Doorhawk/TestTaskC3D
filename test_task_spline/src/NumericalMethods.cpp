#include "NumericalMethods.h"


//-----------------------------
// Точность
//---
double NumericalMethods::answerPrecision     = 1e-9;  // точность ответа
double NumericalMethods::derivativePrecision = 1e-9; // шаг для численой производной
double NumericalMethods::errorPrecision		 = 1e-9; // если при 1/x, x < errorPrecision решение будет прервано


//-----------------------------
// Решение линейной системы методом Гаусса
//---
Vec NumericalMethods::solveLinearSystem(const Mat& A, const Vec& b) {
	int n = A.size();
	Mat M = A;
	Vec x = b;

	for (int i = 0; i < n; i++) {
		int maxRow = i;
		for (int k = i + 1; k < n; k++) {
			if (std::abs(M[k][i]) > std::abs(M[maxRow][i])) {
				maxRow = k;
			}
		}
		std::swap(M[i], M[maxRow]);
		std::swap(x[i], x[maxRow]);

		if (std::abs(M[i][i]) < errorPrecision) {
			throw std::runtime_error("Division by 0");
		}

		double diag = M[i][i];
		for (int j = i; j < n; j++) M[i][j] /= diag;
		x[i] /= diag;

		for (int k = 0; k < n; k++) {
			if (k == i) continue;
			double factor = M[k][i];
			for (int j = i; j < n; j++) {
				M[k][j] -= factor * M[i][j];
			}
			x[k] -= factor * x[i];
		}
	}
	return x;
}


//-----------------------------
// Норма вектора
//---
double NumericalMethods::norm(const Vec& v) {
	double sum = 0.0;
	for (double vi : v)
		sum += vi * vi;
	return std::sqrt(sum);
}


//-----------------------------
// Решение системы методом Ньютона
//---
Vec NumericalMethods::newtonMethod(
	std::function<Vec(const Vec&)> f,			// Система уравнений
	std::function<Mat(const Vec&)> jacobian,	// Якобиан
	Vec x0,										// Начальные условия
	const Vec& minBounds,						// минимальные значения для каждой переменной
	const Vec& maxBounds,						// максимальные значения для каждой переменной
	int maxIter									// маскимум итераций
)
{
	Vec dx;
	Vec x = x0;
	for (int iter = 0; iter < maxIter; ++iter) {
		Mat J = jacobian(x);
		Vec F = f(x);

		dx = solveLinearSystem(J, F);

		for (size_t i = 0; i < x.size(); ++i) {
			x[i] -= dx[i];
			x[i] = std::max(minBounds[i], std::min(maxBounds[i], x[i]));
		}

		if (norm(dx) < answerPrecision)
			return x;
	}
	if (norm(dx) > answerPrecision) {
		std::ostringstream oss;
		oss << "The solution of the system by Newton's method did not converge to accuracy: "
			<< answerPrecision << " "
			<< "in " << maxIter << "steps\n";
		throw std::runtime_error(oss.str());
	}

	return x;
}


//-----------------------------
// Численный якобиан
//---
Mat NumericalMethods::numericalJacobian(
	std::function<Vec(const Vec&)> f,
	const Vec& x)
{

	double h = derivativePrecision;
	int n = x.size();
	Vec f0 = f(x);
	int m = f0.size();

	Mat J(m, Vec(n, 0.0));
	Vec xh = x;

	for (int j = 0; j < n; j++) {
		xh[j] += h;
		Vec f1 = f(xh);
		xh[j] -= h;

		for (int i = 0; i < m; i++) {
			J[i][j] = (f1[i] - f0[i]) / h;
		}
	}
	return J;
}


//-----------------------------
// Решение уравнения методом Ньютона
//---
double NumericalMethods::newtonMethod(
	const std::function<double(double)>& f,		// Функция
	const std::function<double(double)>& df,	// Производная
	double x0,									// Начальное приближение
	double minBound,
	double maxBound,
	size_t maxIter								// маскимум итераций
)
{
	double x = x0;
	double dx = 0;
	for (size_t i = 0; i < maxIter; ++i) {
		double fval = f(x);
		double dfval = df(x);
		if (std::abs(dfval) < errorPrecision) { // производная близка к нулю
			throw std::runtime_error("The derivative is close to zero in Newton's equation");
		}

		dx = fval / dfval;
		x -= dx;
		x = std::max(minBound, std::min(maxBound, x));

		if (std::abs(dx) < answerPrecision) {
			return x;
		}
	}
	if (std::abs(dx) > answerPrecision) {
		std::ostringstream oss;
		oss <<"The solution to the equation by Newton's method did not converge to accuracy: "
			<< answerPrecision << " "
			<< "in " << maxIter << "steps\n";
		throw std::runtime_error(oss.str());
	}
	return x;
	
}