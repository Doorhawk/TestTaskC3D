#pragma once
#include <vector>
#include <functional>
#include <stdexcept>
#include <sstream>

using Vec = std::vector<double>;
using Mat = std::vector<std::vector<double>>;

//-----------------------------
// Численные методы
//---
class NumericalMethods {
private:
	static double	answerPrecision;      // точность ответа
	static double	derivativePrecision;  // шаг для численой производной
	static double	errorPrecision;       // если при 1/x, x < errorPrecision решение будет прервано
public:
	static Vec		solveLinearSystem(const Mat& A, const Vec& b);
	static double	norm(const Vec& v);
	static Vec		newtonMethod(
						std::function<Vec(const Vec&)> f,			// Система уравнений
						std::function<Mat(const Vec&)> jacobian,	// Якобиан
						Vec x0,										// Начальные условия
						const Vec& minBounds,						// минимальные значения для каждой переменной
						const Vec& maxBounds,						// максимальные значения для каждой переменной
						int maxIter	= 30							// маскимум итераций
					);
	static double	newtonMethod(
						const std::function<double(double)>& f,		// Функция
						const std::function<double(double)>& df,	// Производная
						double x0,									// Начальное условие
						double minBound,
						double maxBound,
						size_t maxIter = 30							// маскимум итераций
					);		
	static Mat		numericalJacobian(
						std::function<Vec(const Vec&)> f,
						const Vec& x
					);
};
