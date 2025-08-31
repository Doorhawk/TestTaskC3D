#pragma once
#include <vector>
#include <stdexcept>


//-----------------------------
// Точка
//---
struct Point {
	double x, y;
};


//-----------------------------
// Ось
//---
enum class axis {
	x,
	y
};


//-----------------------------
// Сегмент сплайна
//---
struct SplineSegment {
	double t0, t1;   // границы сегмента
	double y0, y1;   // значения функции в узлах
	double m0, m1;   // вторые производные
};


//-----------------------------
// Сплайн Эрмита
//---
class Spline {
private:
	int n;
	std::vector<double> t;
	std::vector<SplineSegment> splineX, splineY;
	void buildSpline(const std::vector<double>& y, std::vector<SplineSegment>& out);
public:
	Spline();
	Spline(std::vector<Point>& points);
	void setBasePoints(const std::vector<Point>& points);
	double getValue(
					axis ax,	// ось x или y
					int der,    // 0 - значение, 1 - первая производная, 2 - вторая производная
					double tt   // пареметр от 0 до 1.0
					) const;
};