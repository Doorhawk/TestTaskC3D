#pragma once
#include <vector>
#include "Shapes.h"
#include "NumericalMethods.h"


//-----------------------------
// Геометрические операции
//---
class GeometricOperations {
private:
	static double precision;
public:
	static std::tuple<Point, Point, double> distance(
														const Spline& s1, 
														const Spline& s2, 
														bool analit	= true		// использовать аналитический якобиан или численный
													);
	static std::tuple<Point, double>		distance(const Spline& spline, const Point& point);
	static std::vector<Point>				findIntersection(const Spline& s1, const Spline& s2, bool analit = true);
	static std::vector<Point>				findIntersection(Point p1, Point p2, Point q1, Point q2);
};

