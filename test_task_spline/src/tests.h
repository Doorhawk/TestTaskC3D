#pragma once
#include <iostream>
#include "Shapes.h"
#include "NumericalMethods.h"
#include "GeometricOperations.h"


//-----------------------------
// Сравнение по модулю
//---
bool eqval(double val1, double val2, double eps) {
	return abs(val1 - val2) < eps;
}


//-----------------------------
// Тестирование расстояния от точки до сплайна
//---
void testDistanceSplinePoint(Spline spline, Point targetPoints, Point correctPoint, double correctDist) {

	auto [gotPoint, gotDist] = GeometricOperations::distance(spline, targetPoints);

	std::cout << "calc " << gotPoint.x << ", " << gotPoint.y << " d = " << gotDist << '\n';
	std::cout << "need " << correctPoint.x << ", " << correctPoint.y << " d = " << correctDist << '\n';

	bool pass = eqval(gotPoint.x, correctPoint.x, 1e-9) &&
				eqval(gotPoint.y, correctPoint.y, 1e-9) &&
				eqval(gotDist, correctDist, 1e-9);

	std::cout << "\t\t\t\t" << (pass ? "good" : "bad") << '\n';
}


//-----------------------------
// Тестирование расстояния от сплафна до сплайна
//---
void testDistanceSplineSpline(Spline spline1, Spline spline2, Point correctPoint1, Point correctPoint2, double correctDist) {

	auto [gotPoint1, gotPoint2, gotDist] = GeometricOperations::distance(spline1, spline2, true);

	std::cout << "calc1 " << gotPoint1.x << ", " << gotPoint1.y << '\n';
	std::cout << "need1 " << correctPoint1.x << ", " << correctPoint1.y << '\n';

	std::cout << "calc2 " << gotPoint2.x << ", " << gotPoint2.y << '\n';
	std::cout << "need2 " << correctPoint2.x << ", " << correctPoint2.y << '\n';

	std::cout << "calc " << gotDist << '\n';
	std::cout << "need " << correctDist << '\n';

	bool pass = eqval(gotPoint1.x, correctPoint1.x, 1e-9) &&
				eqval(gotPoint1.y, correctPoint1.y, 1e-9) &&
				eqval(gotPoint2.x, correctPoint2.x, 1e-9) &&
				eqval(gotPoint2.y, correctPoint2.y, 1e-9) &&
				eqval(gotDist, correctDist, 1e-9);

	std::cout << "\t\t\t\t" << (pass ? "good" : "bad") << '\n';
}


//-----------------------------
// Тестирование пересечения сплайнов
//---
void testIntersectionSplineSpline(Spline spline1, Spline spline2, std::vector<Point>& correctPoints) {

	auto gotPoints = GeometricOperations::findIntersection(spline1, spline2);

	if (gotPoints.size() != correctPoints.size()) {
		std::cout << "size!=\n";
		std::cout << "\t\t\t\t" << "bad" << '\n';
		return;
	}

	for (int i = 0; i < gotPoints.size(); ++i) {
		std::cout << "calc1 " << gotPoints[i].x << ", " << gotPoints[i].y << '\n';
		std::cout << "need1 " << correctPoints[i].x << ", " << correctPoints[i].y << '\n';

		bool pass = eqval(gotPoints[i].x, correctPoints[i].x, 1e-9) &&
					eqval(gotPoints[i].y, correctPoints[i].y, 1e-9);

		if (!pass) {
			std::cout << "\t\t\t\t" << "bad" << '\n';
			return;
		}

	}
	std::cout << "\t\t\t\t" << "good" << '\n';
}


//-----------------------------
// Проверка функции дистанция от сплайна до точки
//---
void runTestDistanceSplinePoint() {

	// сплайн Прямая
	std::vector<Point> points = { {10, 0}, {0,0}, {-10, 0} };
	Spline spline(points); 
	
	// расстояние до одного края 
	testDistanceSplinePoint(spline, { 10,10 }, { 10,0 }, 10);
	// расстояние до центра
	testDistanceSplinePoint(spline, { 0,10 }, { 0,0 }, 10);
	// расстояние до дургого края 
	testDistanceSplinePoint(spline, { -10,-10 }, { -10,0 }, 10);
	// общий случай
	testDistanceSplinePoint(spline, { 5,5 }, { 5,0 }, 5);
	// расстояние до точки на сплайне
	testDistanceSplinePoint(spline, { 0,0 }, { 0,0 }, 0);

	// сплайн Кривая
	points = { {10, 0}, {0,10}, {-10, 0} };
	spline.setBasePoints(points);

	// расстояние до одного края 
	testDistanceSplinePoint(spline, { 20,0 }, { 10,0 }, 10);
	// расстояние до центра
	testDistanceSplinePoint(spline, { 0,20 }, { 0,10 }, 10);
	// расстояние до другого края
	testDistanceSplinePoint(spline, { -20,0 }, { -10,0 }, 10);
	// общий случай
	testDistanceSplinePoint(spline, { 10,10 }, { 5.5228735443959307,6.2669772664187056 }, 5.8292469521289467);

}


//-----------------------------
// Проверка функции дистанция от сплайна до сплайна
//---
void runTestDistanceSplineSpline() {

	// расстояние между центрами
	//	\ /
	//	/ \ 
	std::vector<Point> points = { {-11, 10}, {-1,0}, {-11, -10} };
	Spline spline1(points);
	points = { {11, 10}, {1,0}, {11, -10} };
	Spline spline2(points);
	testDistanceSplineSpline(spline1, spline2, { -1,0 }, { 1,0 }, 2);

	// расстояние между началами
	//	 / \
	//	/   \ 
	points = { {-1, 10}, {-6,0}, {-11, -10} };
	spline1.setBasePoints(points);
	points = { {1, 10}, {6,0}, {11, -10} };
	spline2.setBasePoints(points);
	testDistanceSplineSpline(spline1, spline2, { -1,10 }, { 1,10 }, 2);

	// расстояние между концами
	//	\   /
	//	 \ /
	points = { {-11, 10}, {-6,0}, {-1, -10} };
	spline1.setBasePoints(points);
	points = { {11, 10}, {6,0}, {1, -10} };
	spline2.setBasePoints(points);
	testDistanceSplineSpline(spline1, spline2, { -1,-10 }, { 1,-10 }, 2);

	// расстояние между концом и началом
	//	\
	//		\
	points = { {-11, 10}, {-6,0} };
	spline1.setBasePoints(points);
	points = { {6,0}, {11, -10}, };
	spline2.setBasePoints(points);
	testDistanceSplineSpline(spline1, spline2, { -6,0 }, { 6,0 }, 12);


	// расстояние между прямой и кривой
	//	\ |
	//	/ |
	points = { {-10, -10}, {0,0}, {-10, 10} };
	spline1.setBasePoints(points);
	points = { {1, 10}, {1,0}, {1, -10} };
	spline2.setBasePoints(points);
	testDistanceSplineSpline(spline1, spline2, { 0,0 }, { 1,0 }, 1);

	// расстояние между краем одного сплайна и серединой другого
	//	\	
	//	/ |
	points = { {-10, -10}, {0,0}, {-10, 10} };
	spline1.setBasePoints(points);
	points = { {1, 10}, {1,0} };
	spline2.setBasePoints(points);
	testDistanceSplineSpline(spline1, spline2, { 0,0 }, { 1,0 }, 1);

	// общий случай
	//	   /
	//	\  \
	//	/	
	points = { {-11, 15}, {-1,5}, {-11, -5} };
	spline1.setBasePoints(points);
	points = { {11, 10}, {1,0}, {11,-10} };
	spline2.setBasePoints(points);
	testDistanceSplineSpline(spline1, spline2, { -1.4519663133616083,3.2099427168089458 }, { 1.4519663133616083, 1.7900572831910535 }, 3.2324757300168350);
}


//-----------------------------
// Проверка функции пересечения сплайнов
//---
void runTestIntersectionSplineSpline() {

// Косание 
//	\ /
//	 |
//	/ \ 
	std::vector<Point> points = { {-10, 10}, {0,0}, {-10, -10} };
	Spline spline1(points);
	points = { {10, 10}, {0,0}, {10, -10} };
	Spline spline2(points);

	std::vector<Point> answer = { {0, 0} };
	testIntersectionSplineSpline(spline1, spline2, answer); // не находит касание

// 1-пересечение 2х пярмых
//	\ /
//	 \
//	/ \ 
	points = { {-10, 10}, {0,0}, {10, -10} };
	spline1.setBasePoints(points);
	points = { {10, 10}, {0,0}, {-10, -10} };
	spline2.setBasePoints(points);

	answer = { {0, 0} };
	testIntersectionSplineSpline(spline1, spline2, answer); 

// 2-пересечение 2х кривых
//	\ /
//	 \
//	/ \
//	\ /
//	 /
//	/ \ 
	points = { {-10, 10}, {5,0}, {-10, -10} };
	spline1.setBasePoints(points);
	points = { {10, 10}, {-5,0}, {10, -10} };
	spline2.setBasePoints(points);

	answer = { {0, 5.1829802948862937},{0,-5.1829802948862937} };
	testIntersectionSplineSpline(spline1, spline2, answer);
}