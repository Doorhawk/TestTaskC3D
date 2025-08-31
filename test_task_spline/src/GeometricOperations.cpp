#include "GeometricOperations.h"

double GeometricOperations::precision = 1e-9;

//-----------------------------
// Пересечение отрезков
//---
std::vector<Point> GeometricOperations::findIntersection(Point p1, Point p2, Point q1, Point q2) {

	double A1 = p2.x - p1.x;
	double B1 = -(q2.x - q1.x);
	double C1 = q1.x - p1.x;

	double A2 = p2.y - p1.y;
	double B2 = -(q2.y - q1.y);
	double C2 = q1.y - p1.y;

	double det = A1 * B2 - A2 * B1;

	if (std::abs(det) < precision) {
		return {};
	}

	double t = (C1 * B2 - C2 * B1) / det;
	double u = (A1 * C2 - A2 * C1) / det;

	if (t < -precision || t > 1.0 + precision) return {};
	if (u < -precision || u > 1.0 + precision) return {};

	Point intersection;
	intersection.x = p1.x + t * (p2.x - p1.x);
	intersection.y = p1.y + t * (p2.y - p1.y);

	return { intersection };
}


//-----------------------------
// Расстояние между сплайнами
//---
std::tuple<Point, Point, double> GeometricOperations::distance(const Spline& s1, const Spline& s2, bool analit) {

	// расстояние
	// ro(t,u)^2 = D(t,u) = ||P1 - P2||^2 = (x1(t) - x2(u))^2 + (y1(t) - y2(u))^2 
	// минимизировать сисиму 
	// {D't = 0, D'u = 0} 

	// основная функция
	auto F = [&](const Vec& v) {

		// F1 = dD/dt = (x1(t)-x2(u))x1'(t) + (y1(t)-y2(u))y1'(t) = 0
		// F2 = dD/du = (x1(t)-x2(u))x2'(u) + (y1(t)-y2(u))y2'(u) = 0

		double t = v[0];
		double u = v[1];

		double x1 = s1.getValue(axis::x, 0, t);
		double y1 = s1.getValue(axis::y, 0, t);
		double x2 = s2.getValue(axis::x, 0, u);
		double y2 = s2.getValue(axis::y, 0, u);

		double dx = x1 - x2;
		double dy = y1 - y2;

		double x1p = s1.getValue(axis::x, 1, t);
		double y1p = s1.getValue(axis::y, 1, t);
		double x2p = s2.getValue(axis::x, 1, u);
		double y2p = s2.getValue(axis::y, 1, u);

		Vec res(2);
		res[0] = dx * x1p + dy * y1p;   // F1
		res[1] = dx * x2p + dy * y2p;   // F2
		return res;
		};
	// Якобиан
	std::function<Mat(const Vec&)> J;
	if (analit)
		J = [&](const Vec& v) -> Mat {
		double t = v[0];
		double u = v[1];

		double x1 = s1.getValue(axis::x, 0, t);
		double y1 = s1.getValue(axis::y, 0, t);
		double x2 = s2.getValue(axis::x, 0, u);
		double y2 = s2.getValue(axis::y, 0, u);

		double dx = x1 - x2;
		double dy = y1 - y2;

		double x1p = s1.getValue(axis::x, 1, t);
		double y1p = s1.getValue(axis::y, 1, t);
		double x2p = s2.getValue(axis::x, 1, u);
		double y2p = s2.getValue(axis::y, 1, u);

		double x1pp = s1.getValue(axis::x, 2, t);
		double y1pp = s1.getValue(axis::y, 2, t);
		double x2pp = s2.getValue(axis::x, 2, u);
		double y2pp = s2.getValue(axis::y, 2, u);

		Mat J(2, Vec(2));

		// dF1/dt
		J[0][0] = x1p * x1p + y1p * y1p + dx * x1pp + dy * y1pp;
		// dF1/du
		J[0][1] = -(x1p * x2p + y1p * y2p);
		// dF2/dt
		J[1][0] = x1p * x2p + y1p * y2p;
		// dF2/du
		J[1][1] = -(x2p * x2p + y2p * y2p + (x2 - x1) * x2pp + (y2 - y1) * y2pp);

		return J;
		};
	else
		J = [&](const Vec& v) -> Mat {
		return NumericalMethods::numericalJacobian(F, v);
		};

	// --- шаг 1: грубый перебор ---
	const int nGrid = 25;
	double minDist2 = std::numeric_limits<double>::max();
	double bestT = 0, bestU = 0;

	for (int i = 0; i < nGrid; ++i) {
		double t = double(i) / (nGrid - 1);
		double x1 = s1.getValue(axis::x, 0, t);
		double y1 = s1.getValue(axis::y, 0, t);

		for (int j = 0; j < nGrid; ++j) {
			double u = double(j) / (nGrid - 1);
			double x2 = s2.getValue(axis::x, 0, u);
			double y2 = s2.getValue(axis::y, 0, u);

			double dist2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
			if (dist2 < minDist2) {
				minDist2 = dist2;
				bestT = t;
				bestU = u;
			}
		}
	}


	// --- шаг 2: уточнение методом Ньютона ---
	Vec start = { bestT, bestU };
	Vec minBounds = { 0.0, 0.0 };
	Vec maxBounds = { 1.0, 1.0 };

	double tFinal = bestT;
	double uFinal = bestU;

	if (tFinal != 1.0 && tFinal != 0.0 && uFinal != 1.0 && uFinal != 0.0) {

		Vec sol = NumericalMethods::newtonMethod(F, J, start, minBounds, maxBounds);

		tFinal = sol[0];
		uFinal = sol[1];

	}
	

	double dist = sqrt(minDist2);
	Point p1 = { s1.getValue(axis::x,0,tFinal), s1.getValue(axis::y,0,tFinal) };
	Point p2 = { s2.getValue(axis::x,0,uFinal), s2.getValue(axis::y,0,uFinal) };
	dist = sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));




	// --- шаг 3: проверка краевых точек ---
	// 1. крайние точки первого сплайна к ближайшей точке второго сплайна
	for (double tEdge : { 0.0, 1.0 }) {
		Point edgeP1{ s1.getValue(axis::x,0, tEdge), s1.getValue(axis::y,0, tEdge) };
		auto [closestOnS2,distS2] = distance(s2,edgeP1);

		if (distS2 < dist) {
			p1 = edgeP1;
			p2 = closestOnS2;
			dist = distS2;
		}
	}

	// 2. крайние точки второго сплайна к ближайшей точке первого сплайна
	for (double uEdge : { 0.0, 1.0 }) {
		Point edgeP2{ s2.getValue(axis::x,0, uEdge), s2.getValue(axis::y,0, uEdge) };
		auto [closestOnS1,distS1] = distance(s1, edgeP2);
		if (distS1 < dist) {
			p1 = closestOnS1;
			p2 = edgeP2;
			dist = distS1;
		}
	}


	return { p1, p2, dist };

}


//-----------------------------
// Расстояние между точкой и сплайном
//---
std::tuple<Point, double> GeometricOperations::distance(const Spline& spline, const Point& point) {

	// ищем точку ближайшую к P0(x0,y0)
	// нужно минимизировать функцию расстояния
	// R(t) = (x(t)-x0)^2 + (y(t)-y0)^2
	// R' = 2(x-x0)x' + 2(y-y0)y'
	// R'' = 2((x')^2 + (x-x0)x'' + (y')^2 + (y-y0)y'')
	// минимум функции растояние это решенеие уравнения R' = 0
	// методом ньютона получаем x = x0 - f/df -> f = R' -> x = x0 - R'/R''

	// --- шаг 1: грубый перебор ---
	const int nGrid = 25; // сетка для грубого поиска
	double minDist2 = std::numeric_limits<double>::max();
	double bestT = 0;


	for (int i = 0; i < nGrid; ++i) {
		double t = double(i) / (nGrid - 1);

		double x = spline.getValue(axis::x, 0, t);
		double y = spline.getValue(axis::y, 0, t);

		double dx = x - point.x;
		double dy = y - point.y;
		double dist2 = dx * dx + dy * dy;

		if (dist2 < minDist2) {
			minDist2 = dist2;
			bestT = t;
		}

	}

	// f'(t)
	auto f = [&](double t) {
		double x = spline.getValue(axis::x, 0, t);
		double y = spline.getValue(axis::y, 0, t);
		double dx = x - point.x;
		double dy = y - point.y;

		double xp = spline.getValue(axis::x, 1, t);
		double yp = spline.getValue(axis::y, 1, t);

		return 2 * (dx * xp + dy * yp);
		};

	// f''(t)
	auto df = [&](double t) {
		double x = spline.getValue(axis::x, 0, t);
		double y = spline.getValue(axis::y, 0, t);

		double dx = x - point.x;
		double dy = y - point.y;

		double xp = spline.getValue(axis::x, 1, t);
		double yp = spline.getValue(axis::y, 1, t);

		double xpp = spline.getValue(axis::x, 2, t);
		double ypp = spline.getValue(axis::y, 2, t);

		return 2 * (xp * xp + dx * xpp + yp * yp + dy * ypp); // f''(t)
		};

	double finalT = bestT;

	// --- шаг 2: уточнение методом Ньютона ---
	if (finalT != 0.0 && finalT != 1.0)
		finalT = NumericalMethods::newtonMethod(f, df, bestT, 0.0, 1.0);

	double xFinal = spline.getValue(axis::x, 0, finalT);
	double yFinal = spline.getValue(axis::y, 0, finalT);
	double dist2 = (xFinal - point.x) * (xFinal - point.x) +
		(yFinal - point.y) * (yFinal - point.y);

	return { Point{xFinal,yFinal}, sqrt(dist2) };
}


//-----------------------------
// Пересечение сплайнов
//---
std::vector<Point> GeometricOperations::findIntersection(const Spline& s1, const Spline& s2, bool analit) {

	auto F = [&](const Vec& v) {

		// F1  = x1(t) - x2(u) = 0
		// F2  = y1(t) - y2(u) = 0

		double t = v[0];
		double u = v[1];

		double x1 = s1.getValue(axis::x, 0, t);
		double y1 = s1.getValue(axis::y, 0, t);
		double x2 = s2.getValue(axis::x, 0, u);
		double y2 = s2.getValue(axis::y, 0, u);

		Vec res(2);
		res[0] = x1 - x2;   // F1
		res[1] = y1 - y2;   // F2
		return res;
		};
	// Якобиан
	std::function<Mat(const Vec&)> J;
	if (analit)
		J = [&](const Vec& v) -> Mat {

		// F1  = x1(t) - x2(u) = 0
		// F2  = y1(t) - y2(u) = 0
		// 
		// dF1/dt = x1'
		// dF1/du = -x2'
		// dF2/dt = y1'
		// dF2/du = -y2'

		double t = v[0];
		double u = v[1];

		double x1 = s1.getValue(axis::x, 0, t);
		double y1 = s1.getValue(axis::y, 0, t);
		double x2 = s2.getValue(axis::x, 0, u);
		double y2 = s2.getValue(axis::y, 0, u);

		double dx = x1 - x2;
		double dy = y1 - y2;

		double x1p = s1.getValue(axis::x, 1, t);
		double y1p = s1.getValue(axis::y, 1, t);
		double x2p = s2.getValue(axis::x, 1, u);
		double y2p = s2.getValue(axis::y, 1, u);

		Mat J(2, Vec(2));

		// dF1/dt
		J[0][0] = x1p;
		// dF1/du
		J[0][1] = -x2p;
		// dF2/dt
		J[1][0] = y1p;
		// dF2/du
		J[1][1] = -y2p;

		return J;
		};
	else
		J = [&](const Vec& v) -> Mat {
		return NumericalMethods::numericalJacobian(F, v);
		};




	// --- грубый перебор ---
	int nGrid = 25;
	std::vector<Point> intersections;

	Vec minBounds = { 0.0, 0.0 };
	Vec maxBounds = { 1.0, 1.0 };

	double t = 0.0;
	double xt = s1.getValue(axis::x, 0, t);
	double yt = s1.getValue(axis::y, 0, t);

	for (int i = 1; i <= nGrid; ++i) {

		double t1 = double(i) / nGrid;
		double xt1 = s1.getValue(axis::x, 0, t1);
		double yt1 = s1.getValue(axis::y, 0, t1);

		double u = 0.0;
		double xu = s2.getValue(axis::x, 0, u);
		double yu = s2.getValue(axis::y, 0, u);

		for (int j = 1; j <= nGrid; ++j) {

			double u1 = double(j) / nGrid;
			double xu1 = s2.getValue(axis::x, 0, u1);
			double yu1 = s2.getValue(axis::y, 0, u1);

			auto inter = findIntersection(Point(xt, yt), Point(xt1, yt1), Point(xu, yu), Point(xu1, yu1));

			// --- Уточнение методом Ньютона---
			if (inter.size() > 0) {

				Vec start = { (t + t1) / 2, (u + u1) / 2 };

				try {
					auto sol = NumericalMethods::newtonMethod(F, J, start, minBounds, maxBounds);
					double tFinal = sol[0];
					double uFinal = sol[1];
					Point p1 = { 
									s1.getValue(axis::x, 0, tFinal),
									s1.getValue(axis::y, 0, tFinal),
								};
					Point p2 = { 
									s2.getValue(axis::x, 0, uFinal),
									s2.getValue(axis::y, 0, uFinal),
								};

					intersections.push_back(p1);
				}
				catch (...) {
					// ничего не делаем
				}

			}

			u = u1;
			xu = xu1;
			yu = yu1;

		}
		t = t1;
		xt = xt1;
		yt = yt1;

	}

	return intersections;
}