#include "Shapes.h"

//-----------------------------
// Заполнить сплайн (решение сисстемы и нахождение производных)
//---
void Spline::buildSpline(const std::vector<double>& y, std::vector<SplineSegment>& out) {
	if (n < 2) {
		out.clear();
		return;
	}

	// шаги по t
	std::vector<double> h(n - 1);
	for (size_t i = 0; i < n - 1; ++i)
		h[i] = t[i + 1] - t[i];

	// коэффициенты для трёхдиагональной системы
	std::vector<double> a(n), b(n), c(n), d(n), m(n);
	b[0] = 1.0; d[0] = 0.0; // граничные условия
	for (size_t i = 1; i < n - 1; ++i) {
		a[i] = h[i - 1];
		b[i] = 2.0 * (h[i - 1] + h[i]);
		c[i] = h[i];
		d[i] = 6.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
	}
	b[n - 1] = 1.0;
	d[n - 1] = 0.0;

	// прогонка
	std::vector<double> c_prim(n), d_prim(n);
	c_prim[0] = c[0] / b[0];
	d_prim[0] = d[0] / b[0];
	for (size_t i = 1; i < n; ++i) {
		double denom = b[i] - a[i] * c_prim[i - 1];
		if (i < n - 1) c_prim[i] = c[i] / denom;
		d_prim[i] = (d[i] - a[i] * d_prim[i - 1]) / denom;
	}

	m[n - 1] = d_prim[n - 1];
	for (int i = n - 2; i >= 0; --i)
		m[i] = d_prim[i] - c_prim[i] * m[i + 1];

	out.resize(n - 1);
	for (size_t i = 0; i < n - 1; ++i) {
		out[i].t0 = t[i];
		out[i].t1 = t[i + 1];
		out[i].y0 = y[i];
		out[i].y1 = y[i + 1];
		out[i].m0 = m[i];
		out[i].m1 = m[i + 1];
	}
}


//-----------------------------
// Конструктор
//---
Spline::Spline() :n(0) {};


//-----------------------------
// Конструктор
//---
Spline::Spline(std::vector<Point>& points) :n(points.size()) {
	setBasePoints(points);
};


//-----------------------------
// Создать сплайн по базовым точкам
//---
void Spline::setBasePoints(const std::vector<Point>& points) {
	n = points.size();

	// 1. Параметризация по длине дуги
	t.resize(n);
	t[0] = 0.0;
	for (size_t i = 1; i < n; ++i) {
		double dx = points[i].x - points[i - 1].x;
		double dy = points[i].y - points[i - 1].y;
		t[i] = t[i - 1] + sqrt(dx * dx + dy * dy);
	}

	// 2. Нормализация
	double total = t.back();
	if (total > 0.0) {
		for (size_t i = 0; i < n; ++i)
			t[i] /= total;
	}

	// 3. Построение сплайна
	std::vector<double> xs(n), ys(n);
	for (size_t i = 0; i < n; ++i) {
		xs[i] = points[i].x;
		ys[i] = points[i].y;
	}

	buildSpline(xs, splineX);
	buildSpline(ys, splineY);
}


//-----------------------------
// Получить знаечние
//---
double Spline::getValue(axis ax, int der, double tt) const {

	std::vector<SplineSegment> splineAx;

	switch (ax)
	{
	case axis::x: splineAx = splineX; break;
	case axis::y: splineAx = splineY; break;
	}

	if (splineAx.empty())
		throw std::runtime_error("Error constructing spline");

	// ищем сегмент соответствующий параметру tt

	SplineSegment seg;
	bool found = false;
	for (size_t i = 0; i < splineAx.size(); ++i) {
		if (tt <= splineAx[i].t1) {
			seg = splineAx[i];
			found = true;
			break;
		}
	}
	if (!found)
		seg = splineAx.back(); // если tt == 1.0 или чуть больше из-за ошибок округления


	// находим значение или производную(1,2) в точке 

	double h = seg.t1 - seg.t0;
	double dx1 = seg.t1 - tt;
	double dx2 = tt - seg.t0;

	double out = 0;

	switch (der)
	{
	case 0:
		out = seg.m0 * dx1 * dx1 * dx1 / (6 * h)
			+ seg.m1 * dx2 * dx2 * dx2 / (6 * h)
			+ (seg.y0 - seg.m0 * h * h / 6) * (dx1 / h)
			+ (seg.y1 - seg.m1 * h * h / 6) * (dx2 / h);
		break;
	case 1:
		out = -seg.m0 * dx1 * dx1 / (2 * h)
			+ seg.m1 * dx2 * dx2 / (2 * h)
			- (seg.y0 - seg.m0 * h * h / 6.0) / h
			+ (seg.y1 - seg.m1 * h * h / 6.0) / h;
		break;
	case 2:
		out = seg.m0 * dx1 / h + seg.m1 * dx2 / h;
		break;
	default:
		throw std::invalid_argument("no derivative num");
	}

	return out;

}
