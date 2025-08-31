#include <fstream>
#include <iomanip>

#include "Shapes.h"
#include "NumericalMethods.h"
#include "GeometricOperations.h"
#include "tests.h"


//-----------------------------
// Прочиать точки сплайна из файла
//---
std::vector<Point> readSpline(std::istream& in) {
    int n;
    in >> n;
    if (n < 2) throw std::runtime_error("Spline size < 2");

    std::vector<Point> pts(n);
    for (int i = 0; i < n; ++i) {
        in >> pts[i].x >> pts[i].y;
    }
    return pts;
}


//-----------------------------
// Вывести информацию о расстоянии между сплайнами
//---
void printDistance(std::ostream& out, const Point& p1, const Point& p2, double dist) {
    out << "Distance\n";
    out << std::fixed << std::setprecision(9);
    out << "p1 = (" << p1.x << ", " << p1.y << ")\n";
    out << "p2 = (" << p2.x << ", " << p2.y << ")\n";
    out << "d  = " << dist << "\n";
}


//-----------------------------
// Вывести точки пересечения между сплайнами
//---
void printIntersections(std::ostream& out, const std::vector<Point>& inters) {
    out << "Intersection\n";
    out << std::fixed << std::setprecision(9);
    for (size_t i = 0; i < inters.size(); ++i) {
        out << "p" << i << " = (" << inters[i].x << ", " << inters[i].y << ")\n";
    }
}


//-----------------------------
// Ищем расстояние между сплайнами, если расстояние = 0, ищем пересечения
// Случай касания сплайнов не учитывается
//---
int main() {
	
	// Тесты
	// runTestDistanceSplinePoint();
	 runTestDistanceSplineSpline();
	// runTestIntersectionSplineSpline();

    std::ifstream fin("data/in.txt");
    std::ofstream fout("data/out.txt");

    try {
        auto pts1 = readSpline(fin);
        auto pts2 = readSpline(fin);

        Spline spline1(pts1);
        Spline spline2(pts2);

        auto [p1, p2, dist] = GeometricOperations::distance(spline1, spline2);

        if (dist > 1e-9) {
            printDistance(fout, p1, p2, dist);
            return 0;
        }

        auto inters = GeometricOperations::findIntersection(spline1, spline2);
        if (inters.empty()) {
            fout << "No intersections found\n";
            return -1;
        }
        printIntersections(fout, inters);
        return 0;
    }
    catch (std::exception& e) {
        fout << "Error: " << e.what() << "\n";
        return -1;
    }
}