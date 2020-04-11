#include <iostream>
#include <iomanip>
#include <cmath>


constexpr double EPS = 1e-7;

struct TPoint {
    double X = 0.0, Y = 0.0;
    TPoint() = default;
    TPoint(double x, double y)
        : X(x), Y(y) {
    }
};

struct TSphere {
    TPoint Position;
    double Direction = 0.0;
    double Velocity = 0.0;

    TSphere() = default;
    TSphere(const TPoint &position, double direction, double velocity)
        : Position(position), Direction(direction), Velocity(velocity) {
        int d = int(Direction / 2 / M_PI);
        Direction -= d * (2 * M_PI);
        if (Direction < 0.0)
            Direction += 2 * M_PI;
    }
};

int PointsOfIntersectionOfTwoSpheres(const TSphere &s0, const TSphere &s1, double r, double &t0, double &t1) {
    t0 = t1 = 0.0;
    double x0 = s0.Position.X, y0 = s0.Position.Y, u0 = cos(s0.Direction) * s0.Velocity, v0 = sin(s0.Direction) * s0.Velocity;
    double x1 = s1.Position.X, y1 = s1.Position.Y, u1 = cos(s1.Direction) * s1.Velocity, v1 = sin(s1.Direction) * s1.Velocity;
    double a = u0 * u0 + u1 * u1 + v0 * v0 + v1 * v1 - 2 * u0 * u1 - 2 * v0 * v1;
    double b = 2 * x0 * u0 - 2 * x0 * u1 - 2 * u0 * x1 + 2 * x1 * u1 + 2 * y0 * v0 - 2 * y0 * v1 - 2 * v0 * y1 + 2 * y1 * v1;
    double c = x0 * x0 + x1 * x1  + y0 * y0 + y1 * y1 - 2 * x0 * x1 - 2 * y0 * y1 - 4 * r * r;
    //std::cout << a << ' '<< b << ' ' << c << std::endl;
    if (fabs(a) < EPS) {
        if (fabs(b) < EPS) {
            if (fabs(c) < EPS)
                return 1;
            else
                return 0;
        } else {
            t0 = -c / b;
            return 1;
        }
    }
    double d2 = b * b - 4 * a * c;
    //std::cout << d2 << std::endl;
    if (d2 < -EPS)
        return 0;
    if (d2 < EPS) {
        t0 = -b / 2 / a;
        return 1;
    }
    double d = sqrt(d2);
    //std::cout << d << std::endl;
    t0 = (-b - d) / 2 / a;
    t1 = (-b + d) / 2 / a;
    return 2;
}

bool TestIntersection1() {
    TSphere s0(TPoint(0, 0), 0, 0);
    TSphere s1(TPoint(0, 0), 0, 1);
    double t0, t1;
    int points = PointsOfIntersectionOfTwoSpheres(s0, s1, 1, t0, t1);
    std::cout << std::fixed << std::setprecision(10) << points << ' ' << t0 << ' ' << t1 << std::endl;
    if (points != 2)
        return false;
    if (fabs(t0 + 2) > EPS || fabs(t1 - 2) > EPS)
        return false;
    return true;
}

bool TestIntersection2() {
    TSphere s0(TPoint(-2, 0), M_PI / 4, 1);
    TSphere s1(TPoint(2, 0), M_PI * 3 / 4, 1);
    double t0, t1;
    int points = PointsOfIntersectionOfTwoSpheres(s0, s1, 1, t0, t1);
    std::cout << std::fixed << std::setprecision(10) << points << ' ' << t0 << ' ' << t1 << std::endl;
    if (points != 2)
        return false;
    if (fabs(t0 - 1.4142135624) > EPS || fabs(t1 - 4.2426406871) > EPS)
        return false;
    return true;
}

int main() {
    //std::cout << TestIntersection1() << std::endl;
    //std::cout << TestIntersection2() << std::endl;
    //double m1 = 1.0, m2 = 1.0, x1 = 2.0, y1 = 1.5, x2 = 0.3, y2 = -y1 * m1 / m2;
    double m1 = 1.0, m2 = 1.0, x1 = 1.0, y1 = 1.0, x2 = 1.0, y2 = -y1 * m1 / m2;
    double e = m1 * (x1 * x1 + y1 * y1) + m2 * (x2 * x2 + y2 * y2), p = m1 * x1 + m2 * x2;
    double a = m1 + m1 * m1 / m2;
    double b = -2 * p * m1 / m2;
    double c = m1 * y1 * y1 + p * p / m2 + y1 * y1 * m1 * m1 / m2 - e;
    double d = sqrt(b * b - 4 * a * c);
    std::cout << a << ' ' << b << ' ' << c << ' ' << d << std::endl;
    double t1 = (-b - d) / 2 / a, t2 = (-b + d) / 2 / a;
    std::cout << x1 << ' ' << y1 << ' ' << (p - x1 * m1) / m2 << ' ' << -y1 * m1 / m2 << std::endl;
    std::cout << t1 << ' ' << y1 << ' ' << (p - t1 * m1) / m2 << ' ' << -y1 * m1 / m2 << std::endl;
    std::cout << t2 << ' ' << y1 << ' ' << (p - t2 * m1) / m2 << ' ' << -y1 * m1 / m2 << std::endl;
    return 0;
}

