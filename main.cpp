#include <iostream>
#include <cmath>


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
        int d = int(Direction / 360);
        Direction -= d * 360;
        if (Direction < 0.0)
            Direction += 360;
    }
};

int PointsOfIntersectionOfTwoSpheres(const TSphere &s0, const TSphere &s1, double r, double &t0, double &t1) {
    t0 = t1 = 0.0;
    double x0 = s0.Position.X, y0 = s0.Position.Y, u0 = cos(s0.Direction) * s0.Velocity, v0 = sin(s0.Direction) * s0.Velocity;
    double x1 = s1.Position.X, y1 = s1.Position.Y, u1 = cos(s1.Direction) * s1.Velocity, v1 = sin(s1.Direction) * s1.Velocity;
    double a = u0 * u0 + u1 * u1 + v0 * v0 + v1 * v1 - 2 * u0 * u1 - 2 * v0 * v1;
    double b = 2 * x0 * u0 - 2 * x0 * u1 - 2 * u0 * x1 + 2 * x1 * u1 + 2 * y0 * v0 - 2 * y0 * v1 - 2 * v0 * y1 + 2 * y1 * v1;
    double c = x0 * x0 + x1 * x1  + y0 * y0 + y1 * y1 - 2 * x0 * x1 - 2 * y0 * y1 - 4 * r * r;
    std::cout << a << ' ' << b << ' ' << c << std::endl;
    return 0;
}

int main() {
    TSphere s0(TPoint(0, 0), 0, 0);
    TSphere s1(TPoint(0, 0), 0, 1);
    double t0, t1;
    PointsOfIntersectionOfTwoSpheres(s0, s1, 1, t0, t1);
    return 0;
}

