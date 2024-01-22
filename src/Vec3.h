//
// Created by i3alumba on 02.06.23.
//

#ifndef RTX_VEC3_H
#define RTX_VEC3_H

#include <iostream>
#include <cmath>

const double EPS = 1e-6;

struct Vec3 {
    double x;
    double y;
    double z;

    Vec3() = default;

    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    Vec3(double value) : x(value), y(value), z(value) {}

    static double dot(Vec3 left, Vec3 right);

    static Vec3 cross(Vec3 left, Vec3 right);

    Vec3 &operator=(const Vec3 &other) = default;

    friend Vec3 operator+(Vec3 left, Vec3 right);

    friend Vec3 operator*(Vec3 left, Vec3 right);

    Vec3 operator-() const;

    friend Vec3 operator-(Vec3 left, Vec3 right);

    friend Vec3 operator/(Vec3 left, Vec3 right);

    friend std::istream &operator>>(std::istream &in, Vec3 &v);

    friend std::ostream &operator<<(std::ostream &out, const Vec3 &v);

    [[nodiscard]] double length() const;

    static double distance(const Vec3 &left, const Vec3 &right);

    [[nodiscard]] Vec3 normalized() const;

    bool operator==(const Vec3 &other) const;

    bool operator!=(const Vec3 &other) const;
};

typedef Vec3 Point;


#endif //RTX_VEC3_H
