//
// Created by i3alumba on 02.06.23.
//

#include "Vec3.h"

double Vec3::dot(Vec3 left, Vec3 right) {
    return left.x * right.x + left.y * right.y + left.z * right.z;
}

Vec3 Vec3::cross(Vec3 left, Vec3 right) {
    return {left.y * right.z - left.z * right.y, left.z * right.x - right.z * left.x,
            left.x * right.y - left.y * right.x};
}

Vec3 operator+(Vec3 left, Vec3 right) {
    return {left.x + right.x, left.y + right.y, left.z + right.z};
}

Vec3 operator*(Vec3 left, Vec3 right) {
    return {left.x * right.x, left.y * right.y, left.z * right.z};
}

Vec3 Vec3::operator-() const {
    return {-x, -y, -z};
}

Vec3 operator-(Vec3 left, Vec3 right) {
    return left + -right;
}

Vec3 operator/(Vec3 left, Vec3 right) {
    return {left.x / right.x, left.y / right.y, left.z / right.z};
}

std::istream &operator>>(std::istream &in, Vec3 &v) {
    return in >> v.x >> v.y >> v.z;
}

std::ostream &operator<<(std::ostream &out, const Vec3 &v) {
    return out << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}

double Vec3::length() const {
    return std::sqrt(x * x + y * y + z * z);
}

double Vec3::distance(const Vec3 &left, const Vec3 &right) {
    return (left - right).length();
}

Vec3 Vec3::normalized() const {
    return *this / length();
}

bool Vec3::operator==(const Vec3 &other) const {
    if (std::abs(x - other.x) > 1e-6) {
        return false;
    }
    if (std::abs(y - other.y) > 1e-6) {
        return false;
    }
    if (std::abs(z - other.z) > 1e-6) {
        return false;
    }
    return true;
}

bool Vec3::operator!=(const Vec3 &other) const {
    return !(*this == other);
}

