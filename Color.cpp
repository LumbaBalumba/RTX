//
// Created by i3alumba on 02.06.23.
//

#include "Color.h"

Color operator+(Color left, Color right) {
    return {left.r + right.r, left.g + right.g, left.b + right.b};
}

Color operator*(Color left, Color right) {
    return {std::min(255.0, std::max(left.r * right.r, 0.0)),
            std::min(255.0, std::max(left.g * right.g, 0.0)),
            std::min(255.0, std::max(left.b * right.b, 0.0))};
}

Color Color::operator-() const {
    return {-r, -g, -b};
}

Color operator-(Color left, Color right) {
    return left + -right;
}

Color operator/(Color left, Color right) {
    return {left.r / right.r, left.g / right.g, left.b / right.b};
}

std::istream &operator>>(std::istream &in, Color &v) {
    return in >> v.r >> v.g >> v.b;
}

std::ostream &operator<<(std::ostream &out, const Color &v) {
    return out << "(" << v.r << ", " << v.g << ", " << v.b << ")";
}

