//
// Created by i3alumba on 02.06.23.
//

#ifndef RTX_COLOR_H
#define RTX_COLOR_H

#include <algorithm>
#include <iostream>

struct Color {
    double r;
    double g;
    double b;

    Color() = default;

    Color(double x, double y, double z) : r(x), g(y), b(z) {}

    Color(double value) : r(value), g(value), b(value) {}


    Color &operator=(const Color &other) = default;

    friend Color operator+(Color left, Color right);

    friend Color operator*(Color left, Color right);

    Color operator-() const;

    friend Color operator-(Color left, Color right);

    friend Color operator/(Color left, Color right);

    friend std::istream &operator>>(std::istream &in, Color &v);

    friend std::ostream &operator<<(std::ostream &out, const Color &v);
};


#endif //RTX_COLOR_H
