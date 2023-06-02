//
// Created by i3alumba on 02.06.23.
//

#ifndef RTX_BOX_H
#define RTX_BOX_H

#include "Vec3.h"
#include "Color.h"
#include "Figure.h"

struct Box : public Figure {
    Vec3 point;
    double a;
    double b;
    double c;

    using Figure::color;

    Box(Color color, Vec3 point, double a, double b, double c) : Figure(color), point(point), a(a), b(b), c(c) {}

    [[nodiscard]] bool contains(Point p) const ;
};


#endif //RTX_BOX_H
