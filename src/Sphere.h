//
// Created by i3alumba on 02.06.23.
//

#ifndef RTX_SPHERE_H
#define RTX_SPHERE_H

#include "Vec3.h"
#include "Color.h"
#include "Figure.h"

struct Sphere : public Figure {
    Vec3 center;
    double radius;

    using Figure::color;

    Sphere(Color c, Vec3 center, double radius) : Figure(c), center(center), radius(radius) {}
};

#endif //RTX_SPHERE_H
