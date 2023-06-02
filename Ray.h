//
// Created by i3alumba on 02.06.23.
//

#ifndef RTX_RAY_H
#define RTX_RAY_H

#include <vector>
#include <array>

#include "Vec3.h"
#include "Color.h"
#include "Figure.h"
#include "Sphere.h"
#include "Box.h"
#include "Tetra.h"

struct Ray {
    Vec3 point;
    Vec3 direction;

    Ray(Vec3 point, Vec3 direction) : point(point), direction(direction) {}

    [[nodiscard]] Point plane_intersection(Point a, Point b, Point c) const;

    [[nodiscard]] Point triangle_intersection(Point a, Point b, Point c) const;

    std::pair<Point, Color> intersect(const Figure *fig, Vec3 light_source) const;
};

#endif //RTX_RAY_H
