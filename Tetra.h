//
// Created by i3alumba on 02.06.23.
//

#ifndef RTX_TETRA_H
#define RTX_TETRA_H

#include "Vec3.h"
#include "Color.h"
#include "Figure.h"

struct Tetra : public Figure {
    Vec3 a;
    Vec3 b;
    Vec3 c;
    Vec3 d;

    using Figure::color;

    Tetra(Color color, Vec3 a, Vec3 b, Vec3 c, Vec3 d) : Figure(color), a(a), b(b), c(c), d(d) {}
};


#endif //RTX_TETRA_H
