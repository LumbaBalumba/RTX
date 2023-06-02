//
// Created by i3alumba on 02.06.23.
//

#include "Box.h"

bool Box::contains(Point p) const {
    return p.x >= point.x - EPS && p.x <= point.x + a + EPS && p.y >= point.y - EPS &&
           p.y <= point.y + b + EPS && p.z >= point.z - EPS &&
           p.z <= point.z + c + EPS;
}
