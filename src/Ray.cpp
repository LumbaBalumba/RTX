//
// Created by i3alumba on 02.06.23.
//

#include "Ray.h"

Point Ray::plane_intersection(Point a, Point b, Point c) const {
    Vec3 normal = Vec3::cross(b - a, c - a);
    double dot = Vec3::dot(direction, normal);
    if (std::abs(dot) < EPS) {
        return point;
    }
    Vec3 displacement = a - point;
    double t = Vec3::dot(displacement, normal) / dot;
    if (t < 0) {
        return point;
    }
    return point + t * direction;
}

Point Ray::triangle_intersection(Point a, Point b, Point c) const {
    Point intersection_point = plane_intersection(a, b, c);
    if (intersection_point == point) {
        return point;
    }
    Vec3 edge1 = b - a;
    Vec3 edge2 = c - a;
    Vec3 norm = Vec3::cross(edge1, edge2);

    Vec3 pa = a - intersection_point;
    Vec3 pb = b - intersection_point;
    Vec3 pc = c - intersection_point;

    Vec3 ppa = Vec3::cross(pa, pb);
    Vec3 ppb = Vec3::cross(pb, pc);
    Vec3 ppc = Vec3::cross(pc, pa);

    double p1 = Vec3::dot(norm, ppa);
    double p2 = Vec3::dot(norm, ppb);
    double p3 = Vec3::dot(norm, ppc);
    if (p1 >= 0 && p2 >= 0 && p3 >= 0) {
        return intersection_point;
    } else {
        return point;
    }
}

std::pair<Point, Color> Ray::intersect(const Figure *fig, Vec3 light_source) const {
    auto *s_p = dynamic_cast<const Sphere *>(fig);
    if (s_p != nullptr) {
        Sphere s = *s_p;
        double a = pow(direction.x, 2) + pow(direction.y, 2) + pow(direction.z, 2);
        double b = 2 * (direction.x * (point.x - s.center.x) + direction.y * (point.y - s.center.y) +
                        direction.z * (point.z - s.center.z));
        double c = pow(point.x - s.center.x, 2) + pow(point.y - s.center.y, 2) + pow(point.z - s.center.z, 2) -
                   pow(s.radius, 2);

        double discriminant = pow(b, 2) - 4 * a * c;

        if (discriminant < 0) {
            return {point, {255, 255, 255}};
        }

        double t1 = (-b + sqrt(discriminant)) / (2 * a);
        double t2 = (-b - sqrt(discriminant)) / (2 * a);

        double t;
        if (t2 > 0) {
            t = t2;
        } else if (t2 <= 0 && t1 > 0) {
            t = t1;
        } else {
            return {point, {255, 255, 255}};
        }
        if (std::abs(t) < EPS) {
            return {point, {255, 255, 255}};
        }
        Point intersection_point = {
                point.x + t * direction.x,
                point.y + t * direction.y,
                point.z + t * direction.z
        };

        return {intersection_point, s.color() * Vec3::dot((light_source - intersection_point).normalized(),
                                                          (intersection_point - s.center).normalized())};
    }
    auto *b_p = dynamic_cast<const Box *>(fig);
    if (b_p != nullptr) {
        Box b = *b_p;
        std::vector<double> params;
        double x1 = b.point.x, x2 = b.point.x + b.a;
        double y1 = b.point.y, y2 = b.point.y + b.b;
        double z1 = b.point.z, z2 = b.point.z + b.c;
        if (std::abs(direction.x) >= EPS) {
            double t = (x1 - point.x) / direction.x;
            Point p = point + t * direction;
            if (t > 0 && b.contains(p)) {
                params.push_back(t);
            }
            t = (x2 - point.x) / direction.x;
            p = point + t * direction;
            if (t > 0 && b.contains(p)) {
                params.push_back(t);
            }
        }
        if (std::abs(direction.y) >= EPS) {
            double t = (y1 - point.y) / direction.y;
            Point p = point + t * direction;
            if (t > 0 && b.contains(p)) {
                params.push_back(t);
            }
            t = (y2 - point.y) / direction.y;
            p = point + t * direction;
            if (t > 0 && b.contains(p)) {
                params.push_back(t);
            }
        }
        if (std::abs(direction.z) >= EPS) {
            double t = (z1 - point.z) / direction.z;
            Point p = point + t * direction;
            if (t > 0 && b.contains(p)) {
                params.push_back(t);
            }
            t = (z2 - point.z) / direction.z;
            p = point + t * direction;
            if (t > 0 && b.contains(p)) {
                params.push_back(t);
            }
        }
        if (params.empty()) {
            return {point, {255, 225, 255}};
        }
        std::sort(params.begin(), params.end());
        double t = params[0];
        Point ref_point = point + t * direction;
        Vec3 norm{};
        if (std::abs(ref_point.x - b.point.x) < EPS) {
            norm = {-1, 0, 0};
        } else if (std::abs(ref_point.x - b.point.x - b.a) < EPS) {
            norm = {1, 0, 0};
        } else if (std::abs(ref_point.y - b.point.y) < EPS) {
            norm = {0, -1, 0};
        } else if (std::abs(ref_point.y - b.point.y - b.b) < EPS) {
            norm = {0, 1, 0};
        } else if (std::abs(ref_point.z - b.point.z) < EPS) {
            norm = {0, 0, -1};
        } else if (std::abs(ref_point.z - b.point.z - b.c) < EPS) {
            norm = {0, 0, 1};
        }
        return {ref_point, b.color() * Vec3::dot(norm.normalized(), (light_source - ref_point).normalized())};
    }
    auto *t_p = dynamic_cast<const Tetra *>(fig);
    if (t_p != nullptr) {
        auto t = *t_p;
        Point intersection_point = point;
        Color intersection_color = {255, 255, 255};
        std::array<Point, 4> points = {t.a, t.b, t.c, t.d};
        for (int i = 0; i < 4; ++i) {
            Point point1 = points[i];
            Point point2 = points[(i + 1) % 4];
            Point point3 = points[(i + 2) % 4];
            Point tmp = triangle_intersection(point1, point2, point3);
            if (tmp == point) {
                continue;
            } else if (intersection_point == point ||
                       Vec3::distance(tmp, point) < Vec3::distance(intersection_point, point)) {
                Vec3 v1 = point2 - point1;
                Vec3 v2 = point3 - point1;
                Vec3 norm = Vec3::cross(v1, v2);
                intersection_point = tmp;
                intersection_color =
                        t.color() * Vec3::dot(norm.normalized(), (light_source - intersection_point).normalized());
            }
        }
        return {intersection_point, intersection_color};
    }
    throw "Invalid figure type";
}