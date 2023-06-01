#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <SFML/Graphics.hpp>
#include <limits>
#include <omp.h>
#include <algorithm>

/*
Аргументы запуска программы - файл с описанием 3D-пространства, файл с описанием фигур, название выходного файла.
Ввод в первом файле осуществляется так же, как описано в ТЗ, за исключением того, что не требуется вводназвания параметров.
Во втором файле вводятся фигуры по следующему формату:
    1) первым вводится вид фигуры
    2) затем вводятся параметры формы (для сферы - центр и радиус, для параллелепипеда - точка и три стороны, для тетраэдра - четыре точки)
    3) последним вводится цвет фигуры, три числа, каждое от 0 до 255
 */


struct Vec3 {
    double x;
    double y;
    double z;

    Vec3() = default;

    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    Vec3(double value) : x(value), y(value), z(value) {}

    static double dot(Vec3 left, Vec3 right) {
        return left.x * right.x + left.y * right.y + left.z * right.z;
    }

    static Vec3 cross(Vec3 left, Vec3 right) {
        return {left.y * right.z - left.z * right.y, left.z * right.x - right.z * left.x,
                left.x * right.y - left.y * right.x};
    }

    Vec3 &operator=(const Vec3 &other) = default;

    friend Vec3 operator+(Vec3 left, Vec3 right) {
        return {left.x + right.x, left.y + right.y, left.z + right.z};
    }

    friend Vec3 operator*(Vec3 left, Vec3 right) {
        return {left.x * right.x, left.y * right.y, left.z * right.z};
    }

    Vec3 operator-() const {
        return {-x, -y, -z};
    }

    friend Vec3 operator-(Vec3 left, Vec3 right) {
        return left + -right;
    }

    friend Vec3 operator/(Vec3 left, Vec3 right) {
        return {left.x / right.x, left.y / right.y, left.z / right.z};
    }

    friend std::istream &operator>>(std::istream &in, Vec3 &v) {
        return in >> v.x >> v.y >> v.z;
    }

    friend std::ostream &operator<<(std::ostream &out, const Vec3 &v) {
        return out << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    }

    [[nodiscard]] double length() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    static double distance(const Vec3 &left, const Vec3 &right) {
        return (left - right).length();
    }

    [[nodiscard]] Vec3 normalized() const {
        return *this / length();
    }

    bool operator==(const Vec3 &other) const {
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

    bool operator!=(const Vec3 &other) const {
        return !(*this == other);
    }
};

typedef Vec3 Point;

struct Color {
    double r;
    double g;
    double b;

    Color() = default;

    Color(double x, double y, double z) : r(x), g(y), b(z) {}

    Color(double value) : r(value), g(value), b(value) {}


    Color &operator=(const Color &other) = default;

    friend Color operator+(Color left, Color right) {
        return {left.r + right.r, left.g + right.g, left.b + right.b};
    }

    friend Color operator*(Color left, Color right) {
        return {std::min(255.0, std::max(left.r * right.r, 0.0)),
                std::min(255.0, std::max(left.g * right.g, 0.0)),
                std::min(255.0, std::max(left.b * right.b, 0.0))};
    }

    Color operator-() const {
        return {-r, -g, -b};
    }

    friend Color operator-(Color left, Color right) {
        return left + -right;
    }

    friend Color operator/(Color left, Color right) {
        return {left.r / right.r, left.g / right.g, left.b / right.b};
    }

    friend std::istream &operator>>(std::istream &in, Color &v) {
        return in >> v.r >> v.g >> v.b;
    }

    friend std::ostream &operator<<(std::ostream &out, const Color &v) {
        return out << "(" << v.r << ", " << v.g << ", " << v.b << ")";
    }
};

struct Figure {
    Color c;

    explicit Figure(Color c) : c(c) {}

    virtual ~Figure() = default;

    Color color() {
        return c;
    }
};

struct Sphere : public Figure {
    Vec3 center;
    double radius;

    using Figure::color;

    Sphere(Color c, Vec3 center, double radius) : Figure(c), center(center), radius(radius) {}
};

struct Box : public Figure {
    Vec3 point;
    double a;
    double b;
    double c;

    using Figure::color;

    Box(Color color, Vec3 point, double a, double b, double c) : Figure(color), point(point), a(a), b(b), c(c) {}

    [[nodiscard]] bool contains(Point p) const {
        return p.x >= point.x && p.x <= point.x + a && p.y >= point.y && p.y <= point.y + b && p.z >= point.z &&
               p.z <= point.z + c;
    }
};

struct Tetrahedron : public Figure {
    Vec3 a;
    Vec3 b;
    Vec3 c;
    Vec3 d;

    using Figure::color;

    Tetrahedron(Color color, Vec3 a, Vec3 b, Vec3 c, Vec3 d) : Figure(color), a(a), b(b), c(c), d(d) {}
};

struct Ray {
    Vec3 point;
    Vec3 direction;

    Ray(Vec3 point, Vec3 direction) : point(point), direction(direction) {}

    [[nodiscard]] Point plane_intersection(Point a, Point b, Point c) const {
        Vec3 normal = Vec3::cross(b - a, c - a);
        double dot = Vec3::dot(direction, normal);
        if (std::abs(dot) < 1e-6) {
            return point;
        }
        Vec3 displacement = a - point;
        double t = Vec3::dot(displacement, normal) / dot;
        if (t < 0) {
            return point;
        }
        return point + t * direction;
    }

    [[nodiscard]] Point triangle_intersection(Point a, Point b, Point c) const {
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

    std::pair<Point, Color> intersect(const Figure *fig, Vec3 light_source) const {
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
            if (std::abs(t) < 1e-6) {
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
            if (std::abs(direction.x) >= 1e-6) {
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
            if (std::abs(direction.y) >= 1e-6) {
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
            if (std::abs(direction.z) >= 1e-6) {
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
                return {point, {255, 255, 255}};
            }
            std::sort(params.begin(), params.end());
            double t = params[0];
            Point ref_point = point + t * direction;
            Vec3 norm{};
            if (std::abs(ref_point.x - b.point.x) < 1e-6) {
                norm = {-1, 0, 0};
            } else if (std::abs(ref_point.x - b.point.x - b.a) < 1e-6) {
                norm = {1, 0, 0};
            } else if (std::abs(ref_point.y - b.point.y) < 1e-6) {
                norm = {0, -1, 0};
            } else if (std::abs(ref_point.y - b.point.y - b.b) < 1e-6) {
                norm = {0, 1, 0};
            } else if (std::abs(ref_point.z - b.point.z) < 1e-6) {
                norm = {0, 0, -1};
            } else if (std::abs(ref_point.z - b.point.z - b.c) < 1e-6) {
                norm = {0, 0, 1};
            }
            return {ref_point, b.color() * Vec3::dot(norm, (light_source - ref_point).normalized())};
        }
        auto *t_p = dynamic_cast<const Tetrahedron *>(fig);
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
};

struct Image {
    int height, width;
    sf::VertexArray matrix;

    Image(int height, int width) : height(height), width(width), matrix({}) {
    }

    sf::Vertex &operator[](size_t i, size_t j) {
        return matrix[i + j * width];
    }
};


int main(int argc, char **argv) {
    char *render_space_file = argv[1];
    char *figures_file = argv[2];
    Vec3 camera{};
    Vec3 norm{};
    Vec3 up{};
    double a0;
    double a1;
    double alpha;
    int width;
    int height;
    Vec3 light_source{};
    {
        std::ifstream in(render_space_file);
        in >> camera >> norm >> up >> a0 >> a1 >> alpha >> width >> height >> light_source;
        norm = norm.normalized();
        up = up.normalized();
        alpha = alpha / 180 * M_PI;
    }
    std::vector<Figure *> figures;
    {
        std::string figure_type;
        std::ifstream in(figures_file);
        while (in >> figure_type) {
            if (figure_type == "sphere") {
                Color color{};
                Point center;
                double radius;
                in >> center >> radius >> color;
                figures.push_back(new Sphere(color, center, radius));
            } else if (figure_type == "box") {
                Color color{};
                Point p;
                double a, b, c;
                in >> p >> a >> b >> c >> color;
                figures.push_back(new Box(color, p, a, b, c));
            } else if (figure_type == "tetra") {
                Color color{};
                Point a, b, c, d;
                in >> a >> b >> c >> d >> color;
                figures.push_back(new Tetrahedron(color, a, b, c, d));
            } else {
                throw "incorrect figure type";
            }
        }
    }
    Image image{height, width};
    double pixel_size = a0 * tan(alpha / 2) * 2 / height;
    double actual_height = height * pixel_size;
    double actual_width = width * pixel_size;
    Point start_point = camera + norm * a0;
    Vec3 tmp = Vec3::cross(norm, up).normalized();
    up = -up;
    start_point = start_point - up * actual_height / 2 - tmp * actual_width / 2;
#pragma opm parallel for
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            Point current_point = start_point + Vec3(pixel_size / 2) + tmp * pixel_size * i + up * pixel_size * j;
            Ray r{camera, (current_point - camera).normalized()};
            Point image_point = camera;
            Color pixel_color{255, 255, 255};
            for (auto &elem: figures) {
                auto p = r.intersect(elem, light_source);
                if (p.first == camera || Vec3::distance(camera, p.first) < a0 ||
                    Vec3::distance(camera, p.first) > a0 + a1) {
                    continue;
                }
                {
                    if (image_point == camera ||
                        Vec3::distance(image_point, camera) > Vec3::distance(p.first, camera)) {
                        image_point = p.first;
                        pixel_color = p.second;
                    }
                }
            }
            sf::Vertex v{};
            v.color = sf::Color(pixel_color.r, pixel_color.g, pixel_color.b);
            v.position = {(float) i, (float) j};
            image.matrix.append(v);
        }
    }
    sf::RenderWindow window(sf::VideoMode(width, height), "Rendered picture");
    while (window.isOpen()) {
        sf::Event event{};
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }
        window.draw(image.matrix);
        window.display();
        window.capture().saveToFile(argv[3]);
        sf::sleep(sf::seconds(1));
    }
    for (auto &elem: figures) {
        delete elem;
    }
    return 0;
}