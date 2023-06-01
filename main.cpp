#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <SFML/Graphics.hpp>
#include <limits>
#include <omp.h>

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
        return {std::max(left.r * right.r, 0.0), std::max(left.g * right.g, 0.0), std::max(left.b * right.b, 0.0)};
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
            double tmin = std::numeric_limits<double>::min();
            double tmax = std::numeric_limits<double>::max();

            if (direction.x != 0) {
                double tx1 = (b.point.x - point.x) / direction.x;
                double tx2 = (b.point.x + b.a - point.x) / direction.x;
                tmin = std::max(tmin, std::min(tx1, tx2));
                tmax = std::min(tmax, std::max(tx1, tx2));
            }

            if (direction.y != 0) {
                double ty1 = (b.point.y - point.y) / direction.y;
                double ty2 = (b.point.y + b.b - point.y) / direction.y;
                tmin = std::max(tmin, std::min(ty1, ty2));
                tmax = std::min(tmax, std::max(ty1, ty2));
            }

            if (direction.z != 0) {
                double tz1 = (b.point.z - point.z) / direction.z;
                double tz2 = (b.point.z + b.c - point.z) / direction.z;
                tmin = std::max(tmin, std::min(tz1, tz2));
                tmax = std::min(tmax, std::max(tz1, tz2));
            }

            if (tmax >= tmin && tmax >= 0) {
                Point intersection_point = {
                        point.x + tmin * direction.x,
                        point.y + tmin * direction.y,
                        point.z + tmin * direction.z
                };
                Vec3 normal{};
                if (std::abs(intersection_point.x - b.point.x) < 1e-6) {
                    normal = {-1, 0, 0};
                } else if (std::abs(intersection_point.x - b.point.x - b.a) < 1e-6) {
                    normal = {1, 0, 0};
                } else if (std::abs(intersection_point.y - b.point.y) < 1e-6) {
                    normal = {0, -1, 0};
                } else if (std::abs(intersection_point.y - b.point.y - b.b) < 1e-6) {
                    normal = {0, 1, 0};
                } else if (std::abs(intersection_point.z - b.point.z) < 1e-6) {
                    normal = {0, 0, -1};
                } else if (std::abs(intersection_point.z - b.point.z - b.c) < 1e-6) {
                    normal = {0, 0, 1};
                }
                return {intersection_point,
                        b.color() * Vec3::dot(normal, (light_source - intersection_point).normalized())};
            }
            return {point, {255, 255, 255}};
        }
        auto *t_p = dynamic_cast<const Tetrahedron *>(fig);
        if (t_p != nullptr) {
            auto tetrahedron = *t_p;
            Point intersection_point;
            double t = -1;
            std::array<Point, 4> points = {tetrahedron.a, tetrahedron.b, tetrahedron.c, tetrahedron.d};
            Vec3 n{};
            for (int i = 0; i < 4; i++) {
                Point p1 = points[i];
                Point p2 = points[(i + 1) % 4];
                Point p3 = points[(i + 2) % 4];
                Point normal = Vec3::cross({p2.x - p1.x, p2.y - p1.y, p2.z - p1.z},
                                           {p3.x - p1.x, p3.y - p1.y, p3.z - p1.z});
                double d = Vec3::dot(normal, {p1.x - point.x, p1.y - point.y, p1.z - point.z}) /
                           Vec3::dot(normal, direction);
                if (d > 0 && (t < 0 || d < t)) {
                    Point point_on_plane = {point.x + d * direction.x, point.y + d * direction.y,
                                            point.z + d * direction.z};
                    double area =
                            Vec3::distance(p1, p2) * Vec3::distance(p2, p3) * Vec3::distance(p3, p1) /
                            (4 * Vec3::distance(normal, {0, 0, 0}));
                    double area1 = Vec3::distance(point_on_plane, p1) * Vec3::distance(p1, p2) *
                                   Vec3::distance(p2, point_on_plane) /
                                   (4 * Vec3::distance(normal, {0, 0, 0}));
                    double area2 = Vec3::distance(point_on_plane, p2) * Vec3::distance(p2, p3) *
                                   Vec3::distance(p3, point_on_plane) /
                                   (4 * Vec3::distance(normal, {0, 0, 0}));
                    double area3 = Vec3::distance(point_on_plane, p3) * Vec3::distance(p3, p1) *
                                   Vec3::distance(p1, point_on_plane) /
                                   (4 * Vec3::distance(normal, {0, 0, 0}));
                    if (std::abs(area - area1 - area2 - area3) <= 1e-6) {
                        intersection_point = point_on_plane;
                        t = d;
                    }
                    n = normal;
                }
            }
            if (intersection_point == point) {
                return {point, {255, 255, 255}};
            }
            return {intersection_point,
                    tetrahedron.color() * Vec3::dot((light_source - intersection_point).normalized(), n.normalized())};
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