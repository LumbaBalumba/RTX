#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <SFML/Graphics.hpp>

enum {
    DEFAULT_PICTURE_WIDTH = 800,
    DEFAULT_PICTURE_HEIGHT = 600
};

struct Vec3 {
    double x;
    double y;
    double z;

    Vec3() = default;

    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    explicit Vec3(double value) : x(value), y(value), z(value) {}

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
        return {std::max(left.x * right.x, 0.0), std::max(left.y * right.y, 0.0), std::max(left.z * right.z, 0.0)};
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

    friend std::istream &operator>>(std::istream &in, Vec3 v) {
        return in >> v.x >> v.y >> v.z;
    }

    friend std::ostream &operator<<(std::ostream &out, const Vec3 &v) {
        return out << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    }

    [[nodiscard]] double length() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    friend double distance(Vec3 left, Vec3 right) {
        return (left - right).length();
    }
};

typedef Vec3 Point;

struct Color {
    double r;
    double g;
    double b;

    Color() = default;

    Color(double x, double y, double z) : r(x), g(y), b(z) {}

    explicit Color(double value) : r(value), g(value), b(value) {}


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

    friend std::istream &operator>>(std::istream &in, Color v) {
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
};

struct Sphere : public Figure {
    Vec3 center;
    double radius;

    Sphere(Color c, Vec3 center, double radius) : Figure(c), center(center), radius(radius) {}
};

struct Box : public Figure {
    Vec3 point;
    double a;
    double b;
    double c;

    Box(Color color, Vec3 point, double a, double b, double c) : Figure(color), point(point), a(a), b(b), c(c) {}
};

struct Tetrahedron : public Figure {
    Vec3 a;
    Vec3 b;
    Vec3 c;
    Vec3 d;

    Tetrahedron(Color color, Vec3 a, Vec3 b, Vec3 c, Vec3 d) : Figure(color), a(a), b(b), c(c), d(d) {}

    [[nodiscard]] bool contains(Vec3 p) const {
        Vec3 v0 = {a.x - d.x, a.y - d.y, a.z - d.z};
        Vec3 v1 = {b.x - d.x, b.y - d.y, b.z - d.z};
        Vec3 v2 = {c.x - d.x, c.y - d.y, c.z - d.z};
        Vec3 n = Vec3::cross(v1, v0);
        double dot = Vec3::dot(n, v2);
        if (dot == 0) {
            return false;
        }
        double inv_dot = 1.0 / dot;
        Vec3 w = {p.x - d.x, p.y - d.y, p.z - d.z};
        double u = Vec3::dot(Vec3::cross(v1, v2), w) * inv_dot;
        double v = Vec3::dot(Vec3::cross(v0, v2), w) * inv_dot;
        double t1 = Vec3::dot(Vec3::cross(v0, v1), w) * inv_dot;
        if (u < 0 || v < 0 || t1 < 0 || u + v + t1 > 1) {
            return false;
        }
        return true;
    }
};


struct Ray {
    Vec3 point;
    Vec3 direction;

    Ray(Vec3 point, Vec3 direction) : point(point), direction(direction) {}

    Vec3 intersect(const Figure *fig) const {
        auto *s_p = dynamic_cast<const Sphere *>(fig);
        if (s_p != nullptr) {
            /* sphere intersection */
        }
        auto *b_p = dynamic_cast<const Box *>(fig);
        if (b_p != nullptr) {
            /* box intersection */
        }
        auto *t_p = dynamic_cast<const Tetrahedron *>(fig);
        if (t_p != nullptr) {
            /* tetra intersection */
        }
        throw "Invalid figure type";
    }
};

struct Image {
    size_t height, width;
    sf::VertexArray matrix;

    Image(size_t height, size_t width) : height(height), width(width) {
        matrix.resize(height * width);
    }

    sf::Vertex &operator[](size_t i, size_t j) {
        return matrix[i * width + j];
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
    size_t width;
    size_t height;
    Vec3 light_source{};
    {
        std::ifstream in(render_space_file);
        in >> camera >> norm >> up >> a0 >> a1 >> alpha >> width >> height >> light_source;
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
    sf::RenderWindow window(sf::VideoMode(DEFAULT_PICTURE_WIDTH, DEFAULT_PICTURE_HEIGHT), "Rendered picture");
    /* main stuff */
    Image image{height, width};
    while (window.isOpen()) {
        sf::Event event{};
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }
        window.draw(image.matrix);
        window.display();
    }
    /* save window to image*/
    return 0;
}