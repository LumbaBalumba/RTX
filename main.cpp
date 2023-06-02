#include <iostream>
#include <fstream>
#include <vector>
#include <SFML/Graphics.hpp>

#include "Vec3.h"
#include "Color.h"
#include "Figure.h"
#include "Sphere.h"
#include "Box.h"
#include "Tetra.h"
#include "Ray.h"
#include "Image.h"

/*
Аргументы запуска программы - файл с описанием 3D-пространства, файл с описанием фигур, название выходного файла.
Ввод в первом файле осуществляется так же, как описано в ТЗ, за исключением того, что не требуется вводназвания параметров.
Во втором файле вводятся фигуры по следующему формату:
    1) первым вводится вид фигуры
    2) затем вводятся параметры формы (для сферы - центр и радиус, для параллелепипеда - точка и три стороны, для тетраэдра - четыре точки)
    3) последним вводится цвет фигуры, три числа, каждое от 0 до 255
 */


int main(int argc, char **argv) {
    if (argc != 4) {
        throw "Incorrect arguments";
    }
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
        if (!in.is_open()) {
            throw "Could not open space description file";
        }
        in >> camera >> norm >> up >> a0 >> a1 >> alpha >> width >> height >> light_source;
        norm = norm.normalized();
        up = up.normalized();
        alpha = alpha / 180 * M_PI;
    }
    std::vector<Figure *> figures;
    {
        std::string figure_type;
        std::ifstream in(figures_file);
        if (!in.is_open()) {
            throw "Could not find open file";
        }
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
                figures.push_back(new Tetra(color, a, b, c, d));
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
                if (Vec3::distance(camera, p.first) < a0 || Vec3::distance(camera, p.first) > a0 + a1) {
                    continue;
                }
                if (image_point == camera ||
                    Vec3::distance(image_point, camera) > Vec3::distance(p.first, camera)) {
                    image_point = p.first;
                    pixel_color = p.second;
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