//
// Created by i3alumba on 02.06.23.
//

#ifndef RTX_IMAGE_H
#define RTX_IMAGE_H

#include <SFML/Graphics.hpp>

struct Image {
    int height, width;
    sf::VertexArray matrix;

    Image(int height, int width) : height(height), width(width), matrix({}) {}

    sf::Vertex &operator[](size_t i, size_t j);
};


#endif //RTX_IMAGE_H
