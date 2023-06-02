//
// Created by i3alumba on 02.06.23.
//

#include "Image.h"

sf::Vertex &Image::operator[](size_t i, size_t j) {
    return matrix[i + j * width];
}
