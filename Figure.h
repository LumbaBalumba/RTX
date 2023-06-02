//
// Created by i3alumba on 02.06.23.
//

#ifndef RTX_FIGURE_H
#define RTX_FIGURE_H

#include "Color.h"

struct Figure {
    Color c;

    explicit Figure(Color c) : c(c) {}

    virtual ~Figure() = default;

    [[nodiscard]] Color color() const;
};

#endif //RTX_FIGURE_H
