//
// Created by unkwis on 03/12/16.
//

#ifndef RT_HEADER_H
#define RT_HEADER_H

#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>

typedef struct {
    double x, y, z;
} T3;

typedef struct {
    double r, g, b;
} Color;

typedef struct {
    double x0, y0, x1, y1;
} Ortho;

typedef struct {
    int w, h;
} Size;

typedef struct light{
    T3* coords;
    double intensity;

    light(){
        this->coords = new T3();
    }
    ~light(){
        delete this->coords;
    }
} Light;

typedef struct object{
    double a, b, c, d, e ,f, g, h, j, k, ka, kd, ks, n, KS, KT, ir;
    Color* color;

    object(){
        this->color = new Color();
    }
    ~object(){
        delete this->color;
    }
} Object;

#include "SDL.cpp"
#endif //RT_HEADER_H
