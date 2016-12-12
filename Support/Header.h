//
// Created by unkwis on 03/12/16.
//

#ifndef RT_HEADER_H
#define RT_HEADER_H
#define INF 1e10
#include <utility>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cmath>

struct T3 {
    double x, y, z;

    T3 operator+ (T3 t) { //Soma de vetores
        T3 ret;
        ret.x = x + t.x;
        ret.y = y + t.y;
        ret.z = z + t.z;
        return ret;
    }

    T3 operator* (double p) { //Produto por um escalar
        T3 ret;
        ret.x = x*p;
        ret.y = y*p;
        ret.z = x*p;
        return ret;
    }

    T3 operator- (T3 t) { //Subtracao de vetores
        T3 ret;
        ret.x = x - t.x;
        ret.y = y - t.y;
        ret.z = z - t.z;
        return ret;
    }

    double operator% (T3 t) { //Produto escalar
        return x*t.x + y*t.y + z*t.z;
    }

    T3 operator^ (T3 t) { //Produto vetorial
        T3 ret;
        ret.x = y*t.z - z*t.y;
        ret.y = z*t.x - x*t.z;
        ret.z = x*t.z - y*t.x;
        return ret;
    }

    double norma() {    //Norma de um vetor
        return sqrt(x*x + y*y + z*z);
    }

    T3 normalizar() { //Normalizacao do vetor
        T3 ret;
        ret.x = x;
        ret.y = y;
        ret.z = z;
        double n = this->norma();
        ret.x /= n;
        ret.y /=n;
        ret.z /=n;
        return ret;
    }
};

struct Color{
    double r, g, b;
};

struct Ortho{
    double x0, y0, x1, y1;
};

struct Size {
    int w, h;
};

struct Light {
    T3 coords;
    double intensity;
};

struct Camera {
    T3 coords, dir, u, v;
    double dist, hx, hy;
    Camera(T3& coods, T3& u, T3& v, double dist, double hx, double hy) {
        dir = T3{0,0,1};    //The ortho plain is fixed on the z plan.
        this->coords = coords;
        this->u = u;
        this->v = v;
        this->dist = dist;
        this->hx = hx;
        this->hy = hy;
    }
};

struct Object{
    double a, b, c, d, e ,f, g, h, j, k, ka, kd, ks, n, KS, KT, ir;
    Color color;

    Object(double a, double b, double c, double d, double e, double f, double g, double h, double j, double k,
           double ka, double kd, double ks, double n, double KS, double KT, double ir, double red, double green, double blue) {
        this->a = a;
        this->b = b;
        this->c = c;
        this->d = d;
        this->e = e;
        this->f = f;
        this->g = g;
        this->h = h;
        this->j = j;
        this->k = k;
        this->ka = ka;
        this->kd = kd;
        this->ks = ks;
        this->n = n;
        this->KS = ks;
        this->KT = KT;
        this->ir = ir;
        color.r = red;
        color.g = green;
        color.b = blue;
    }
};

struct Ray {
    T3 org, dir;
    int depth;
    Ray(T3 org, T3 dir, int depth){
        this->org = org;
        this->dir = dir;
        this->depth = depth;
    }
};

double intersect(Ray ray, Object *obj )

//  Compute the intersection point, if it exists, between the given ray
//  and the given object
//
//  Victor uses a strange formula from Watt & Watt for quadrics:
//
//     Ax^2 + Ey^2 + Hz^2 + Bxy + Fyz + Cxz + Dx + Gy + Jz + K
//
//  rather than the standard formula:
//
//     Ax^2 + By^2 + Cz^2 + 2Dxy + 2Eyz + 2Fxz + 2Gx + 2Hy + 2Jz + K
//
//  ray:  Ray being shot into scene
//  obj:  Object to test for intersection
{
    double  a, b, c, d, e;// Coefficents of equation of..
    double  f, g, h, j, k;// ..quadric surface
    double  acoef, bcoef, ccoef;// Intersection coefficents
    double  dx, dy, dz;// Direction - origin coordinates
    double  disc;// Distance to intersection
    double  root;// Root of distance to intersection
    double  t;// Distance along ray to intersection
    double  x0, y0, z0;// Origin coordinates


    a = obj->a;
    b = obj->b;
    c = obj->c;
    d = obj->d;
    e = obj->e;
    f = obj->f;
    g = obj->g;
    h = obj->h;
    j = obj->j;
    k = obj->k;

    dx = ray.dir.x - ray.org.x;
    dy = ray.dir.y - ray.org.y;
    dz = ray.dir.z - ray.org.z;

    x0 = ray.org.x;
    y0 = ray.org.y;
    z0 = ray.org.z;

    acoef = 2 * f * dx * dz + 2 * e * dy * dz + c * dz * dz + b * dy * dy +
            a * dx * dx + 2 * d * dx * dy;

    bcoef = 2 * b * y0 * dy + 2 * a * x0 * dx + 2 * c * z0 * dz +
            2 * h * dy + 2 * g * dx + 2 * j * dz + 2 * d * x0 * dy +
            2 * e * y0 * dz + 2 * e * z0 * dy + 2 * d * y0 * dx +
            2 * f * x0 * dz + 2 * f * z0 * dx;

    ccoef = a * x0 * x0 + 2 * g * x0 + 2 * f * x0 * z0 + b * y0 * y0 +
            2 * e * y0 * z0 + 2 * d * x0 * y0 + c * z0 * z0 + 2 * h * y0 +
            2 * j * z0 + k;

    //  The following was modified by David J. Brandow to allow for planar
    //  quadrics

    if ( 1.0 + acoef == 1.0 ) {
        if ( 1.0 + bcoef == 1.0 ) {
            return -1.0;
        }

        t = ( -ccoef ) / ( bcoef );

    } else {
        disc = bcoef * bcoef - 4 * acoef * ccoef;
        if ( 1.0 + disc < 1.0 ) {
            return -1.0;
        }

        root = sqrt( disc );
        t = ( -bcoef - root ) / ( acoef + acoef );
        if ( t < 0.0 ) {
            t = ( -bcoef + root ) / ( acoef + acoef );
        }
    }

    if ( t < 0.001 )
        return -1.0;

    return t;
}// End procedure intersect

T3 getDir(Ortho ortho, Size size, int i, int j) {
    T3 ret;
    ret.z = 0;
    double pw, ph;
    //Calculando o comprimento e a largura de cada pixel.
    pw = (ortho.x0 > ortho.x1) ? (ortho.x0 - ortho.x1)/ size.h : (ortho.x1 - ortho.x0)/ size.w;
    ph = (ortho.y0 > ortho.y1) ? (ortho.y0 - ortho.y1)/ size.h : (ortho.y1 - ortho.y0)/ size.h;

    //Calculando as coordenadas do pixel
    ret.x = (ortho.x0 + pw/2) + (pw*j);
    ret.y = (ortho.y0 + ph/2) + (ph*i);
    return ret;
}

double rayTracing() {

}


#include "SDL.cpp"
#endif //RT_HEADER_H
