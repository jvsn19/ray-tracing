//
// Created by unkwis on 03/12/16.
//

#ifndef RT_HEADER_H
#define RT_HEADER_H
#define INF 1e10
#include <algorithm>
#include <utility>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cmath>

using namespace std;

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
        ret.z = z*p;
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

    T3 normal() { //Normalizacao do vetor
        T3 ret;
        double n = this->norma();
        ret.x = x/n;
        ret.y = y/n;
        ret.z = z/n;
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
    T3 dir;
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
    double a, b, c, d, e ,f, g, h, j, k, ka, kd, ks, KS, KT, ir;
    int n;
    Color color;

    Object(double a, double b, double c, double d, double e, double f, double g, double h, double j, double k,
           double ka, double kd, double ks, int n, double KS, double KT, double ir, double red, double green, double blue) {
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
        this->KS = KS;
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

Size size;
Ortho ortho;
vector<Object> objects;
vector<Light> lights;
T3 background;
double ambient;

T3  normalize( T3 v ) {
    T3 retorno;
    double denom;            // Temporary denominator

    //  Absolute value of vector's coordinates

    double x = ( v.x > 0.0 ) ? v.x : - v.x;
    double y = ( v.y > 0.0 ) ? v.y : - v.y;
    double z = ( v.z > 0.0 ) ? v.z : - v.z;


    if ( x > y ) {
        if ( x > z ) {
            y = y / x;
            z = z / x;
            denom = 1.0 / ( x * sqrt( 1.0 + y * y + z * z ) );

        } else {                 // z > x > y
            if ( 1.0 + z > 1.0 ) {
                y = y / z;
                x = x / z;
                denom = 1.0 / ( z * sqrt( 1.0 + y * y + x * x ) );
            }
        }

    } else {
        if ( y > z ) {
            z = z / y;
            x = x / y;
            denom = 1.0 / ( y * sqrt( 1.0 + z * z + x * x ) );

        } else {                 // x < y < z
            if ( 1.0 + z > 1.0 ) {
                y = y / z;
                x = x / z;
                denom = 1.0 / ( z * sqrt( 1.0 + y * y + x * x ) );
            }
        }
    }

    if ( 1.0 + x + y + z > 1.0 ) {
        retorno = v*denom;
    }
    return retorno;
}

double intersect( Ray ray, Object *obj )

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
    double  a, b, c, d, e;        // Coefficents of equation of..
    double  f, g, h, j, k;        // ..quadric surface
    double  acoef, bcoef, ccoef;        // Intersection coefficents
    double  dx, dy, dz;            // Direction - origin coordinates
    double  disc;            // Distance to intersection
    double  root;            // Root of distance to intersection
    double  t;                // Distance along ray to intersection
    double  x0, y0, z0;            // Origin coordinates


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

    T3 temp;
    temp.x = ray.dir.x;
    temp.y = ray.dir.y;
    temp.z = ray.dir.z;
    temp = normalize(temp);

//    cout << temp.x << " " << temp.y << " " << temp.z << endl;

    dx = temp.x;
    dy = temp.y;
    dz = temp.z;

    //dx = ray.dir.x/* - ray.org.x*/;
    //dy = ray.dir.y/* - ray.org.y*/;
    //dz = ray.dir.z/* - ray.org.z*/;

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
}                    // End procedure intersect

bool isShadow(Ray ray) {
    bool retorno = false;
    double dist = INF;
    double distTemp = -1;

    for(int i = 0; i < objects.size(); i++) {
        distTemp = intersect(ray, &objects[i]);
//        cout << distTemp << endl;
        if(distTemp > 0 && distTemp < dist) {
            dist = distTemp;
        }
    }
    if(dist > 0 && dist < INF) retorno = true;

    return retorno;
}

T3 getDir(int i, int j) {
    T3 ret;
    ret.z = 0;
    double pw, ph;
    //Calculando o comprimento e a largura de cada pixel.
    pw = (ortho.x1 > ortho.x0) ? (ortho.x1 - ortho.x0)/ size.w : (ortho.x0 - ortho.x1)/ size.w;
    ph = (ortho.y1 > ortho.y0) ? (ortho.y1 - ortho.y0)/ size.h : (ortho.y0 - ortho.y1)/ size.h;

    //Calculando as coordenadas do pixel
    ret.x = (ortho.x0 + pw/2) + (pw*j);
    ret.y = (ortho.y1 - ph/2) - (ph*i);
    return ret;
}

//Retorna o indice do objeto mais próximo ou -1 caso nao exista interseccao
int nextObject(Ray ray, vector<Object> objects){
    int ret = -1;
    double dist = INF, distAux = -1;
    for(int i = 0; i < objects.size(); ++i){
        distAux = intersect(ray, &objects[i]);
        if(distAux > 0 && distAux < dist) {
            dist = distAux;
            ret = i;
        }
    }
    if (dist == INF) return -1;
//    cout << ret << endl;
    return ret;
}

T3 intersectionPoint(Ray ray, Object object) {
    T3 ret;
    ret.x = ray.dir.x;
    ret.y = ray.dir.y;
    ret.z = ray.dir.z;
    ret = normalize(ret);
    double distance = intersect(ray, &object);
    ret = ret * distance;
    ret = ret + ray.org;
    return ret;
}

T3 normalQuadric(Object object, T3 point) {
    T3 ret;
    double x = point.x;
    double y = point.y;
    double z = point.z;
    double a = object.a;
    double b = object.b;
    double c = object.c;
    double d = object.d;
    double e = object.e;
    double f = object.f;
    double g = object.g;
    double h = object.h;
    double j = object.j;
    ret.x = (2*a*x) + (2*d*y) + (2*f*z) + (2*g);
    ret.y = (2*b*y) + (2*d*x) + (2*e*z) + (2*h);
    ret.z = (2*c*z) + (2*e*y) + (2*f*x) + (2*j);

    return ret;
}

//r = i - 2*n(<n,i>)
double calcSpecularIntensity(Ray ray, Object object, Light light){
    //Finding the normal vector(n)
    T3 intersectPoint = intersectionPoint(ray, object);
    T3 normalObject = normalQuadric(object, intersectPoint);
    normalObject = normalize(normalObject);

    //Vector pointing to the light font(i)
    T3 lightVector = light.dir;
    lightVector = lightVector*(-1);
    lightVector = normalize(lightVector);

    double cosR = lightVector%normalObject;
    if(cosR < 0) cosR = 0;
    T3 r = normalObject*(2*cosR) - lightVector;
    r = normalize(r);

    Ray lightRay = Ray(intersectPoint, lightVector, 0);
    if(!isShadow(lightRay)){
        T3 v = ray.dir;
        v = v*(-1);
        v = normalize(v);
        double cos = r%v;
        if(cos < 0) cos = 0;
        return light.intensity * object.ks * pow(cos, object.n);
    }
    return 0;
}

double calcDifusalIntensity(Ray ray, Object object, Light light) {
    //Finding the normal vector normalized
    T3 point = intersectionPoint(ray, object);
    T3 normalVec = normalQuadric(object, point);
    normalVec = normalize(normalVec);

    //Finding the light vector, pointing to the light font, normalized
    T3 aux = light.dir;
    aux = aux * (-1);
    aux = normalize(aux);

    //Test if this light ray intersect some object. If yes, there will be shadow
    Ray lightRay = Ray(point, aux, 0);
    if(!isShadow(lightRay)){
        double cos = aux%normalVec;
        if(cos < 0) {
            cos = 0;
        };
        //Ip * kd * cos
        return light.intensity * object.kd * cos;
    }
    return 0;
}

T3 shade (Object object, Ray ray){
    T3 ret;
    double intensidadeAmbiente = object.ka*ambient, intensidadeDifusa = 0.0, intensidadeEspecular = 0.0;
    for(int i = 0; i < lights.size(); ++i){
        intensidadeEspecular += calcSpecularIntensity(ray, object, lights[i]);
        intensidadeDifusa += calcDifusalIntensity(ray, object, lights[i]);
    }
    ret.x = object.color.r * (intensidadeAmbiente + intensidadeDifusa) + intensidadeEspecular;
    ret.y = object.color.g * (intensidadeAmbiente + intensidadeDifusa) + intensidadeEspecular;
    ret.z = object.color.b * (intensidadeAmbiente + intensidadeDifusa) + intensidadeEspecular;

    ret.x = min(ret.x, 1.0);
    ret.y = min(ret.y, 1.0);
    ret.z = min(ret.z, 1.0);

    return ret;
}

T3 merge(Object object, T3 color) {
    T3 ret;
    ret.x = ((1 - object.KS - object.KT) * color.x);
    ret.y = ((1 - object.KS - object.KT) * color.y);
    ret.z = ((1 - object.KS - object.KT) * color.z);
    return ret;
}

T3 rayTracing(Ray &ray) {
    T3 ret;
    int objectIndex = nextObject(ray, objects);
    if (objectIndex < 0) {
        ret.x = background.x;
        ret.y = background.y;
        ret.z = background.z;
        return ret;
    }
    Object object = objects[objectIndex];
    T3 color;
    color = shade(object, ray);
//    cout << color.x << " " << color.y << " " << color.z << endl;
    ret = merge(object, color);

    return ret;
}


#include "SDL.cpp"
#endif //RT_HEADER_H
