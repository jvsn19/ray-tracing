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
    Color(){};
    Color(double r, double g, double b) {
        this->r = r;
        this->g = g;
        this->b = b;
    }
    Color operator *(double val){
        Color color = Color(r*val, g*val, b*val);
        return color;
    }
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
           double ka, double kd, double ks, int n, double KS, double KT, double ir, Color color) {
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
        this->color = color;
    }
};


struct Ray {
    T3 org, dir;
    int depth;
    Ray();
    Ray(T3 org, T3 dir, int depth){
        this->org = org;
        this->dir = dir;
        this->depth = depth;
    }
};

struct Texture {
    int wres, hres, maxGrey;
    string magicValue;
    vector< vector<Color> > texMatrix;

    Texture();
    Texture (string &filePath) {
        int lineCounter = 1;
        ifstream ifs;
        ifs.open(filePath, std::ifstream::in);
        while(ifs.eof()) {
            if(ifs.peek() == '\n' || ifs.peek() == ' ') {
                ifs.get();
            }
            else if(ifs.peek() == '#'){
                ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            else {
                if(ifs.peek() == 'P'){
                    ++lineCounter;
                    ifs >> magicValue;
                }
                else if(lineCounter == 2) {
                    ifs >> wres >> hres;
                    ++lineCounter;
                }
                else {
                    for(int i = 0; i < hres; ++i) {
                        for(int j = 0; j < wres; ++j) {
                            double r, g, b;
                            ifs >> r >> g >> b;
                            Color color = Color(r,g,b);
                            texMatrix[i][j] = color;
                        }
                    }
                }
            }
        }
    }

};

Size size;
Ortho ortho;
vector<Object> objects;
vector<Light> lights;
T3 background;
double ambient;
bool supersampling;
int depth;
double pw;
double ph;

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

    dx = temp.x;
    dy = temp.y;
    dz = temp.z;

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
}




bool isShadow(Ray ray) {
    for(int i = 0; i < objects.size(); i++) {
        if (intersect(ray, &objects[i]) != -1) return true;
    }
    return false;
}

//Metodo complementar para formar o vetor diretor do raio.
T3 getDir(int i, int j) {
    T3 ret;
    ret.z = 0;  //A janela esta no plano z = 0
    /*pw e ph calculados anteriormente.
     *pw = pixel weight. Calculado dividindo o comprimento total da janela fabs(x1 - x0) pela resolucao
     *ph = pixel height. Mesmo calculo, porem utilizando y1 - y0
    */
    ret.x = (ortho.x0 + pw/2) + (pw*j);
    ret.y = (ortho.y1 - ph/2) - (ph*i);
    return ret;
}

//Retorna o indice do objeto mais próximo ou -1 caso nao exista interseccao
int nextObject(Ray ray, vector<Object> objects){
    int ret = -1;
    double dist = INF, distAux = -1;
    //Simplesmente um for que percorre todos os objetos do meio
    for(int i = 0; i < objects.size(); ++i){
        distAux = intersect(ray, &objects[i]);
        if(distAux > 0 && distAux < dist) { //Caso eu ache algum objeto, testo se ele eh mais proximo que o meu atual
            //Caso seja, substituo
            dist = distAux;
            ret = i;
        }
    }
    if (dist == INF) return -1; //Caso a distancia ao obj mais proximo for INF, nenhum obj foi encontrado
    return ret;
}

//Metodo para calcular o ponto de interseccao
T3 intersectionPoint(Ray ray, Object object) {
    /*
     * O metodo consiste em, a partir do vetor diretor, encontrar a intersecao multiplicando o mesmo pela
     * distancia ate o objeto. Dessa forma eu sei que esse vetor tocara o objeto
     * Utilizando a equacao P = Q + tr
     * IntersectionPoint = Origin + t*(ray.dir)
     */
    T3 ret = ray.dir;
    ret = normalize(ret);
    double distance = intersect(ray, &object);
    ret = ret * distance;
    ret = ret + ray.org;
    return ret;
}

//Metodo para retornar a normal de uma quadrica em um determinado ponto
T3 normalQuadric(Object object, T3 point) {
    /*
     * f(x,y,z) = ax^2 + by^2 + cz^2 + 2dxy + 2eyz + 2fxz + 2gx + 2hy + 2jz + k = 0 forma quadratica
     * Equacao da normal encontrada no site http://euclid.nmu.edu/~bpeterso/CS446-Handouts/Notes/CS446_Note_7.pdf
     *
     */
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
    ret.x = (2*a*x) + (2*e*z) + (2*f*y) + (2*g);
    ret.y = (2*b*y) + (2*d*z) + (2*f*x) + (2*h);
    ret.z = (2*c*z) + (2*d*y) + (2*e*x) + (2*j);

    return ret;
}

// Calcular o raio refletido para realizar a recursao do raytracing
Ray reflectionRay(Object object, Ray ray){
    /*
     * 1) Ponto de interseccao entre o raio e o objeto
     * 2) No ponto encontrado deveremos encontrar o vetor normal (N)
     * 3) Calculamos o vetor que atingiu o objeto (V). Como o sentido de V eh ray.org -> intersectPoint, devemos
     *    inverter o seu sentido: V = - V
     * 4) Tendo em posse o vetor normal e o vetor de incidencia de luz, calculamos o vetor refletido a partir da
     *    equacao: R = 2*<N,V>*N - V
     *    OBS: Quando temos <N,L> < 0, significa que a normal esta dentro do objeto. Logo N = -N
     * 5) Tendo entao o ponto de origem e o vetor de direcao do raio, podemos fazer o nosso Ray rayReflected, aumen-
     *    Tando em 1 a sua profundidade.
     */
    T3 intersectPoint = intersectionPoint(ray, object);
    T3 N = normalQuadric(object, intersectPoint);
    N = normalize(N);

    T3 V = ray.dir;
    V = V*(-1);
    V = normalize(V);

    double cosNV = N%V;
    if(cosNV < 0) N = N*(-1);
    T3 R = N*(2*cosNV)-V;
    Ray rayReflected = Ray(intersectPoint, R, ray.depth+1);

    return rayReflected;
}

//Metodo para calcular a componente especular da equacao de Phong
double calcSpecularIntensity(Ray ray, Object object, Light light){
    /*
     * Is = ks*Il*(<R,V>)^n
     * 1) Calculamos o ponto de interseccao entre o raio e o objeto
     * 2) Vejo se existe algum objeto entre este mesmo e o foco de luz. Caso exista, nao deveremos colocar nenhuma
     *    cor
     * 3) Tendo em posse esse ponto, podemos calcular o vetor normal ao objeto naquele local (N)
     * 4) Calculamos entao o vetor que aponta para o foco de luz atual (tendo em vista que uma cena pode conter
     *    varios focos de luz. Como a direcao da luz esta no sentido luz->objeto, devemos inverte-la: L = - L
     * 5) Com o L e N, podemos calcular o vetor refletido pelo objeto usando a equacao
     *    R = 2*<N,L>*N - 1.
     *    OBS: Caso o <N,L> < 0, nao teremos componente especular. Existe outro caso que sera explicitado mais a frente
     * 6) Calculamos o vetor incidente de visao(V), no caso o ray.dir. Da mesma forma que no metodo de reflexao, o sen-
     *    tido desse vetor esta invertido. V = -V
     * 7) De posse do V e R, podemos calcular a componente especular:
     *    Is = ks*Il*(<R,V>)^n
     *    OBS: Caso <R,V> < 0, nao existe componente especular
     */
    T3 intersectPoint = intersectionPoint(ray, object);
    T3 L = light.dir;
    L = L*(-1);
    L = normalize(L);
    Ray shadow = Ray(intersectPoint, L, 0);
    if(isShadow(shadow)) return 0;

    T3 N = normalQuadric(object, intersectPoint);
    N = normalize(N);

    double cosNL = N%L;
    if(cosNL < 0) return 0;
    T3 R = N*(2*cosNL) - L;
    R = normalize(R);

    T3 V = ray.dir;
    V = V*(-1);
    V = normalize(V);

    double cosRV = R%V;
    if(cosRV < 0) return 0;
    double Is = object.ks*light.intensity*pow(cosRV, object.n);
    return Is;

}

double calcDifusalIntensity(Ray ray, Object object, Light light) {
    /*
     * Id = kd*Il*<N,L>
     * 1) Calculando o ponto de interseccao entre o raio e o objeto
     * 2) Vemos se a luz consegue alcancar nessa superficie. Caso exista algum objeto entre a luz e este mesmo,
     *    nao deve haver cor aqui
     * 3) Calculando o vetor normal ao objeto no ponto especificado (N)
     * 4) Calculamos o vetor que aponta para a luz (L)
     * 5) Aplicamos a equacao acima
     *    OBS: se <N,L> < 0, nao teremos componente difusa
     */
    T3 intersectPoint = intersectionPoint(ray, object);
    T3 L = light.dir;
    L = L*(-1);
    L = normalize(L);
    Ray shadow = Ray(intersectPoint, L, 0);
    if(isShadow(shadow)) return 0;

    T3 N = normalQuadric(object, intersectPoint);
    N = normalize(N);

    double cosNL = N%L;
    if(cosNL < 0) return 0;

    double Id = object.kd*light.intensity*cosNL;
    return Id;
}

Color shading (Object object, Ray ray){
    Color ret;
    double iAmbient = object.ka*ambient, iDifusal = 0.0, iSpecular = 0.0;
    for(int i = 0; i < lights.size(); ++i){
        iSpecular += calcSpecularIntensity(ray, object, lights[i]);
        iDifusal += calcDifusalIntensity(ray, object, lights[i]);
    }
    ret.r = object.color.r * (iAmbient + iDifusal) + iSpecular;
    ret.g = object.color.g * (iAmbient + iDifusal) + iSpecular;
    ret.b = object.color.b * (iAmbient + iDifusal) + iSpecular;

    ret.r = min(ret.r, 1.0);
    ret.g = min(ret.g, 1.0);
    ret.b = min(ret.b, 1.0);

    return ret;
}

Color merge(Object object, Color color, Color refColor) {
    Color ret;
    ret.r = ((1 - object.KS - object.KT) * color.r) + (object.KS*refColor.r);
    ret.g = ((1 - object.KS - object.KT) * color.g) + (object.KS*refColor.g);
    ret.b = ((1 - object.KS - object.KT) * color.b) + (object.KS*refColor.b);
    return ret;
}

void normalizeColor(Color *color) {
    color->r = min(color->r, 1.0);
    color->g = min(color->g, 1.0);
    color->b = min(color->b, 1.0);
}

Color rayTracing(Ray &ray) {
    Color ret;
    int objectIndex = nextObject(ray, objects);
    if (objectIndex < 0) {
        if(ray.depth == 0){
            ret.r = background.x;
            ret.g = background.y;
            ret.b = background.z;
            return ret;

        }
        ret.r = 0.0;
        ret.g = 0.0;
        ret.b = 0.0;
        return ret;
    }
    Object object = objects[objectIndex];
    Color color;
    color = shading(object, ray);

    Color refColor = {0.0, 0.0, 0.0};

    if(ray.depth < depth) {
        if(object.KS > 0) {
            Ray rayRef = reflectionRay(object, ray);
            refColor=rayTracing(rayRef);
        }
    }

    ret = merge(object, color, refColor);

    normalizeColor(&ret);

    return ret;
}

void printImage(const vector<vector<Color>> &objMatrix) {
    cout << "P3" << endl;
    cout << size.w << " " << size.h << endl;
    cout << 255 << endl;
    for(int i = 0; i < size.h; ++i){
        for(int j = 0; j < size.w; ++j){
            cout << objMatrix[i][j].r << " " << objMatrix[i][j].g << " " << objMatrix[i][j].b << endl;
        }
    }
}


#include "SDL.cpp"
#endif //RT_HEADER_H
