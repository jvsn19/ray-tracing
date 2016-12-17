#include "Support/Header.h"

void printImage(const vector<vector<Color>> &objMatrix);

void getPNMTexture(int objectIndex);

using namespace std;

//O vetor diretor da Camera sempre é (0,0,1) pois esse vetor deve ser paralelo ao vetor normal ao grid.

int main(void){
    SDL sdl = SDL("Files/onesphere.sdl");

    //Recuperando os elementos da cena
    size = sdl.getSize();
    ortho = sdl.getOrtho();
    objects = sdl.getObjects();
    lights = sdl.getLights();
    background = sdl.getBackground();
    ambient = sdl.getAmbient();
    supersampling = sdl.getSuperSampling();
    depth = sdl.getDepth();
    ph = fabs(ortho.x1 - ortho.x0) / size.w;    //Pixel's weight
    pw = fabs(ortho.y1 - ortho.y0) / size.h;    //Pixel's height

    //Formalizando a direcao da luz
    for(int i = 0; i < lights.size(); ++i){
        lights[i].dir = lights[i].dir*(-1);
    }

    vector< vector<Color> > objMatrix((unsigned long) size.h);
    for(int i = 0; i < size.h; ++i) {
        objMatrix[i] = vector<Color>((unsigned long) size.w);
    }

    for(int i = 0; i < size.h; ++i){
        for(int j = 0; j < size.w; ++j){
            T3 direction = (getDir(i, j) - sdl.getEye());
            Ray ray = Ray(sdl.getEye(), direction, 0);
            objMatrix[i][j] = rayTracing(ray)*255;
        }
    }
    printImage(objMatrix);
    //getPNMTexture(0);

    return 0;
}



