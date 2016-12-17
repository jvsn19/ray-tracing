#include "Support/Header.h"

void printImage(const vector<vector<Color>> &objMatrix);

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
        for(int j = 0; j < size.w; ++j) {
            objMatrix[i][j] = Color(0.0, 0.0, 0.0);
        }
    }
    for(int i = 0; i < size.h; ++i){
        for(int j = 0; j < size.w; ++j){
            //Supersampling
            for(double frag_i = i; frag_i < i+1.0; frag_i += 0.5){
                for(double frag_j = j; frag_j <= j + 1.0; frag_j += 0.5){
                    T3 direction = (getDir(frag_i, frag_j) - sdl.getEye());
                    Ray ray = Ray(sdl.getEye(), direction, 0);
                    objMatrix[i][j] = objMatrix[i][j] + (rayTracing(ray))*0.15; //Diminui o coeficiente de cada raio
                }
            }
            normalizeColor(&objMatrix[i][j]);
            objMatrix[i][j] = objMatrix[i][j]*255;
        }
    }
    printImage(objMatrix);
    //getPNMTexture(0);

    return 0;
}



