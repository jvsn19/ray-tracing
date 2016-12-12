#include "Support/Header.h"

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

    cout << "P3" << endl;
    cout << size.w << " " << size.h << endl;
    cout << 255 << endl;

    for(int i = 0; i < size.h; ++i){
        for(int j = 0; j < size.w; ++j){
            Ray ray = Ray(sdl.getEye(), getDir(ortho, size, i, j) - sdl.getEye(), sdl.getDepth());
            T3 cor = rayTracing(ray);
            cout << cor.x*255 << " " << cor.y*255 << " " << cor.z *255<< endl;
        }
    }
    return 0;
}
