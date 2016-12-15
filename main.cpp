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
    supersampling = sdl.getSuperSampling();
    depth = sdl.getDepth();
    ph = fabs(ortho.x1 - ortho.x0) / size.w;    //Pixel's weight
    pw = fabs(ortho.y1 - ortho.y0) / size.h;    //Pixel's height

    for(int i = 0; i < lights.size(); ++i){
        lights[i].dir = lights[i].dir*(-1);
    }

    cout << "P3" << endl;
    cout << size.w << " " << size.h << endl;
    cout << 255 << endl;
    for(int i = 0; i < size.h; ++i){
        for(int j = 0; j < size.w; ++j){
            T3 direction = (getDir(i, j) - sdl.getEye());  //Direction of the ray vector
            Ray ray = Ray(sdl.getEye(), direction, 0);
            T3 cor = rayTracing(ray);
           /* //Melhorar
            if(supersampling) {
                Ray extra = ray;
                extra.dir.x -= pw/3.0;
                extra.dir.y -= ph/3.0;
                cor = cor + rayTracing(extra);
                extra.dir.x += pw/3.0;
                extra.dir.y -= ph/3.0;
                cor = cor + rayTracing(extra);
                extra.dir.x -= pw/3.0;
                extra.dir.y += ph/3.0;
                cor = cor + rayTracing(extra);
                extra.dir.x += pw/3.0;
                extra.dir.y += ph/3.0;
                cor = cor + rayTracing(extra);
                cor = cor*(1.0/4.0);
            }*/
            cout << cor.x*255 << " " << cor.y*255 << " " << cor.z *255<< endl;
        }
    }
    return 0;
}
