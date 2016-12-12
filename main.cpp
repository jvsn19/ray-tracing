#include "Support/Header.h"

using namespace std;

//O vetor diretor da Camera sempre é (0,0,1) pois esse vetor deve ser paralelo ao vetor normal ao grid.

int main(void){
    SDL sdl = SDL("Files/onesphere.sdl");

    //Recuperando os elementos da cena
    Size size = sdl.getSize();
    Ortho ortho = sdl.getOrtho();
    vector<Object> objects = sdl.getObjects();
    vector<Light> lights = sdl.getLights();

    cout << "P3" << endl;
    cout << size.w << " " << size.h << endl;
    cout << 255 << endl;

    for(int i = 0; i < size.h; ++i){
        for(int j = 0; j < size.w; ++j){
            Ray ray = Ray(sdl.getEye(), getDir(ortho, size, i, j) - sdl.getEye(), 0);
            for(int k = 0; k < objects.size(); ++k){
                if(intersect(ray, &objects[k]) != -1) cout << "200" << endl;
                else cout << 0 << endl;
            }
        }
    }
    return 0;
}
