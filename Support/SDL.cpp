//
// Created by unkwis on 03/12/16.
//

#include "SDL.h"
SDL::SDL(std::string path) {
    std::ifstream ifs;
    ifs.open (path, std::ifstream::in);
    this->supersampling = false;
    while(!ifs.eof()){
        if(ifs.peek() == '\n' || ifs.peek() == ' ') {
            ifs.get();
        }
        else if(ifs.peek() == '#'){
            ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        else {
            std::string tag;
            ifs >> tag;
            if(!tag.compare("output")){
                ifs >> this->output;
            }
            else if(!tag.compare(("eye"))){
                T3 eye = T3();
                ifs >> eye.x
                    >> eye.y
                    >> eye.z;
                this->eye = eye;
            }
            else if(!tag.compare("ortho")){
                Ortho ortho = Ortho();
                ifs >> ortho.x0
                    >> ortho.y0
                    >> ortho.x1
                    >> ortho.y1;
                this->ortho = ortho;
            }
            else if(!tag.compare(("size"))){
                Size size = Size();
                ifs >> size.w
                    >> size.h;
                this->size = size;
            }
            else if(!tag.compare("background")){
                T3 background = T3();
                ifs >> background.x
                    >> background.y
                    >> background.z;
                this->background = background;
            }
            else if(!tag.compare("ambient")){
                ifs >> this->ambient;
            }
            else if(!tag.compare("light")){
                Light light = Light();
                ifs >> light.coords.x
                    >> light.coords.y
                    >> light.coords.z;
                this->lights.push_back(light);
            }
            else if(!tag.compare("supersample")) {
                std::string decisor;
                ifs >> decisor;
                decisor.compare("on") ? this->supersampling = false : this->supersampling = true;
            }
            else if(!tag.compare("profundidade")){
                ifs >> this->depth;
            }
            else if(!tag.compare("object")){
                double a, b, c, d, e, f, g, h, j, k, red, green, blue, ka, kd, ks, n, KS, KT, ir;
                ifs >> a
                    >> b
                    >> c
                    >> d
                    >> e
                    >> f
                    >> g
                    >> h
                    >> j
                    >> k
                    >> red
                    >> green
                    >> blue
                    >> ka
                    >> kd
                    >> ks
                    >> n
                    >> KS
                    >> KT
                    >> ir;
                Object object = Object(a, b, c, d, e, f, g, h, j, k, ka, kd, ks, n, KS, KT, ir, red, green, blue);
                this->objects.push_back(object);
            }
        }
    }
    ifs.close();
}

std::string SDL::getOutput() {
    return this->output;
}

T3 SDL::getEye(){
    return this->eye;
}

Ortho SDL::getOrtho(){
    return this->ortho;
}

Size SDL::getSize(){
    return this->size;
}

T3 SDL::getBackground(){
    return this->background;
}

double SDL::getAmbient(){
    return this->ambient;
}

std::vector<Light> SDL::getLights(){
    return this->lights;
}

bool SDL::getSuperSampling(){
    return this->supersampling;
}

double SDL::getDepth(){
    return this->depth;
}

std::vector<Object> SDL::getObjects(){
    return this->objects;
}