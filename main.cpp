#include "Support/Header.h"

using namespace std;

int main(void){
    SDL sdl = SDL("Files/onesphere.sdl");
    Camera cam = Camera(0, 0 ,10, 0, 0, -1,0 ,1, 0,1, 0, 0,16,1,1);
    Object obj = Object(1, 1, 1, 0, 0, 0, 0, 0, 0, -300, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    cout<<"P1\n100 100\n";
    for(int i=0;i<100;i++)for(int j=0;j<100;j++){
            if(intersect(cam.coords,(cam.dir*cam.dist+cam.v*(cam.hx*(i-50)/100.0*2.0)+cam.u*(cam.hy*(j-50)/100.0*2.0)).normalizacao(),&obj)>0)cout<<0<<"\n";            else cout<<1<<"\n";
        }

    return 0;
}
