#include <GLFW/glfw3.h>
#include <utility>
#include "Support/Header.h"

using namespace std;
//eventos
void mbpressed(GLFWwindow* window, int button, int action, int mods){
}
void mmove(GLFWwindow* window, double xpos, double ypos){
}

int main(void){
    SDL* sdl = new SDL("Files/onesphere.sdl");
    cout << sdl->getOutput() << endl;
    cout << sdl->getObjects().at(0)->KT << endl;
    GLFWwindow* window;
    if (!glfwInit())return -1;
    window = glfwCreateWindow(640, 480,  "Ray Tracing", NULL, NULL);
    if (!window){
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetMouseButtonCallback(window, mbpressed);
    glfwSetCursorPosCallback(window, mmove);
    while (!glfwWindowShouldClose(window)){
        glClearColor(0,0,0, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);


        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    glfwTerminate();
    delete sdl;
    return 0;
}
