# rt
A implementation of a ray tracing for my college

To manipulate spheres just modify the sdl file, on the object field
g = -x0
h = -y0
j = -z0
R^2 = x0² + y0² + z0² - k

sphere: (x-x0)² + (y-y0)² + (z-z0)² = R²

To run, just compile the main.cpp (g++ main.cpp -o out -std=c++11 && ./out)
