/*
 *  Linhas apos # sao comentarios:
 *  # ...
 *
 *  Nome do arquivo de saida eh definido com output
 *  output name
 *
 *  A posicao do observador virtual eh definido com eye
 *  eye x y z
 *  onde z > 0
 *
 *  A janela do universo eh definida com ortho
 *  ortho x0 y0 x1 y1
 *  onde x0,y0 sao as coordenadas do canto inferior esquerdo e x1, y1 do canto superior direito
 *
 *  A especificacao de como dividir a janela nos pixels individuais eh dada como
 *  size w h
 *  onde a janela dada pelo ortho deve ser dividida em w pixels de largura e h pixels de altura
 *
 *  A cor do fundo do cenario eh dada como
 *  background r g b
 *  onde os valores pertencem ao intervalo [0,1]
 *
 *  a intensidade da luz ambiente eh dada por
 *  ambiente Ia
 *  onde 0 <= Ia <= 1
 *
 *  A fonte de luz eh especificada como
 *  light x y z Ip
 *  onde x y z eh a direcao e Ip sua intensidade
 *  Uma mesma cena pode conter varias fontes de luz
 *
 *  Para especificar se a cena ultilizara supersampling
 *  supersample [on/off]
 *
 *  A profundidade recursiva eh dada por
 *  profundidade n
 *  exemplo:
 *      n = 0 : os raios param no primeiro objeto
 *      n = 1 : os raios percorrem uma vez a cena recursivamente
 *
 *  object a b c d e f g h j k red green blue ka kd ks n KS KT ir
 *  -Especifica uma superficie quadrica atraves dos coeficientes de a ate k.
 *  -A cor do objeto eh especificada pelas componentes red green blue, valores no intervalo [0, 1]
 *  -ka, kd, ks e n sao os coeficientes de reflexao do objeto
 *  -KS e KT sao os coeficientes de reflexao e transmissao globais do objeto, no intervalo [0, 1]
 *  -para objetos transparentes o coeficiente ir especifica o indice de refracao do objeto. O ir do ar eh 1
 *  -uma cena pode ter varios objetos
 */

#ifndef RT_SDLREADER_H
#define RT_SDLREADER_H


class SDL {
    std::string output;
    T3* eye;
    Ortho* ortho;
    Size* size;
    T3* background;
    double ambient;
    std::vector<Light*> lights;
    bool supersampling;
    double depth;
    std::vector<Object*> objects;
    void read();    //Method to read the file and fill the attribute fields

public:
    SDL(std::string path);
    ~SDL();
    std::string getOutput();
    T3* getEye();
    Ortho* getOrtho();
    Size* getSize();
    T3* getBackground();
    double getAmbient();
    std::vector<Light*> getLights();
    bool getSuperSampling();
    double getDepth();
    std::vector<Object*> getObjects();
};


#endif //RT_SDLREADER_H
