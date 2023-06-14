// ######################## VTK LIBRARIES #######################
// Points
#include <vtkActor.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSimplePointsReader.h>
// Circle
#include <vtkPolyData.h>
#include <vtkRegularPolygonSource.h>

// ######################## MAIN LIBRARIES #######################
#include <iostream>
#include <vector>

#include <sstream>
#include <string>
#include <fstream>

#include "pca.h"

using namespace std;

vector<vector<float>> allCircleInformation;

template <int ndim>
struct Point{
    float puntos[ndim];

    Point(){
        for(int i = 0; i < ndim; i++){
            puntos[i] = 0;
        }
    }

    Point( float _puntos[] ){
        for(int i = 0; i < ndim; i++){
            puntos[i] = _puntos[i];
        }
    }

    void printPoint(){
        for(int i = 0; i < ndim; i++){
            cout << puntos[i] << '\t';
        }
        cout << endl;
    }

    float &operator[](int posicion){
        return puntos[posicion];
    }
};

template <typename T, int ndim>
struct Node{
    int nPuntos;
    int M;
    bool isLeaf;

    vector< Point<ndim> > grupoRecords;
    vector< Node<T, ndim> > grupoHijos;

    Point<ndim> puntoCentral;
    float radio;

    Node( int _maxRecords = 0 ){
        nPuntos = 0;
        M = _maxRecords;
        isLeaf = true;

        puntoCentral;
        radio = 0;
    }

    // ################### PRINT FUNCTIONS ###################################
    void printPointsNode( Node<T, ndim> *ptr ){
        for(int i = 0; i < ptr->nPuntos; i++){
            cout << i+1 << ". " << '\t';
            ptr->grupoRecords[i].printPoint();
        }
    }

    // ################### PRINT FUNCTIONS ###################################
    float getEuclideanDistance(Point<ndim> punto1, Point<ndim> punto2){
        float res = 0;
        for(int i = 0; i < ndim; i++){
            res = res + pow( punto2[i]-punto1[i] , 2);
        }
        return sqrt(res);
    }

    // ################### MAIN FUNCTIONS ###################################
    void getCircleInfo(Node<T, ndim> *ptr){
        // Obtener el punto central
        for(int i = 0; i < ndim; i++){
            float pntoCentralAux = 0;
            for(int j = 0; j < ptr->nPuntos; j++){
                pntoCentralAux = pntoCentralAux + ptr->grupoRecords[j][i];
            }
            pntoCentralAux = pntoCentralAux/ptr->nPuntos;
            puntoCentral[i] = pntoCentralAux;
        }

        // Obtener el radio
        float maxDistance = getEuclideanDistance( ptr->grupoRecords[0], ptr->puntoCentral);
        for(int i = 1; i < ptr->nPuntos; i++){
            float maxDistanceAux = getEuclideanDistance( ptr->grupoRecords[i],ptr->puntoCentral);
            if(maxDistanceAux > maxDistance){
                maxDistance = maxDistanceAux;
            }
        }
        ptr->radio = maxDistance;
        cout << "El radio es: " << ptr->radio << endl;
    }

    void insertToNode( Point<ndim> record, Node<T, ndim> *ptr ){
        ptr->grupoRecords.push_back(record);
        ptr->nPuntos++;
        ptr->isLeaf = false;
    }

    void createEstructure(Node *ptr){
        // --> Print all Points
        cout << "All points" << endl;
        printPointsNode(ptr);

        // 1. GetCircleInfo
        getCircleInfo(ptr);
        cout << "El punto centroal es: ";
        ptr->puntoCentral.printPoint();

        // 2. Split

    }

    void knnSearch( Point<ndim> recordBase , Node<T, ndim> *ptr, int Nvecinos){

    }

    void getMBRCircle( Node<T, ndim> *ptr ){
        if(ptr->isLeaf){
            return;
        }
        else{
            vector<float> circleAux;
            circleAux.push_back( ptr->puntoCentral[0] );
            circleAux.push_back( ptr->puntoCentral[1] );
            circleAux.push_back( 0 );
            circleAux.push_back( ptr->radio );

            allCircleInformation.push_back(circleAux);

            for(int i = 0; i < ptr->grupoHijos.size(); i++){
                cout << "Recorriendo hijo" << i+1 << endl;
                Node<T, ndim> *pntauxPrint;
                pntauxPrint = &ptr->grupoHijos[i];
                getMBRCircle(pntauxPrint);
            }
        }
    }
};

template <typename T, int ndim>
class BallStarTree {
private:
    Node<T, ndim> root;
    Node<T, ndim> *ptr;

public:
    BallStarTree(int _M){
        Node<T, ndim> aux(_M);
        root = aux;

        ptr = &root;
    }

    void insertToRoot( Point<ndim> newRecord ){
        root.insertToNode(newRecord, ptr);
    }

    void _createEstructure(){
        root.createEstructure(ptr);
    }

    void _knnSearch( string song_name, int Nvecinos ){
        Point<ndim> recordBase;
        // TO DO...
        // recordBase = getdata(song_name)
        root.knnSearch(recordBase, ptr, Nvecinos);
    }

    void getMBRCircle2(){
        root.getMBRCircle(ptr);
    }
};

vtkNew<vtkActor> drawCircle( float pointX, float pointY, float pointZ, float radius){
    vtkNew<vtkNamedColors> colors;

    // Create a circle
    vtkNew<vtkRegularPolygonSource> polygonSource;
    // Comment this line to generate a disk instead of a circle.
    polygonSource->GeneratePolygonOff();
    polygonSource->SetNumberOfSides(50);
    polygonSource->SetRadius(radius);
    polygonSource->SetCenter(pointX, pointY, pointZ);

    // Visualize
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(polygonSource->GetOutputPort());

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());

    return actor;
}

int main(int argc, char* argv[])
{
    // Verify input arguments
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " Filename(.xyz) e.g coords.txt"
                  << std::endl;
        return EXIT_FAILURE;
    }

    // ######################### BALL STAR TREE #################################
    const int M = 10; // maximo numero de puntos por nodo
    const int ndim = 3;
    //const int ndim = 12;

    cout << endl << "############################################### " << endl;
    cout << "               BALL * TREE " << endl;
    cout << "############################################## " << endl;
    BallStarTree<Point<ndim>, ndim> test1(M);

    // ########################## LECTURA .SVG ##################################

    string filename = "testData.csv";
    ifstream in(filename);
    string song_name;
    float pntoX, pntoY, pntoZ;

    while(  in >> pntoX >> pntoY >> pntoZ){
        float recordAux[] = {pntoX, pntoY, pntoZ};
        Point<ndim> record(recordAux);

        test1.insertToRoot(record);
    }

    /*
    string filename = "genres.csv";
    ifstream in(filename);
    string song_name;
    float   danceability, energy, key, loudness, speechiness,
            acousticness, instrumentalness, liveness, valence, tempo,
            duration_ms, time_signature;

    while(  in >> danceability >> energy >> key >> loudness >> speechiness >>
            acousticness >> instrumentalness >> liveness >> valence >> tempo >>
            duration_ms >>  time_signature){
        float recordAux[] = {danceability, energy, key, loudness, speechiness,
                             acousticness, instrumentalness, liveness, valence, tempo,
                             duration_ms, time_signature};
        Point<ndim> record(recordAux);

        test1.insertToRoot(record);
    }
    */

    // ######################## CREAR ESTRUCTURA ###############################
    cout << "Puntos insertados en el root ..." << endl;
    test1._createEstructure();

    // ###################### BUSCAR KNN MAS CERCANO ##########################
    int Nvecinos = 500;
    test1._knnSearch(song_name, Nvecinos);

    // ###################### CONFIG VTK #####################################
    vtkNew<vtkNamedColors> colors;

    // Read the file Point
    vtkNew<vtkSimplePointsReader> reader;
    reader->SetFileName(argv[1]);
    reader->Update();

    // Visualize
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(reader->GetOutputPort());

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetPointSize(10);
    actor->GetProperty()->SetColor(colors->GetColor3d("Tomato").GetData());

    vtkNew<vtkRenderer> renderer;
    renderer->AddActor(actor);

    test1.getMBRCircle2();
    for(int i = 0; i < allCircleInformation.size(); i++){
        float x = allCircleInformation[i][0];
        float y = allCircleInformation[i][1];
        float z = 0;
        float radius = allCircleInformation[i][3];

        renderer->AddActor(drawCircle(x,y,z,radius));
    }
    //renderer->AddActor(drawCircle(10,10,0,10));

    renderer->SetBackground(colors->GetColor3d("White").GetData());

    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->AddRenderer(renderer);
    renderWindow->SetWindowName("SimplePointsReader");

    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    renderWindowInteractor->SetRenderWindow(renderWindow);

    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}
