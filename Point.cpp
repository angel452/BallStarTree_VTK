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

// ######################## PCA LIBRARIES #######################
#include "pca.h"
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

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
            //cout << puntos[i] << '\t';
            cout << puntos[i] << "  ";
        }
        //cout << endl;
    }

    float &operator[](int posicion){
        return puntos[posicion];
    }
};

template <typename T, int ndim>
struct Node{
    int nPuntos; // numero de puntos en el nodo
    int M; // Maximo de modos por nodo
    bool isLeaf; // si es hoja o no

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
            //cout << i+1 << ". " << '\t';
            cout << "( " ;
            ptr->grupoRecords[i].printPoint();
            cout << ")";
        }
        cout << endl;
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
            float pntoCentralAux;
            float minPoint = ptr->grupoRecords[0][i];
            float maxPoint = ptr->grupoRecords[0][i];
            for(int j = 1; j < ptr->nPuntos; j++){
                float minPointAux = ptr->grupoRecords[j][i];
                float maxPointAux = ptr->grupoRecords[j][i];

                if(minPointAux < minPoint){
                    minPoint = minPointAux;
                }
                if(maxPointAux > maxPoint){
                    maxPoint = maxPointAux;
                }
                //pntoCentralAux = pntoCentralAux + ptr->grupoRecords[j][i];
            }
            pntoCentralAux = (maxPoint + minPoint)/2;
            //pntoCentralAux = pntoCentralAux/ptr->nPuntos;
            ptr->puntoCentral[i] = pntoCentralAux;
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

        /*
        // Mostrar la informacion
        cout << "El punto central es: "; ptr->puntoCentral.printPoint(); cout << endl;
        cout << "El radio es: " << ptr->radio << endl;
        cout << "El numero de puntos dentro es: " << ptr->nPuntos << endl;
        */
     }

    void insertToNode( Point<ndim> record, Node<T, ndim> *ptr ){
        ptr->grupoRecords.push_back(record);
        ptr->nPuntos++;
        ptr->isLeaf = false;
    }

    // ################### ALGORITMO ORDENAMIENTO QUICKSORT ###################################
    int particion(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> &reprojectPoints, int inicio, int final, Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> &OriginalPoints){
        float pivote = reprojectPoints(inicio, 0); // El pivote sera el primer elemento
        int aux = inicio + 1;
        // Recorremos desde el pivote hasta el final
        for(int i = aux; i <= final; i++){
            if(reprojectPoints(i, 0) < pivote){
                // Intercambiamos los valores de la matriz de puntos proyectados
                reprojectPoints.row(aux).swap(reprojectPoints.row(i));
                // Interacambiamos los valores de la matriz de puntos originales
                OriginalPoints.row(aux).swap(OriginalPoints.row(i));
                aux++;
            }
        }
        reprojectPoints.row(inicio).swap(reprojectPoints.row(aux-1));
        OriginalPoints.row(inicio).swap(OriginalPoints.row(aux-1));
        return aux-1;
    }

    void quickSort(Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> &reprojectPoints, int inicio, int final, Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> &OriginalPoints ){
        if(inicio < final){
            int pivote = particion(reprojectPoints, inicio, final, OriginalPoints);
            quickSort(reprojectPoints, inicio, pivote-1, OriginalPoints);
            quickSort(reprojectPoints, pivote+1, final, OriginalPoints);
        }
    }
    // ##########################################################################################

    pair<Node, Node> splitPCA(Node<T, ndim> *ptr){

        pair<Node, Node> resultSplit;

        // --> Variables:
        const int nPuntosAux = ptr->nPuntos;
        Eigen::MatrixXf pca_data_matrix(nPuntosAux, ndim);
        Eigen::MatrixXf pca_center_matrix(nPuntosAux, ndim);

        // --> Input Matrix:
        for(int i = 0; i < nPuntosAux; i++){
            for(int j = 0; j < ndim; j++){
                pca_data_matrix(i, j) = ptr->grupoRecords[i][j];
                //pca_data_matrix << ptr->grupoRecords[i][j];
            }
        }
        // --> PCA:
        pca_t<float> pca;
        pca.set_input(pca_data_matrix);
        pca.compute();

        /*
        std::cout
                << "Input Matrix:		\n" << pca.get_input_matrix() << std::endl << std::endl
                << "Centered Matrix:	\n" << pca.get_centered_matrix() << std::endl << std::endl
                << "Covariance Matrix:	\n" << pca.get_covariance_matrix() << std::endl << std::endl
                << "Projection Matrix:	\n" << pca.get_projection_matrix() << std::endl << std::endl
                << "Mean Matrix:		\n" << pca.get_mean() << std::endl << std::endl
                << "Eigen Values:		\n" << pca.get_eigen_values() << std::endl << std::endl
                << "Eigen Vectors:		\n" << pca.get_eigen_vectors() << std::endl << std::endl;
        */

        //const auto& reprojection = pca.reprojection2();
        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> reprojection = pca.reprojection2();
        //std::cout << "Reprojected Matrix:	\n" << reprojection << std::endl << std::endl;

        // --> Ordenar los puntos segun la primera columna de la matriz reprojection
        quickSort(reprojection, 0, (ptr->nPuntos)-1, pca_data_matrix);

        /*
        std::cout
                << "Reprojected Matrix Ordered:	\n" << reprojection << std::endl << std::endl
                << "Original Matrix Ordered:	\n" << pca_data_matrix << std::endl << std::endl;
        */

        // --> Separar hijos
        // 1. Crear hijo Izquierdo
        Node<T, ndim> leftChild;
        Node<T, ndim> *ptrLeftChild;
        Node<T, ndim> auxNode1(M);
        leftChild = auxNode1;
        ptrLeftChild = &leftChild;

        // 2. Crear hijo Derecho
        Node<T, ndim> rightChild;
        Node<T, ndim> *ptrRightChild;
        Node<T, ndim> auxNode2(M);
        rightChild = auxNode2;
        ptrRightChild = &rightChild;

        // --> Insertar puntos en los hijos
        // 1. Insertar puntos en el hijo izquierdo
        for(int i = 0; i < (ptr->nPuntos)/2; i++){
            Point<ndim> auxPoint;
            for(int j = 0; j < ndim; j++){
                auxPoint[j] = pca_data_matrix(i, j);
            }
            ptrLeftChild->grupoRecords.push_back(auxPoint);
            ptrLeftChild->nPuntos++;
            //insertToNode(auxPoint, ptrLeftChild);
        }

        // 2. Insertar puntos en el hijo derecho
        for(int i = (ptr->nPuntos)/2; i < ptr->nPuntos; i++){
            Point<ndim> auxPoint;
            for(int j = 0; j < ndim; j++){
                auxPoint[j] = pca_data_matrix(i, j);
            }
            ptrRightChild->grupoRecords.push_back(auxPoint);
            ptrRightChild->nPuntos++;
            //insertToNode(auxPoint, ptrRightChild);
        }

        resultSplit.first = leftChild;
        resultSplit.second = rightChild;

        return resultSplit;
    }

    void createEstructure(Node *ptr){
        // --> Print all Points
        cout << "Points in node..." << endl;
        printPointsNode(ptr);

        // 1. GetCircleInfo
        getCircleInfo(ptr);
        cout << "El punto central es: "; ptr->puntoCentral.printPoint(); cout << endl;
        cout << "El radio es: " << ptr->radio << endl;
        cout << "El numero de puntos dentro es: " << ptr->nPuntos << endl;

        // 2. Split
        cout << endl;

        if((ptr->nPuntos)/2 >= M){
            cout << "--- El nodo actual se puede dividir ---" << endl << endl;
            pair<Node, Node> resultSplit = splitPCA(ptr);

            // --> Insertamos los hijos al nodo actual
            ptr->grupoHijos.push_back(resultSplit.first);
            ptr->grupoHijos.push_back(resultSplit.second);

            // --> Creamos punteros a los hijos
            Node<T, ndim> *ptrLeftChild = &ptr->grupoHijos[0];
            Node<T, ndim> *ptrRightChild = &ptr->grupoHijos[1];

            ptrLeftChild->isLeaf = false;
            createEstructure(ptrLeftChild);

            ptrRightChild->isLeaf = false;
            createEstructure(ptrRightChild);

        }
        else{
            cout << "El nodo actual no se puede dividir mas " << endl;
        }

        cout << endl;
    }

    void knnSearch( Point<ndim> recordBase , Node<T, ndim> *ptr, int Nvecinos){

    }

    void getMBRCircle( Node<T, ndim> *ptr ){
        if(ptr->isLeaf){
            //cout << "Es hoja" << endl;
            vector<float> circleAux2;
            circleAux2.push_back( ptr->puntoCentral[0] );
            circleAux2.push_back( ptr->puntoCentral[1] );
            circleAux2.push_back( 0 );
            circleAux2.push_back( ptr->radio );

            allCircleInformation.push_back(circleAux2);

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
                //cout << "Recorriendo hijo" << i+1 << endl;
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
    const int M = 18; // maximo numero de puntos por nodo
    const int ndim = 3;
    //const int ndim = 12;

    cout << endl << "############################################### " << endl;
    cout << "               BALL * TREE " << endl;
    cout << "############################################## " << endl;
    BallStarTree<Point<ndim>, ndim> test1(M);

    // ########################## LECTURA .SVG ##################################

    string filename = "testData.csv";
    //string filename = "testData2.csv";
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
    cout << "La cantidad de circulos es: " << allCircleInformation.size() << endl;
    for(int i = 0; i < allCircleInformation.size(); i++){
        float x = allCircleInformation[i][0];
        float y = allCircleInformation[i][1];
        float z = 0;
        float radius = allCircleInformation[i][3];

        cout << "Info Circulo " << i+1 << ": " << x << " " << y << " " << z << " " << radius << endl;

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

/*
 * Canciones para buscar
 * Where Is The Love?
 * Bitch Better Have My Money
 * FourFiveSeconds
 * We Own It (Fast & Furious)
 */