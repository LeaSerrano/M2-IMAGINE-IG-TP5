// -------------------------------------------
// gMini : a minimal OpenGL/GLUT application
// for 3D graphics.
// Copyright (C) 2006-2008 Tamy Boubekeur
// All rights reserved.
// -------------------------------------------

// -------------------------------------------
// Disclaimer: this code is dirty in the
// meaning that there is no attention paid to
// proper class attribute access, memory
// management or optimisation of any kind. It
// is designed for quick-and-dirty testing
// purpose.
// -------------------------------------------

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <GL/glut.h>
#include <float.h>
#include "src/Vec3.h"
#include "src/Camera.h"

enum DisplayMode{ WIRE=0, SOLID=1, LIGHTED_WIRE=2, LIGHTED=3 };

struct Triangle {
    inline Triangle () {
        v[0] = v[1] = v[2] = 0;
    }
    inline Triangle (const Triangle & t) {
        v[0] = t.v[0];   v[1] = t.v[1];   v[2] = t.v[2];
    }
    inline Triangle (unsigned int v0, unsigned int v1, unsigned int v2) {
        v[0] = v0;   v[1] = v1;   v[2] = v2;
    }
    unsigned int & operator [] (unsigned int iv) { return v[iv]; }
    unsigned int operator [] (unsigned int iv) const { return v[iv]; }
    inline virtual ~Triangle () {}
    inline Triangle & operator = (const Triangle & t) {
        v[0] = t.v[0];   v[1] = t.v[1];   v[2] = t.v[2];
        return (*this);
    }
    // membres indices des sommets du triangle:
    unsigned int v[3];
};


//Transformation made of a rotation and translation
struct Transformation {
    Mat3 rotation;
    Vec3 translation;
};

bool contain(std::vector<unsigned int> const & i_vector, unsigned int element) {
    for (unsigned int i = 0; i < i_vector.size(); i++) {
        if (i_vector[i] == element) return true;
    }
    return false;
}

void collect_one_ring (std::vector<Vec3> const & i_vertices,
                       std::vector< Triangle > const & i_triangles,
                       std::vector<std::vector<unsigned int> > & o_one_ring) {
    o_one_ring.clear();
    o_one_ring.resize(i_vertices.size()); //one-ring of each vertex, i.e. a list of vertices with which it shares an edge
    //Parcourir les triangles et ajouter les voisins dans le 1-voisinage
    //Attention verifier que l'indice n'est pas deja present
    for (unsigned int i = 0; i < i_triangles.size(); i++) {
        //Tous les points opposés dans le triangle sont reliés
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                if (j != k) {
                    if (!contain(o_one_ring[i_triangles[i][j]], i_triangles[i][k])) {
                        o_one_ring[i_triangles[i][j]].push_back(i_triangles[i][k]);
                    }
                }
            }
        }
    }
}


struct Mesh {
    std::vector< Vec3 > vertices; //array of mesh vertices positions
    std::vector< Vec3 > normals; //array of vertices normals useful for the display
    std::vector< Triangle > triangles; //array of mesh triangles
    std::vector< Vec3 > triangle_normals; //triangle normals to display face normals
    std::vector<float> vunicurvature_;
    std::vector<Vec3> Laplacien;
    std::vector<float> vcurvature_;
    std::vector<Vec3> LaplacienBeltrami;
    std::vector<std::vector<float>> edge_weights;
    std::vector<float> vertex_weights;
    std::vector<float> tshape_;
    std::vector<float> vgausscurvature_;

    //Compute face normals for the display
    void computeTrianglesNormals(){

        //A faire : implémenter le calcul des normales par face
        //Attention commencer la fonction par triangle_normals.clear();
        //Iterer sur les triangles

        //La normal du triangle i est le resultat du produit vectoriel de deux ses arêtes e_10 et e_20 normalisé (e_10^e_20)
        //L'arete e_10 est représentée par le vecteur partant du sommet 0 (triangles[i][0]) au sommet 1 (triangles[i][1])
        //L'arete e_20 est représentée par le vecteur partant du sommet 0 (triangles[i][0]) au sommet 2 (triangles[i][2])

        //Normaliser et ajouter dans triangle_normales

        triangle_normals.clear();
        for( unsigned int i = 0 ; i < triangles.size() ;i++ ){
            const Vec3 & e0 = vertices[triangles[i][1]] - vertices[triangles[i][0]];
            const Vec3 & e1 = vertices[triangles[i][2]] - vertices[triangles[i][0]];
            Vec3 n = Vec3::cross( e0, e1 );
            n.normalize();
            triangle_normals.push_back( n );
        }
    }

    //Compute vertices normals as the average of its incident faces normals
    void computeVerticesNormals(  ){
        //Utiliser weight_type : 0 uniforme, 1 aire des triangles, 2 angle du triangle

        //A faire : implémenter le calcul des normales par sommet comme la moyenne des normales des triangles incidents
        //Attention commencer la fonction par normals.clear();
        //Initializer le vecteur normals taille vertices.size() avec Vec3(0., 0., 0.)
        //Iterer sur les triangles

        //Pour chaque triangle i
        //Ajouter la normal au triangle à celle de chacun des sommets en utilisant des poids
        //0 uniforme, 1 aire du triangle, 2 angle du triangle

        //Iterer sur les normales et les normaliser
        normals.clear();
        normals.resize( vertices.size(), Vec3(0., 0., 0.) );
        for( unsigned int i = 0 ; i < triangles.size() ;i++ ){
            for( unsigned int t = 0 ; t < 3 ; t++ )
                normals[ triangles[i][t] ] += triangle_normals[i];
        }
        for( unsigned int i = 0 ; i < vertices.size() ;i++ )
            normals[ i ].normalize();


    }

    void computeNormals(){
        computeTrianglesNormals();
        computeVerticesNormals();
    }

    void calc_uniform_mean_curvature() {
        vunicurvature_.clear();
        Laplacien.clear();

        std::vector<std::vector<unsigned int> > oneRing;
        collect_one_ring (vertices, triangles, oneRing);

        Vec3 somme;
        Vec3 L;

        for (int i = 0; i < vertices.size(); i++) {
            somme = Vec3(0, 0, 0);

            for (int j = 0; j < oneRing[i].size(); j++) {

                somme += vertices[oneRing[i][j]];
            }

            L = (somme/oneRing[i].size()) - vertices[i];
            vunicurvature_.push_back((float)(L.length()/2.0));
            Laplacien.push_back(L);
        }

    }

    void uniform_smooth(unsigned int _iters) {
        for (int i = 0; i < _iters; i++) {
            std::vector<Vec3> vPrime;
            vPrime.resize(vertices.size());
            calc_uniform_mean_curvature();

            for (int j = 0; j < vertices.size(); j++) {
                vPrime[j] = (vertices[j] + ((1.0/2.0)*Laplacien[j]));
            }

            vertices = vPrime;
        }
        computeNormals();
    }

    void taubinSmooth(float lambda, float mu, unsigned int _iters) {
        for (int i = 0; i < _iters; i++) {
            std::vector<Vec3> vPrime;
            vPrime.resize(vertices.size());
            calc_uniform_mean_curvature();

            bool iter = true;
            float mult;
    
            for (int j = 0; j < vertices.size(); j++) {
                if (iter) {
                    mult = lambda;
                } else {
                    mult = mu;
                }

                vPrime[j] = (vertices[j] + mult * Laplacien[j]);
            }

            vertices = vPrime;
            iter = !iter;
        }
        computeNormals();
    }


    void addNoise() {
        for (unsigned int i = 0; i < vertices.size(); i++) {
            float factor = 0.03;

            const Vec3 & p = vertices[i];
            const Vec3 & n = normals[i];

            vertices[i] = Vec3(p[0] + factor * ((double)(rand())/(double)(RAND_MAX))*n[0], p[1] + factor * ((double)(rand())/(double)(RAND_MAX))*n[1], p[2] + factor * ((double)(rand())/(double)(RAND_MAX))*n[2]);
        }
    }

    void calc_triangle_quality() {

        tshape_.resize(triangles.size());
        tshape_.clear();

        for (int i = 0; i < triangles.size(); ++i) {
            Vec3 a = vertices[triangles[i].v[0]];
            Vec3 b = vertices[triangles[i].v[1]];
            Vec3 c = vertices[triangles[i].v[2]];

            Vec3 ab, ac;
            ab[0] = b[0] - a[0];
            ab[1] = b[1] - a[1];
            ab[2] = b[2] - a[2];
    
            ac[0] = c[0] - a[0];
            ac[1] = c[1] - a[1];
            ac[2] = c[2] - a[2];
    
            Vec3 normal = Vec3::cross(ab, ac);

            double area;
            if (normal.norm() < 1e-6) {
                area = 1e6;
            } else {
                area = normal.norm()/2.0;
            }

            tshape_.push_back(area);
        }
    }

    void calc_weights() {

        edge_weights.clear();
        vertex_weights.clear();
        edge_weights.resize(vertices.size(), std::vector<float>(vertices.size(), 0.0f));
        vertex_weights.resize(vertices.size(), 0.0f);


        for( unsigned int t = 0 ; t < triangles.size() ; ++t )
        {
            unsigned int v0 = triangles[t][0];
            unsigned int v1 = triangles[t][1];
            unsigned int v2 = triangles[t][2];

            Vec3 p0( vertices[v0] );
            Vec3 p1( vertices[v1] );
            Vec3 p2( vertices[v2] );

            double p0p1_slength = (p1-p0).sqrnorm();
            double p1p2_slength = (p2-p1).sqrnorm();
            double p2p0_slength = (p0-p2).sqrnorm();

            double dot0 = Vec3::dot(p1-p0,p2-p0);
            double dot1 = Vec3::dot(p0-p1,p2-p1);
            double dot2 = Vec3::dot(p0-p2,p1-p2);

            if( dot0 < 0.0 )
            {
                Vec3 const & fakeCircumcenter = (p1+p2)/2.0;

                double edge02Weight = sqrt( ( (p0+p2)/2.0 - fakeCircumcenter ).sqrnorm()  /  p2p0_slength );

                edge_weights[v0][v2] += edge02Weight;
                edge_weights[v2][v0] += edge02Weight;

                double edge01Weight = sqrt( ( (p0+p1)/2.0 - fakeCircumcenter ).sqrnorm()  /  p0p1_slength );
                edge_weights[v0][v1] += edge01Weight;
                edge_weights[v1][v0] += edge01Weight;

                double t_area = Vec3::cross( p1 - p0 , p2 - p0 ).norm() / 2.0;

                vertex_weights[v0] += t_area/2.0;
                vertex_weights[v1] += t_area/4.0;
                vertex_weights[v2] += t_area/4.0;

            }
            else if(dot1 < 0.0)
            {
                Vec3 const & fakeCircumcenter = (p0+p2)/2.0;

                double edge12Weight = sqrt( ( (p2+p1)/2.0 - fakeCircumcenter ).sqrnorm()  /  p1p2_slength );
                edge_weights[v1][v2] += edge12Weight;
                edge_weights[v2][v1] += edge12Weight;

                double edge01Weight = sqrt( ( (p0+p1)/2.0 - fakeCircumcenter ).sqrnorm()  /  p0p1_slength );
                edge_weights[v0][v1] += edge01Weight;
                edge_weights[v1][v0] += edge01Weight;

                double t_area = Vec3::cross( p1 - p0 , p2 - p0 ).norm() / 2.0;

                vertex_weights[v0] += t_area/4.0;
                vertex_weights[v1] += t_area/2.0;
                vertex_weights[v2] += t_area/4.0;
            }
            else if(dot2 < 0.0)
            {
                Vec3 const & fakeCircumcenter = (p0+p1)/2.0;

                double edge12Weight = sqrt( ( (p2+p1)/2.0 - fakeCircumcenter ).sqrnorm()  /  p1p2_slength );
                edge_weights[v1][v2] += edge12Weight;
                edge_weights[v2][v1] += edge12Weight;

                double edge02Weight = sqrt( ( (p0+p2)/2.0 - fakeCircumcenter ).sqrnorm()  /  p2p0_slength );
                edge_weights[v0][v2] += edge02Weight;
                edge_weights[v2][v0] += edge02Weight;

                double t_area = Vec3::cross( p1 - p0 , p2 - p0 ).norm() / 2.0;

                vertex_weights[v0] += t_area/4.0;
                vertex_weights[v1] += t_area/4.0;
                vertex_weights[v2] += t_area/2.0;
            }
            else
            {
                double cotW0_by_2 = dot0  / ( 2.0 * sqrt(p0p1_slength*p2p0_slength-dot0*dot0) );
                double cotW1_by_2 = dot1  / ( 2.0 * sqrt(p1p2_slength*p0p1_slength-dot1*dot1) );
                double cotW2_by_2 = dot2  / ( 2.0 * sqrt(p1p2_slength*p2p0_slength-dot2*dot2) );

                edge_weights[v1][v2] += cotW0_by_2;
                edge_weights[v2][v1] += cotW0_by_2;

                edge_weights[v0][v2] += cotW1_by_2;
                edge_weights[v2][v0] += cotW1_by_2;

                edge_weights[v1][v0] += cotW2_by_2;
                edge_weights[v0][v1] += cotW2_by_2;

                vertex_weights[v1] += cotW0_by_2 * p1p2_slength/2.0;
                vertex_weights[v2] += cotW0_by_2 * p1p2_slength/2.0;

                vertex_weights[v0] += cotW1_by_2 * p2p0_slength/2.0;
                vertex_weights[v2] += cotW1_by_2 * p2p0_slength/2.0;

                vertex_weights[v0] += cotW2_by_2 * p0p1_slength/2.0;
                vertex_weights[v1] += cotW2_by_2 * p0p1_slength/2.0;
            }
        }

        double total_weight_sum = 0.0;
        for (unsigned int v = 0; v < vertices.size(); ++v) {
                total_weight_sum += vertex_weights[v];
        }

        for (unsigned int v = 0; v < vertices.size(); ++v) {
                vertex_weights[v] /= total_weight_sum;
        }

        total_weight_sum = 0.0;
        for (unsigned int v0 = 0; v0 < vertices.size(); ++v0) {
            for (unsigned int v1 = 0; v1 < vertices.size(); ++v1) {
                total_weight_sum += edge_weights[v0][v1];
            }
        }

        for (unsigned int v0 = 0; v0 < vertices.size(); ++v0) {
            for (unsigned int v1 = 0; v1 < vertices.size(); ++v1) {
                edge_weights[v0][v1] /= total_weight_sum;
            }
        }
    }

    void calc_mean_curvature() {
        vcurvature_.clear();
        LaplacienBeltrami.clear();
        vertex_weights.clear();

        std::vector<std::vector<unsigned int>> oneRing;
        collect_one_ring (vertices, triangles, oneRing);

        Vec3 mean;
        float sum;

        calc_weights();

        for (int i = 0; i < vertices.size(); i++) {
            mean = Vec3(0, 0, 0);
            sum = 0;

            for (int j = 0; j < oneRing[i].size(); j++) {
                //for (int k = 0; k < oneRing[i].size(); k++) {

                mean += vertex_weights[oneRing[i][j]]*(vertices[oneRing[i][j]]-vertices[i]);
                sum += vertex_weights[oneRing[i][j]];

                //mean += edge_weights[oneRing[i][j]][oneRing[i][k]]*(vertices[oneRing[i][j]]-vertices[i]);
                //sum += edge_weights[oneRing[i][j]][oneRing[i][k]];
                //}
            }

            vcurvature_.push_back(mean.length()/2.0);
            LaplacienBeltrami.push_back((1.0/sum)*mean);
        }
    }

    void uniform_smooth_LaplaceBeltrami(unsigned int _iters) {
        for (int i = 0; i < _iters; i++) {
            std::vector<Vec3> vPrime;
            vPrime.resize(vertices.size());
            calc_mean_curvature();

            for (int j = 0; j < vertices.size(); j++) {
               vPrime[j] = (vertices[j] + ((1.0/2.0)*LaplacienBeltrami[j]));
            }

            vertices = vPrime;
        }
        computeNormals();
    }

    void calc_gauss_curvature() {
        vgausscurvature_.clear();

        std::vector<std::vector<unsigned int> > oneRing;
        collect_one_ring (vertices, triangles, oneRing);

        float somme;
        Vec3 L;

        for (int i = 0; i < vertices.size(); i++) {
            somme = 0.0;

            for (int j = 0; j < oneRing[i].size(); j+=2) {
                float angle;

                if (j == oneRing[i].size() - 1) {
                    float length1 = vertices[oneRing[i][j]].norm();
                    float length2 = vertices[oneRing[i][0]].norm();

                    angle = std::acos((vertices[oneRing[i][j]][0] * vertices[oneRing[i][0]][0] + vertices[oneRing[i][j]][1] * vertices[oneRing[i][0]][1] + vertices[oneRing[i][j]][2] * vertices[oneRing[i][0]][2]) / (length1 * length2));
                }
                else {
                    float length1 = vertices[oneRing[i][j]].norm();
                    float length2 = vertices[oneRing[i][j+1]].norm();

                    angle = std::acos((vertices[oneRing[i][j]][0] * vertices[oneRing[i][j+1]][0] + vertices[oneRing[i][j]][1] * vertices[oneRing[i][j+1]][1] + vertices[oneRing[i][j]][2] * vertices[oneRing[i][j+1]][2]) / (length1 * length2));
                
                }

                somme += angle;
            }

            vgausscurvature_.push_back(2*M_PI - somme);
        }

    }
};

void getMinMax(float &min, float &max, std::vector<float> p) {
    min = FLT_MAX;
    max = -FLT_MAX;

    for (int i = 0; i < p.size(); i++) {
        if (p[i] > max) {
            max = p[i];
        }
        else if (p[i] < min) {
            min = p[i];
        }
    }
}

//Input mesh loaded at the launch of the application
Mesh mesh;
std::vector< float > current_field; //normalized filed of each vertex

bool display_normals;
bool display_smooth_normals;
bool display_mesh;

DisplayMode displayMode;
int weight_type;

void updateCurvature(std::vector<float> curvature) {
    float min, max;
    getMinMax(min, max, curvature);
    current_field.clear();

    for (int i = 0; i < curvature.size(); i++) {
        current_field.push_back((curvature[i]-min)/(max-min));
    }
}

// -------------------------------------------
// OpenGL/GLUT application code.
// -------------------------------------------

static GLint window;
static unsigned int SCREENWIDTH = 900;
static unsigned int SCREENHEIGHT = 900;
static Camera camera;
static bool mouseRotatePressed = false;
static bool mouseMovePressed = false;
static bool mouseZoomPressed = false;
static int lastX=0, lastY=0, lastZoom=0;
static bool fullScreen = false;

// ------------------------------------
// File I/O
// ------------------------------------
bool saveOFF( const std::string & filename ,
              std::vector< Vec3 > const & i_vertices ,
              std::vector< Vec3 > const & i_normals ,
              std::vector< Triangle > const & i_triangles,
              std::vector< Vec3 > const & i_triangle_normals ,
              bool save_normals = false ) {
    std::ofstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open()) {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    myfile << "OFF" << std::endl ;

    unsigned int n_vertices = i_vertices.size() , n_triangles = i_triangles.size();
    myfile << n_vertices << " " << n_triangles << " 0" << std::endl;

    for( unsigned int v = 0 ; v < n_vertices ; ++v ) {
        myfile << i_vertices[v][0] << " " << i_vertices[v][1] << " " << i_vertices[v][2] << " ";
        if (save_normals) myfile << i_normals[v][0] << " " << i_normals[v][1] << " " << i_normals[v][2] << std::endl;
        else myfile << std::endl;
    }
    for( unsigned int f = 0 ; f < n_triangles ; ++f ) {
        myfile << 3 << " " << i_triangles[f][0] << " " << i_triangles[f][1] << " " << i_triangles[f][2]<< " ";
        if (save_normals) myfile << i_triangle_normals[f][0] << " " << i_triangle_normals[f][1] << " " << i_triangle_normals[f][2];
        myfile << std::endl;
    }
    myfile.close();
    return true;
}

void openOFF( std::string const & filename,
              std::vector<Vec3> & o_vertices,
              std::vector<Vec3> & o_normals,
              std::vector< Triangle > & o_triangles,
              std::vector< Vec3 > & o_triangle_normals,
              bool load_normals = true )
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return;
    }

    std::string magic_s;

    myfile >> magic_s;

    if( magic_s != "OFF" )
    {
        std::cout << magic_s << " != OFF :   We handle ONLY *.off files." << std::endl;
        myfile.close();
        exit(1);
    }

    int n_vertices , n_faces , dummy_int;
    myfile >> n_vertices >> n_faces >> dummy_int;

    o_vertices.clear();
    o_normals.clear();

    for( int v = 0 ; v < n_vertices ; ++v )
    {
        float x , y , z ;

        myfile >> x >> y >> z ;
        o_vertices.push_back( Vec3( x , y , z ) );

        if( load_normals ) {
            myfile >> x >> y >> z;
            o_normals.push_back( Vec3( x , y , z ) );
        }
    }

    o_triangles.clear();
    o_triangle_normals.clear();
    for( int f = 0 ; f < n_faces ; ++f )
    {
        int n_vertices_on_face;
        myfile >> n_vertices_on_face;

        if( n_vertices_on_face == 3 )
        {
            unsigned int _v1 , _v2 , _v3;
            myfile >> _v1 >> _v2 >> _v3;

            o_triangles.push_back(Triangle( _v1, _v2, _v3 ));

            if( load_normals ) {
                float x , y , z ;
                myfile >> x >> y >> z;
                o_triangle_normals.push_back( Vec3( x , y , z ) );
            }
        }
        else if( n_vertices_on_face == 4 )
        {
            unsigned int _v1 , _v2 , _v3 , _v4;
            myfile >> _v1 >> _v2 >> _v3 >> _v4;

            o_triangles.push_back(Triangle(_v1, _v2, _v3 ));
            o_triangles.push_back(Triangle(_v1, _v3, _v4));
            if( load_normals ) {
                float x , y , z ;
                myfile >> x >> y >> z;
                o_triangle_normals.push_back( Vec3( x , y , z ) );
            }

        }
        else
        {
            std::cout << "We handle ONLY *.off files with 3 or 4 vertices per face" << std::endl;
            myfile.close();
            exit(1);
        }
    }

}

// ------------------------------------
// Application initialization
// ------------------------------------
void initLight () {
    GLfloat light_position1[4] = {22.0f, 16.0f, 50.0f, 0.0f};
    GLfloat direction1[3] = {-52.0f,-16.0f,-50.0f};
    GLfloat color1[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat ambient[4] = {0.3f, 0.3f, 0.3f, 0.5f};

    glLightfv (GL_LIGHT1, GL_POSITION, light_position1);
    glLightfv (GL_LIGHT1, GL_SPOT_DIRECTION, direction1);
    glLightfv (GL_LIGHT1, GL_DIFFUSE, color1);
    glLightfv (GL_LIGHT1, GL_SPECULAR, color1);
    glLightModelfv (GL_LIGHT_MODEL_AMBIENT, ambient);
    glEnable (GL_LIGHT1);
    glEnable (GL_LIGHTING);
}

void init () {
    camera.resize (SCREENWIDTH, SCREENHEIGHT);
    initLight ();
    glCullFace (GL_BACK);
    glDisable (GL_CULL_FACE);
    glDepthFunc (GL_LESS);
    glEnable (GL_DEPTH_TEST);
    glClearColor (0.2f, 0.2f, 0.3f, 1.0f);
    glEnable(GL_COLOR_MATERIAL);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    display_normals = false;
    display_mesh = true;
    display_smooth_normals = true;
    displayMode = LIGHTED;

}


// ------------------------------------
// Rendering.
// ------------------------------------

void drawVector( Vec3 const & i_from, Vec3 const & i_to ) {

    glBegin(GL_LINES);
    glVertex3f( i_from[0] , i_from[1] , i_from[2] );
    glVertex3f( i_to[0] , i_to[1] , i_to[2] );
    glEnd();
}

void drawAxis( Vec3 const & i_origin, Vec3 const & i_direction ) {

    glLineWidth(4); // for example...
    drawVector(i_origin, i_origin + i_direction);
}

void drawReferenceFrame( Vec3 const & origin, Vec3 const & i, Vec3 const & j, Vec3 const & k ) {

    glDisable(GL_LIGHTING);
    glColor3f( 0.8, 0.2, 0.2 );
    drawAxis( origin, i );
    glColor3f( 0.2, 0.8, 0.2 );
    drawAxis( origin, j );
    glColor3f( 0.2, 0.2, 0.8 );
    drawAxis( origin, k );
    glEnable(GL_LIGHTING);

}


typedef struct {
    float r;       // ∈ [0, 1]
    float g;       // ∈ [0, 1]
    float b;       // ∈ [0, 1]
} RGB;



RGB scalarToRGB( float scalar_value ) //Scalar_value ∈ [0, 1]
{
    RGB rgb;
    float H = scalar_value*360., S = 1., V = 0.85,
            P, Q, T,
            fract;

    (H == 360.)?(H = 0.):(H /= 60.);
    fract = H - floor(H);

    P = V*(1. - S);
    Q = V*(1. - S*fract);
    T = V*(1. - S*(1. - fract));

    if      (0. <= H && H < 1.)
        rgb = (RGB){.r = V, .g = T, .b = P};
    else if (1. <= H && H < 2.)
        rgb = (RGB){.r = Q, .g = V, .b = P};
    else if (2. <= H && H < 3.)
        rgb = (RGB){.r = P, .g = V, .b = T};
    else if (3. <= H && H < 4.)
        rgb = (RGB){.r = P, .g = Q, .b = V};
    else if (4. <= H && H < 5.)
        rgb = (RGB){.r = T, .g = P, .b = V};
    else if (5. <= H && H < 6.)
        rgb = (RGB){.r = V, .g = P, .b = Q};
    else
        rgb = (RGB){.r = 0., .g = 0., .b = 0.};

    return rgb;
}

void drawSmoothTriangleMesh( Mesh const & i_mesh , bool draw_field = false ) {
    glBegin(GL_TRIANGLES);
    for(unsigned int tIt = 0 ; tIt < i_mesh.triangles.size(); ++tIt) {

        for(unsigned int i = 0 ; i < 3 ; i++) {
            const Vec3 & p = i_mesh.vertices[i_mesh.triangles[tIt][i]]; //Vertex position
            const Vec3 & n = i_mesh.normals[i_mesh.triangles[tIt][i]]; //Vertex normal

            if( draw_field && current_field.size() > 0 ){
                RGB color = scalarToRGB( current_field[i_mesh.triangles[tIt][i]] );
                glColor3f( color.r, color.g, color.b );
            }
            glNormal3f( n[0] , n[1] , n[2] );
            glVertex3f( p[0] , p[1] , p[2] );
        }
    }
    glEnd();

}

void drawTriangleMesh( Mesh const & i_mesh , bool draw_field = false  ) {
    glBegin(GL_TRIANGLES);
    for(unsigned int tIt = 0 ; tIt < i_mesh.triangles.size(); ++tIt) {
        const Vec3 & n = i_mesh.triangle_normals[ tIt ]; //Triangle normal
        for(unsigned int i = 0 ; i < 3 ; i++) {
            const Vec3 & p = i_mesh.vertices[i_mesh.triangles[tIt][i]]; //Vertex position

            if( draw_field ){
                RGB color = scalarToRGB( current_field[i_mesh.triangles[tIt][i]] );
                glColor3f( color.r, color.g, color.b );
            }
            glNormal3f( n[0] , n[1] , n[2] );
            glVertex3f( p[0] , p[1] , p[2] );
        }
    }
    glEnd();

}

void drawMesh( Mesh const & i_mesh , bool draw_field = false ){
    if(display_smooth_normals)
        drawSmoothTriangleMesh(i_mesh, draw_field) ; //Smooth display with vertices normals
    else
        drawTriangleMesh(i_mesh, draw_field) ; //Display with face normals
}

void drawVectorField( std::vector<Vec3> const & i_positions, std::vector<Vec3> const & i_directions ) {
    glLineWidth(1.);
    for(unsigned int pIt = 0 ; pIt < i_directions.size() ; ++pIt) {
        Vec3 to = i_positions[pIt] + 0.02*i_directions[pIt];
        drawVector(i_positions[pIt], to);
    }
}

void drawNormals(Mesh const& i_mesh){

    if(display_smooth_normals){
        drawVectorField( i_mesh.vertices, i_mesh.normals );
    } else {
        std::vector<Vec3> triangle_baricenters;
        for ( const Triangle& triangle : i_mesh.triangles ){
            Vec3 triangle_baricenter (0.,0.,0.);
            for( unsigned int i = 0 ; i < 3 ; i++ )
                triangle_baricenter += i_mesh.vertices[triangle[i]];
            triangle_baricenter /= 3.;
            triangle_baricenters.push_back(triangle_baricenter);
        }

        drawVectorField( triangle_baricenters, i_mesh.triangle_normals );
    }
}

//Draw fonction
void draw () {
    glShadeModel(GL_FLAT);

    if(displayMode == LIGHTED || displayMode == LIGHTED_WIRE){

        glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
        glEnable(GL_LIGHTING);

    }  else if(displayMode == WIRE){

        glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
        glDisable (GL_LIGHTING);

    }  else if(displayMode == SOLID ){
        glDisable (GL_LIGHTING);
        glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

    }

    glColor3f(0.8,1,0.8);
    drawMesh(mesh, true);

    if(displayMode == SOLID || displayMode == LIGHTED_WIRE){
        glEnable (GL_POLYGON_OFFSET_LINE);
        glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
        glLineWidth (1.0f);
        glPolygonOffset (-2.0, 1.0);

        glColor3f(0.,0.,0.);
        drawMesh(mesh, false);

        glDisable (GL_POLYGON_OFFSET_LINE);
        glEnable (GL_LIGHTING);
    }



    glDisable(GL_LIGHTING);
    if(display_normals){
        glColor3f(1.,0.,0.);
        drawNormals(mesh);
    }

    glEnable(GL_LIGHTING);


}

void changeDisplayMode(){
    if(displayMode == LIGHTED)
        displayMode = LIGHTED_WIRE;
    else if(displayMode == LIGHTED_WIRE)
        displayMode = SOLID;
    else if(displayMode == SOLID)
        displayMode = WIRE;
    else
        displayMode = LIGHTED;
}

void display () {
    glLoadIdentity ();
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera.apply ();
    draw ();
    glFlush ();
    glutSwapBuffers ();
}

void idle () {
    glutPostRedisplay ();
}

// ------------------------------------
// User inputs
// ------------------------------------
//Keyboard event
int it = 0;

void key (unsigned char keyPressed, int x, int y) {
    switch (keyPressed) {
    case 'f':
        if (fullScreen == true) {
            glutReshapeWindow (SCREENWIDTH, SCREENHEIGHT);
            fullScreen = false;
        } else {
            glutFullScreen ();
            fullScreen = true;
        }
        break;


    case 'w':
        changeDisplayMode();
        break;

    case 'a' :
        mesh.uniform_smooth(30);
        break;

    case 'b' :
        mesh.taubinSmooth(0.330, -0.331, 30);
        break;

    case 'c' : 
        mesh.uniform_smooth_LaplaceBeltrami(30);
        break;

    case 'd' : 
        mesh.calc_uniform_mean_curvature();
        updateCurvature(mesh.vunicurvature_);
        break;

    case 'e' : 
        mesh.calc_mean_curvature();
        updateCurvature(mesh.vcurvature_);
        break;

    case 'g' : 
        mesh.calc_triangle_quality();
        updateCurvature(mesh.tshape_);
        break;

    case 'h' : 
        mesh.calc_gauss_curvature();
        updateCurvature(mesh.vgausscurvature_);
        break;

    case 'n': //Press n key to display normals
        display_normals = !display_normals;
        break;

    case '1': //Toggle loaded mesh display
        display_mesh = !display_mesh;
        break;

    case 's': //Switches between face normals and vertices normals
        display_smooth_normals = !display_smooth_normals;
        break;

    case '+': //Changes weight type: 0 uniforme, 1 aire des triangles, 2 angle du triangle
        weight_type ++;
        if(weight_type == 3) weight_type = 0;
        mesh.computeVerticesNormals(); //recalcul des normales avec le type de poids choisi
        break;

    default:
        break;
    }
    idle ();
}

//Mouse events
void mouse (int button, int state, int x, int y) {
    if (state == GLUT_UP) {
        mouseMovePressed = false;
        mouseRotatePressed = false;
        mouseZoomPressed = false;
    } else {
        if (button == GLUT_LEFT_BUTTON) {
            camera.beginRotate (x, y);
            mouseMovePressed = false;
            mouseRotatePressed = true;
            mouseZoomPressed = false;
        } else if (button == GLUT_RIGHT_BUTTON) {
            lastX = x;
            lastY = y;
            mouseMovePressed = true;
            mouseRotatePressed = false;
            mouseZoomPressed = false;
        } else if (button == GLUT_MIDDLE_BUTTON) {
            if (mouseZoomPressed == false) {
                lastZoom = y;
                mouseMovePressed = false;
                mouseRotatePressed = false;
                mouseZoomPressed = true;
            }
        }
    }

    idle ();
}

//Mouse motion, update camera
void motion (int x, int y) {
    if (mouseRotatePressed == true) {
        camera.rotate (x, y);
    }
    else if (mouseMovePressed == true) {
        camera.move ((x-lastX)/static_cast<float>(SCREENWIDTH), (lastY-y)/static_cast<float>(SCREENHEIGHT), 0.0);
        lastX = x;
        lastY = y;
    }
    else if (mouseZoomPressed == true) {
        camera.zoom (float (y-lastZoom)/SCREENHEIGHT);
        lastZoom = y;
    }
}


void reshape(int w, int h) {
    camera.resize (w, h);
}

// ------------------------------------
// Start of graphical application
// ------------------------------------
int main (int argc, char ** argv) {
    if (argc > 2) {
        exit (EXIT_FAILURE);
    }
    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize (SCREENWIDTH, SCREENHEIGHT);
    window = glutCreateWindow ("TP HAI917I");

    init ();
    glutIdleFunc (idle);
    glutDisplayFunc (display);
    glutKeyboardFunc (key);
    glutReshapeFunc (reshape);
    glutMotionFunc (motion);
    glutMouseFunc (mouse);
    key ('?', 0, 0);

    //Mesh loaded with precomputed normals
    openOFF("data/elephant_n.off", mesh.vertices, mesh.normals, mesh.triangles, mesh.triangle_normals);

    mesh.computeNormals();
    mesh.addNoise();

    // A faire : normaliser les champs pour avoir une valeur flotante entre 0. et 1. dans current_field
    //***********************************************//

    current_field.clear();

    //mesh.calc_uniform_mean_curvature();
    //mesh.calc_triangle_quality();
    //mesh.calc_gauss_curvature();
    //updateCurvature(mesh.vgausscurvature_);

    glutMainLoop ();
    return EXIT_SUCCESS;
}

