#include "shapes/MeshUnstructured_t.h"
#include "entities/MeshGeometryUnstructured_t.h"
#include "shapes/TriangleUnstructured_t.h"
#include <iostream>
#include <limits>

MeshUnstructured_t::MeshUnstructured_t(APTracer::Entities::Material_t *material, APTracer::Entities::TransformMatrix_t *transform_matrix, MeshGeometryUnstructured_t* geom) 
    : material_(material), transformation_(transform_matrix), geom_(geom), n_tris_(geom->n_elements_) {
    
    createTriangles();
}       

MeshUnstructured_t::~MeshUnstructured_t(){
    if (triangles_ != nullptr){
        for (unsigned int i = 0; i < n_tris_; i++){
            if (triangles_[i] != nullptr){
                delete triangles_[i];
            }
        }
        delete [] triangles_;
    }
}

void MeshUnstructured_t::createTriangles(){    
    triangles_ = new APTracer::Entities::Shape_t*[n_tris_];
    for (unsigned int i = 0; i < n_tris_; i++){
        triangles_[i] = new TriangleUnstructured_t(material_, transformation_, geom_, i);
    }
}

void MeshUnstructured_t::update(){
    for (unsigned int i = 0; i < n_tris_; i++){
        triangles_[i]->update();
    }
}

Vec3f MeshUnstructured_t::mincoord() const {
    Vec3f coord = Vec3f(std::numeric_limits<double>::infinity());
    for (unsigned int i = 0; i < n_tris_; i++){
        coord.min(triangles_[i]->mincoord());
    }    
    return coord;
}

Vec3f MeshUnstructured_t::maxcoord() const {
    Vec3f coord = Vec3f(-std::numeric_limits<double>::infinity());
    for (unsigned int i = 0; i < n_tris_; i++){
        coord.max(triangles_[i]->maxcoord());
    }
    return coord;
}