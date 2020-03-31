#include "shapes/MeshUnstructured_t.h"
#include "entities/MeshGeometry_t.h"
#include "entities/MaterialMap_t.h"
#include "shapes/TriangleMesh_t.h"
#include <iostream>
#include <limits>

MeshUnstructured_t::MeshUnstructured_t(APTracer::Entities::Material_t *material, APTracer::Entities::TransformMatrix_t *transform_matrix, APTracer::Entities::MeshGeometry_t* geom) 
    : MeshTop_t(material, transform_matrix, geom) {
    
    createTriangles();
}

MeshUnstructured_t::MeshUnstructured_t(APTracer::Entities::MaterialMap_t *materialmap, APTracer::Entities::TransformMatrix_t *transform_matrix, APTracer::Entities::MeshGeometry_t* geom) 
    : MeshTop_t(materialmap->getFirst(), transform_matrix, geom) {

    createTriangles(materialmap);
}        

MeshUnstructured_t::~MeshUnstructured_t(){}

void MeshUnstructured_t::createTriangles(){    
    triangles_ = new APTracer::Entities::Shape_t*[n_tris_];
    for (unsigned int i = 0; i < n_tris_; i++){
        triangles_[i] = new APTracer::Shapes::TriangleMesh_t(material_, transformation_, geom_, i);
    }
}

void MeshUnstructured_t::createTriangles(APTracer::Entities::MaterialMap_t *materialmap){
    triangles_ = new APTracer::Entities::Shape_t*[n_tris_];
    for (unsigned int i = 0; i < n_tris_; i++){
        triangles_[i] = new APTracer::Shapes::TriangleMesh_t(materialmap->getMaterial(geom_->mat_[i]), transformation_, geom_, i);
    }
}