#ifndef MESH_T_H
#define MESH_T_H

#include "entities/Vec3f.h"

namespace APTracer { namespace Entities {
    class Material_t;
    class TransformMatrix_t;
    class Shape_t;
}}

class MeshGeometryUnstructured_t;

using APTracer::Entities::Material_t;
using APTracer::Entities::TransformMatrix_t;
using APTracer::Entities::Vec3f;
using APTracer::Entities::Shape_t;

class MeshUnstructured_t {
    public:
        MeshUnstructured_t(Material_t *material, TransformMatrix_t *transform_matrix, MeshGeometryUnstructured_t* geom);
        virtual ~MeshUnstructured_t() final;

        Material_t *material_;
        TransformMatrix_t *transformation_;
        MeshGeometryUnstructured_t* geom_;
        unsigned int n_tris_;
        Shape_t** triangles_; // Maybe should be triangle**?

        virtual void update();
        virtual Vec3f mincoord() const;
        virtual Vec3f maxcoord() const;
        virtual void createTriangles() final;
};
#endif