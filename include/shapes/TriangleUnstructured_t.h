#ifndef TRIANGLEUNSTRUCTURED_T_H
#define TRIANGLEUNSTRUCTURED_T_H

#include <entities/Ray_t.h>
#include <entities/Vec3f.h>
#include <entities/Shape_t.h>

namespace APTracer { namespace Entities {
    class TransformMatrix_t;
    class Material_t;
}}

class MeshGeometryUnstructured_t;

using APTracer::Entities::Ray_t;
using APTracer::Entities::Vec3f;
using APTracer::Entities::Material_t;
using APTracer::Entities::TransformMatrix_t;
using APTracer::Entities::Shape_t;

class TriangleUnstructured_t final : public Shape_t{
    public:
        TriangleUnstructured_t(Material_t *material, TransformMatrix_t *transform_matrix, MeshGeometryUnstructured_t* geom, unsigned int index);
        virtual ~TriangleUnstructured_t() final;

        Vec3f points_[3];
        Vec3f normals_[3];
        Vec3f v0v1_;
        Vec3f v0v2_;
        MeshGeometryUnstructured_t* geom_;
        unsigned int index_;

        virtual void update() final;
        virtual bool intersection(const Ray_t &ray, double &t, std::array<double, 2> &uv) const final; 
        virtual Vec3f normaluv(double time, std::array<double, 2> uv, std::array<double, 2> &tuv) const final;
        virtual Vec3f normal(double time, std::array<double, 2> uv) const final;
        virtual Vec3f normal_uv_tangent(double time, std::array<double, 2> uv, std::array<double, 2> &tuv, Vec3f &tangentvec) const final;
        virtual Vec3f normal_face(double time) const final;
        virtual Vec3f mincoord() const final;
        virtual Vec3f maxcoord() const final;
};
#endif