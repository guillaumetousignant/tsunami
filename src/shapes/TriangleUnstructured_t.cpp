#include "shapes/TriangleUnstructured_t.h"
#include <entities/TransformMatrix_t.h>
#include <entities/Material_t.h>
#include <cmath>
#include <limits>
#include <algorithm>
#include "entities/MeshGeometryUnstructured_t.h"

#define PI 3.141592653589793238463
#define EPSILON 0.00000001

using APTracer::Entities::Vec3f;

TriangleUnstructured_t::TriangleUnstructured_t(APTracer::Entities::Material_t *material, APTracer::Entities::TransformMatrix_t *transform_matrix, MeshGeometryUnstructured_t* geom, unsigned int index) 
    : Shape_t(material, transform_matrix), geom_(geom), index_(index) {

    const APTracer::Entities::TransformMatrix_t transform_norm = transformation_->transformDir();

    points_[0] = transformation_->multVec(geom_->points_[geom_->elements_[3 * index_]]);
    points_[1] = transformation_->multVec(geom_->points_[geom_->elements_[3 * index_ + 1]]);
    points_[2] = transformation_->multVec(geom_->points_[geom_->elements_[3 * index_ + 2]]);
    normals_[0] = transform_norm.multDir(geom_->normals_[geom_->element_normals_[3 * index_]]);
    normals_[1] = transform_norm.multDir(geom_->normals_[geom_->element_normals_[3 * index_ + 1]]);
    normals_[2] = transform_norm.multDir(geom_->normals_[geom_->element_normals_[3 * index_ + 2]]);

    v0v1_ = points_[1] - points_[0];
    v0v2_ = points_[2] - points_[0];
}

TriangleUnstructured_t::~TriangleUnstructured_t(){}

void TriangleUnstructured_t::update() {
    const APTracer::Entities::TransformMatrix_t transform_norm = transformation_->transformDir();

    points_[0] = transformation_->multVec(geom_->points_[geom_->elements_[3 * index_]]);
    points_[1] = transformation_->multVec(geom_->points_[geom_->elements_[3 * index_ + 1]]);
    points_[2] = transformation_->multVec(geom_->points_[geom_->elements_[3 * index_ + 2]]);
    normals_[0] = transform_norm.multDir(geom_->normals_[geom_->element_normals_[3 * index_]]);
    normals_[1] = transform_norm.multDir(geom_->normals_[geom_->element_normals_[3 * index_ + 1]]);
    normals_[2] = transform_norm.multDir(geom_->normals_[geom_->element_normals_[3 * index_ + 2]]);

    v0v1_ = points_[1] - points_[0];
    v0v2_ = points_[2] - points_[0];
}

bool TriangleUnstructured_t::intersection(const Ray_t &ray, double &t, std::array<double, 2> &uv) const {
    const Vec3f pvec = ray.direction_.cross(v0v2_);
    const double det = v0v1_.dot(pvec);

    if (std::abs(det) < EPSILON){
        t = std::numeric_limits<double>::infinity();
        uv[0] = NAN;
        uv[1] = NAN;
        return false;
    }

    const double invdet = 1.0/det;
    const Vec3f tvec = ray.origin_ - points_[0];
    const double u = tvec.dot(pvec) * invdet;
    uv[0] = u;

    if ((u < 0.0) || (u > 1.0)){
        t = std::numeric_limits<double>::infinity();
        uv[1] = NAN;
        return false;
    }

    const Vec3f qvec = tvec.cross(v0v1_);
    const double v = ray.direction_.dot(qvec) * invdet;
    uv[1] = v;

    if ((v < 0.0) || ((u+v) > 1.0)){
        t = std::numeric_limits<double>::infinity();
        return false;
    }

    t = v0v2_.dot(qvec) * invdet;

    if (t < 0.0){
        t = std::numeric_limits<double>::infinity();
        return false;
    }

    return true;
}

Vec3f TriangleUnstructured_t::normaluv(double time, std::array<double, 2> uv, std::array<double, 2> &tuv) const {
    const Vec3f distance = Vec3f(1.0 - uv[0] - uv[1], uv[0], uv[1]);

    // Matrix multiplication, optimise.
    tuv[0] = 0.0; // Not used
    tuv[1] = 0.0; // Not used

    return Vec3f(distance[0] * normals_[0][0] + distance[1] * normals_[1][0] + distance[2] * normals_[2][0], 
        distance[0] * normals_[0][1] + distance[1] * normals_[1][1] + distance[2] * normals_[2][1],
        distance[0] * normals_[0][2] + distance[1] * normals_[1][2] + distance[2] * normals_[2][2]);
}

Vec3f TriangleUnstructured_t::normal(double time, std::array<double, 2> uv) const {
    const Vec3f distance = Vec3f(1.0 - uv[0] - uv[1], uv[0], uv[1]);
    return Vec3f(distance[0] * normals_[0][0] + distance[1] * normals_[1][0] + distance[2] * normals_[2][0], 
        distance[0] * normals_[0][1] + distance[1] * normals_[1][1] + distance[2] * normals_[2][1],
        distance[0] * normals_[0][2] + distance[1] * normals_[1][2] + distance[2] * normals_[2][2]);
    // Matrix multiplication, optimise.
}

Vec3f TriangleUnstructured_t::normal_uv_tangent(double time, std::array<double, 2> uv, std::array<double, 2> &tuv, Vec3f &tangentvec) const {
    const Vec3f distance = Vec3f(1.0 - uv[0] - uv[1], uv[0], uv[1]);
    
    tuv[0] = 0.0; // Not used
    tuv[1] = 0.0; // Not used

    tangentvec = Vec3f(0.0, 0.0, 1.0);
    // Matrix multiplication, optimise.
    return Vec3f(distance[0] * normals_[0][0] + distance[1] * normals_[1][0] + distance[2] * normals_[2][0], 
        distance[0] * normals_[0][1] + distance[1] * normals_[1][1] + distance[2] * normals_[2][1],
        distance[0] * normals_[0][2] + distance[1] * normals_[1][2] + distance[2] * normals_[2][2]);
}  

Vec3f TriangleUnstructured_t::normal_face(double time) const{
    return v0v1_.cross(v0v2_).normalize_inplace();
}

Vec3f TriangleUnstructured_t::mincoord() const {
    return points_[0].getMin(points_[1]).min(points_[2]);
}

Vec3f TriangleUnstructured_t::maxcoord() const {
    return points_[0].getMax(points_[1]).max(points_[2]);
}