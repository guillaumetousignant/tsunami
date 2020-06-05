#ifndef MESHGEOMETRYUNSTRUCTURED_T_H
#define MESHGEOMETRYUNSTRUCTURED_T_H

#include <string>
#include <vector>
#include <another_path_tracer/entities/Vec3f.h>

class MeshGeometryUnstructured_t{
    public:
        MeshGeometryUnstructured_t(const std::string &filename);
        ~MeshGeometryUnstructured_t();

        unsigned int n_points_;
        unsigned int n_elements_;
        unsigned int n_normals_;
        unsigned int n_walls_;
        unsigned int* n_wall_;
        unsigned int sum_n_wall_;
        unsigned int n_farfield_;
        APTracer::Entities::Vec3f* points_;
        unsigned int* elements_;
        unsigned int* element_normals_;
        APTracer::Entities::Vec3f* normals_;
        unsigned int** walls_;
        unsigned int* farfield_;
        std::vector<std::vector<unsigned int>> point_to_elements_;

        void computeNormals(unsigned int n_points);
        void verify();
    
    private:
        void readSU2(const std::string &filename);
        void computeNodeToFace();
};

#endif