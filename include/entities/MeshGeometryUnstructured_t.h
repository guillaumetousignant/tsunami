#ifndef MESHGEOMETRYUNSTRUCTURED_T_H
#define MESHGEOMETRYUNSTRUCTURED_T_H

#include <string>
#include <vector>
#include "entities/Vec3f.h"

class MeshGeometryUnstructured_t{
    public:
        MeshGeometryUnstructured_t(const std::string &filename);
        ~MeshGeometryUnstructured_t();

        unsigned int n_points_;
        unsigned int n_elements_;
        unsigned int n_wall_;
        unsigned int n_farfield_;
        APTracer::Entities::Vec3f* points_;
        unsigned int* elements_;
        APTracer::Entities::Vec3f* normals_;
        unsigned int* wall_;
        unsigned int* farfield_;
        std::vector<std::vector<unsigned int>> point_to_elements_;
    
    private:
        void readSU2(const std::string &filename);
        void computeNodeToFace();
        void computeNormals();
};

#endif