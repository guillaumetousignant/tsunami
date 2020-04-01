#include "another_path_tracer.h"
#include "entities/MeshGeometryUnstructured_t.h"
#include "shapes/MeshUnstructured_t.h"

#include <string>

using APTracer::Entities::Vec3f;

int main(int argc, char **argv){
    if (argc < 3) {
        std::cerr << "Error: input arguments with format 'mesh_file data_file'. Exiting." << std::endl;
    }
    std::string mesh_file = argv[1];
    std::string data_file = argv[2];

    MeshGeometryUnstructured_t water_mesh_geometry(mesh_file);
    MeshGeometryUnstructured_t sand_mesh_geometry(mesh_file);

    APTracer::Materials::Absorber_t water_scatterer(Vec3f(0.0, 0.0, 0.0), Vec3f(0.0, 0.0, 0.0), 1000, 1000);

    APTracer::Materials::ReflectiveRefractive_t water(Vec3f(0.0, 0.0, 0.0), Vec3f(1.0, 1.0, 1.0), 1.4, 10, &water_scatterer);
    APTracer::Materials::Diffuse_t sand(Vec3f(0.0, 0.0, 0.0), Vec3f(1.0, 1.0, 1.0), 1.0);

    MeshUnstructured_t water_mesh(&water, new TransformMatrix_t(), &water_mesh_geometry);
    MeshUnstructured_t sand_mesh(&sand, new TransformMatrix_t(), &sand_mesh_geometry);

    

    return 0;
}