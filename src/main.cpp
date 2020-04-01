#include "another_path_tracer.h"
#include "entities/MeshGeometryUnstructured_t.h"
#include "shapes/MeshUnstructured_t.h"
#include "entities/RandomGenerator_t.h"

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

    // Setting the water height to 0 for now
    for (unsigned int i = 0; i < water_mesh_geometry.n_points_; ++i) {
        water_mesh_geometry.points_[i][2] = 0.0;
    }

    // Negating the sand height
    for (unsigned int i = 0; i < sand_mesh_geometry.n_points_; ++i) {
        sand_mesh_geometry.points_[i][2] = -sand_mesh_geometry.points_[i][2];
    }

    APTracer::Materials::Absorber_t water_scatterer(Vec3f(0.0, 0.0, 0.0), Vec3f(0.0, 0.0, 0.0), 1000, 1000);

    APTracer::Materials::ReflectiveRefractive_t water(Vec3f(0.0, 0.0, 0.0), Vec3f(1.0, 1.0, 1.0), 1.4, 10, &water_scatterer);
    APTracer::Materials::Diffuse_t sand(Vec3f(0.0, 0.0, 0.0), Vec3f(1.0, 1.0, 1.0), 1.0);

    APTracer::Entities::TransformMatrix_t water_transform;
    APTracer::Entities::TransformMatrix_t sand_transform;

    MeshUnstructured_t water_mesh(&water, &water_transform, &water_mesh_geometry);
    MeshUnstructured_t sand_mesh(&sand, &sand_transform, &sand_mesh_geometry);

    APTracer::Entities::TransformMatrix_t sun_transform;

    APTracer::Entities::DirectionalLight_t sun(Vec3f(0.0, 0.0, 0.0), &sun_transform);
    sun.transformation_->scale(0.95);
    sun.transformation_->rotateZ(-0.7854);
    sun.transformation_->rotateX(-1.1781);
    sun.update();

    APTracer::Skyboxes::SkyboxFlatSun_t sky(Vec3f(0.7176, 0.9569, 1.0), &sun);

    APTracer::Entities::Scene_t scene;
    scene.add(water_mesh.triangles_, water_mesh.n_tris_);
    scene.add(sand_mesh.triangles_, sand_mesh.n_tris_);

    APTracer::Entities::TransformMatrix_t camera_transform;
    APTracer::Cameras::Cam_t camera(&camera_transform, );

    return 0;
}