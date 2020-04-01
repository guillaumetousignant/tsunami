#include "another_path_tracer.h"
#include "entities/MeshGeometryUnstructured_t.h"
#include "shapes/MeshUnstructured_t.h"
#include "entities/RandomGenerator_t.h"

#include <limits>
#include <cmath>
#include <string>

using APTracer::Entities::Vec3f;

double get_max_depth(MeshGeometryUnstructured_t* mesh_geometry);
void extrude_farfield(MeshGeometryUnstructured_t* mesh_geometry, double height);

int main(int argc, char **argv){
    if (argc < 3) {
        std::cerr << "Error: input arguments with format 'mesh_file data_file'. Exiting." << std::endl;
    }
    std::string mesh_file = argv[1];
    std::string data_file = argv[2];

    MeshGeometryUnstructured_t water_mesh_geometry(mesh_file);
    MeshGeometryUnstructured_t sand_mesh_geometry(mesh_file);

    unsigned int n_grid_points = water_mesh_geometry.n_points_;

    // Setting the water height to 0 for now
    for (unsigned int i = 0; i < n_grid_points; ++i) {
        water_mesh_geometry.points_[i][2] = 0.0;
    }

    // Negating the sand height
    for (unsigned int i = 0; i < n_grid_points; ++i) {
        sand_mesh_geometry.points_[i][2] = -sand_mesh_geometry.points_[i][2];
    }

    double max_depth = get_max_depth(&sand_mesh_geometry); // Is negative
    extrude_farfield(&sand_mesh_geometry, 2 * max_depth);

    // Render stuff
    APTracer::Materials::Absorber_t water_scatterer(Vec3f(0.0, 0.0, 0.0), Vec3f(0.0, 0.0, 0.0), 1000, 1000);
    APTracer::Materials::NonAbsorber_t air_scatterer;

    APTracer::Materials::Refractive_t water(Vec3f(0.0, 0.0, 0.0), Vec3f(1.0, 1.0, 1.0), 1.33, 10, &water_scatterer);
    APTracer::Materials::Diffuse_t sand(Vec3f(0.0, 0.0, 0.0), Vec3f(1.0, 0.9217, 0.7098), 1.0);
    APTracer::Materials::Transparent_t air(0, &air_scatterer);

    APTracer::Entities::TransformMatrix_t water_transform;
    APTracer::Entities::TransformMatrix_t sand_transform;

    MeshUnstructured_t water_mesh(&water, &water_transform, &water_mesh_geometry);
    MeshUnstructured_t sand_mesh(&sand, &sand_transform, &sand_mesh_geometry);

    APTracer::Entities::TransformMatrix_t sun_transform;

    APTracer::Entities::DirectionalLight_t sun(Vec3f(10.0, 10.0, 8.0), &sun_transform);
    sun.transformation_->scale(0.95);
    sun.transformation_->rotateZ(-0.7854);
    sun.transformation_->rotateX(-1.1781);
    sun.update();

    APTracer::Skyboxes::SkyboxFlatSun_t sky(Vec3f(0.9020, 0.9725, 1.0), &sun);

    APTracer::Entities::Scene_t scene;
    scene.add(water_mesh.triangles_, water_mesh.n_tris_);
    scene.add(sand_mesh.triangles_, sand_mesh.n_tris_);

    APTracer::Entities::ImgBufferOpenGL_t imgbuffer(1800, 1200);

    APTracer::Entities::TransformMatrix_t camera_transform;
    double fov[2] = {2.0/3.0 * 80.0 * M_PI/180.0, 80.0 * M_PI/180.0};
    unsigned int subpix[2] = {1, 1};
    std::list<Medium_t*> medium_list = {&air, &air};
    APTracer::Cameras::Cam_t camera(&camera_transform, "images/output.png", Vec3f(0.0, 0.0, 1.0), fov, subpix, &imgbuffer, medium_list, &sky, 8, 1.0);
    camera.transformation_->translate(Vec3f(0.0, -2000.0, 0.0));
    camera.transformation_->rotateXAxis(-30.0 * M_PI/180);
    camera.update();

    scene.build_acc();

    APTracer::Entities::OpenGLRenderer_t opengl_renderer(&scene, &camera, &imgbuffer);
    opengl_renderer.initialise();
    opengl_renderer.render();

    return 0;
}

double get_max_depth(MeshGeometryUnstructured_t* mesh_geometry) {
    double depth = std::numeric_limits<double>::infinity();
    for (unsigned int i = 0; i < mesh_geometry->n_points_; ++i) {
        depth = std::min(depth, mesh_geometry->points_[i][2]);
    }
    return depth;
}

void extrude_farfield(MeshGeometryUnstructured_t* mesh_geometry, double height) {
    unsigned int new_n_points = mesh_geometry->n_points_ + mesh_geometry->n_farfield_;
    unsigned int new_n_elements = mesh_geometry->n_elements_ + mesh_geometry->n_farfield_ * 2;
    unsigned int new_n_normals = mesh_geometry->n_normals_ + mesh_geometry->n_farfield_;

    // Putting back old stuff
    APTracer::Entities::Vec3f* new_points = new Vec3f[new_n_points];
    for (unsigned int i = 0; i < mesh_geometry->n_points_; ++i) {
        new_points[i] = mesh_geometry->points_[i];
    }
    APTracer::Entities::Vec3f* new_normals = new Vec3f[new_n_normals];
    for (unsigned int i = 0; i < mesh_geometry->n_normals_; ++i) {
        new_normals[i] = mesh_geometry->normals_[i];
    }
    unsigned int* new_elements = new unsigned int[3 * new_n_elements];
    for (unsigned int i = 0; i < mesh_geometry->n_elements_; ++i) {
        new_elements[3 * i] = mesh_geometry->elements_[3 * i];
        new_elements[3 * i + 1] = mesh_geometry->elements_[3 * i + 1];
        new_elements[3 * i + 2] = mesh_geometry->elements_[3 * i + 2];
    }
    unsigned int* new_element_normals = new unsigned int[3 * new_n_elements];
    for (unsigned int i = 0; i < mesh_geometry->n_elements_; ++i) {
        new_element_normals[3 * i] = mesh_geometry->element_normals_[3 * i];
        new_element_normals[3 * i + 1] = mesh_geometry->element_normals_[3 * i + 1];
        new_element_normals[3 * i + 2] = mesh_geometry->element_normals_[3 * i + 2];
    }

    for (unsigned int i = 0; i < mesh_geometry->n_farfield_; ++i){
        new_points[i + mesh_geometry->n_points_] = mesh_geometry->points_[mesh_geometry->farfield_[2 * i]] + Vec3f(0.0, 0.0, height);
    }

    for (unsigned int i = 0; i < mesh_geometry->n_farfield_; ++i){
        new_normals[i + mesh_geometry->n_normals_] = mesh_geometry->points_[mesh_geometry->farfield_[2 * i]].normalize();
    }

    // Adds two elements per boundary, created with the new points
    for (unsigned int i = 0; i < mesh_geometry->n_farfield_ - 1; ++i){
        new_elements[3 * mesh_geometry->n_elements_ + 6 * i] = mesh_geometry->farfield_[2 * i];
        new_elements[3 * mesh_geometry->n_elements_ + 6 * i + 1] = mesh_geometry->n_points_ + i;
        new_elements[3 * mesh_geometry->n_elements_ + 6 * i + 2] = mesh_geometry->n_points_ + i + 1;
        new_elements[3 * mesh_geometry->n_elements_ + 6 * i + 3] = mesh_geometry->farfield_[2 * i + 1];
        new_elements[3 * mesh_geometry->n_elements_ + 6 * i + 4] = mesh_geometry->farfield_[2 * i];
        new_elements[3 * mesh_geometry->n_elements_ + 6 * i + 5] = mesh_geometry->n_points_ + i + 1;
    }
    new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_ - 1)] = mesh_geometry->farfield_[2 * (mesh_geometry->n_farfield_ - 1)];
    new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_ - 1) + 1] = mesh_geometry->n_points_ + (mesh_geometry->n_farfield_ - 1);
    new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_ - 1) + 2] = mesh_geometry->n_points_;
    new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_ - 1) + 3] = mesh_geometry->farfield_[2 * (mesh_geometry->n_farfield_ - 1) + 1];
    new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_ - 1) + 4] = mesh_geometry->farfield_[2 * (mesh_geometry->n_farfield_ - 1)];
    new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_ - 1) + 5] = mesh_geometry->n_points_ ;

    for (unsigned int i = 0; i < mesh_geometry->n_farfield_ - 1; ++i){
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * i] = mesh_geometry->n_normals_ + i;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * i + 1] = mesh_geometry->n_normals_ + i;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * i + 2] = mesh_geometry->n_normals_ + i + 1;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * i + 3] = mesh_geometry->n_normals_ + i;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * i + 4] = mesh_geometry->n_normals_ + i + 1;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * i + 5] = mesh_geometry->n_normals_ + i + 1;
    }
    new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_ - 1)] = mesh_geometry->n_normals_ + (mesh_geometry->n_farfield_ - 1);
    new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_ - 1) + 1] = mesh_geometry->n_normals_ + (mesh_geometry->n_farfield_ - 1);
    new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_ - 1) + 2] = mesh_geometry->n_normals_;
    new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_ - 1) + 3] = mesh_geometry->n_normals_ + (mesh_geometry->n_farfield_ - 1);
    new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_ - 1) + 4] = mesh_geometry->n_normals_;
    new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_ - 1) + 5] = mesh_geometry->n_normals_ ;

    std::swap(mesh_geometry->points_, new_points);
    std::swap(mesh_geometry->normals_, new_normals);
    std::swap(mesh_geometry->elements_, new_elements);
    std::swap(mesh_geometry->element_normals_, new_element_normals);
    mesh_geometry->n_points_ = new_n_points;
    mesh_geometry->n_normals_ = new_n_normals;
    mesh_geometry->n_elements_ = new_n_elements;

    delete [] new_points;
    delete [] new_normals;
    delete [] new_elements;
    delete [] new_element_normals;
}