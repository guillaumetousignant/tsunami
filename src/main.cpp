#include "another_path_tracer.h"
#include "entities/MeshGeometryUnstructured_t.h"
#include "shapes/MeshUnstructured_t.h"
#include "entities/RandomGenerator_t.h"

#include <limits>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <complex>
#include <iomanip>
#include <chrono>

using APTracer::Entities::Vec3f;

double get_max_depth(MeshGeometryUnstructured_t* mesh_geometry);
void extrude_farfield(MeshGeometryUnstructured_t* mesh_geometry, double height);
void extrude_wall(MeshGeometryUnstructured_t* mesh_geometry, double height);
std::vector<std::complex<double>> get_eta(std::string filename, double &amplitude, double &omega);
void openGL_accumulate();

namespace Rendering {
    APTracer::Entities::OpenGLRenderer_t* renderer = nullptr;
    double time = 0.0;
    double delta_time = 1;
    MeshGeometryUnstructured_t* mesh_geometry = nullptr;
    MeshUnstructured_t* mesh = nullptr;
    APTracer::Entities::AccelerationStructure_t* acc = nullptr;
    std::vector<std::vector<std::complex<double>>> etas;
    unsigned int n_points = 0;
    std::vector<double> omegas;
    unsigned int n_timestep = 0;
    unsigned int write_interval = 1000;
}


int main(int argc, char **argv){
    if (argc < 3) {
        std::cerr << "Error: input arguments with format 'mesh_file data_file'. Exiting." << std::endl;
    }
    std::string mesh_file = argv[1];

    unsigned int n_waves = argc - 2;
    std::vector<std::string> data_files(n_waves, "");
    for (unsigned int i = 0; i < n_waves; ++i) {
        data_files[i] = argv[2 + i];
    }

    std::vector<double> amplitudes(n_waves, 0);
    std::vector<double> omegas(n_waves, 0);
    std::vector<std::vector<std::complex<double>>> etas(n_waves, std::vector<std::complex<double>>());
    for (unsigned int i = 0; i < n_waves; ++i) {
        etas[i] = get_eta(data_files[i], amplitudes[i], omegas[i]);
    }

    MeshGeometryUnstructured_t water_mesh_geometry(mesh_file);
    MeshGeometryUnstructured_t sand_mesh_geometry(mesh_file);

    unsigned int n_grid_points = water_mesh_geometry.n_points_;

    // Setting the water height at t = 0
    for (unsigned int i = 0; i < n_grid_points; ++i) {
        water_mesh_geometry.points_[i][2] = 0.0;
        for (unsigned int j = 0; j < etas.size(); ++j) {
            water_mesh_geometry.points_[i][2] += std::real(etas[j][i] * std::exp(std::complex<double>(0.0, -1.0) * omegas[j] * 0.0));
        }
    }

    water_mesh_geometry.computeNormals(n_grid_points);

    // Negating the sand height
    for (unsigned int i = 0; i < n_grid_points; ++i) {
        sand_mesh_geometry.points_[i][2] = -sand_mesh_geometry.points_[i][2];
    }
    sand_mesh_geometry.computeNormals(n_grid_points);

    double max_depth = get_max_depth(&sand_mesh_geometry); // Is negative
    extrude_farfield(&sand_mesh_geometry, 4 * max_depth);
    extrude_farfield(&water_mesh_geometry, max_depth);
    extrude_wall(&sand_mesh_geometry, -max_depth);

    // Render stuff
    APTracer::Materials::Absorber_t water_scatterer(Vec3f(0.0, 0.0, 0.0), Vec3f(0.92, 0.97, 0.99), 1000, 32);
    APTracer::Materials::NonAbsorber_t air_scatterer;

    APTracer::Materials::ReflectiveRefractive_t water(Vec3f(0.0, 0.0, 0.0), Vec3f(1.0, 1.0, 1.0), 1.33, 10, &water_scatterer);
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

    APTracer::Entities::ImgBufferOpenGL_t imgbuffer(1920, 1080);

    APTracer::Entities::TransformMatrix_t camera_transform;
    double fov[2] = {9.0/16.0 * 80.0 * M_PI/180.0, 80.0 * M_PI/180.0};
    unsigned int subpix[2] = {1, 1};
    std::list<Medium_t*> medium_list = {&air, &air};
    APTracer::Cameras::Cam_t camera(&camera_transform, "images/output.png", Vec3f(0.0, 0.0, 1.0), fov, subpix, &imgbuffer, medium_list, &sky, 16, 1.0);
    camera.transformation_->translate(Vec3f(0.0, -12000.0, 0.0));
    camera.transformation_->rotateXAxis(-30.0 * M_PI/180);
    camera.transformation_->translate(Vec3f(0.0, 0.0, -2500.0));
    camera.update();

    scene.build_acc();

    APTracer::Entities::OpenGLRenderer_t opengl_renderer(&scene, &camera, &imgbuffer);
    Rendering::renderer = &opengl_renderer;
    // We need to overload the render function with our own to apply movement etc
    opengl_renderer.render_function_ = openGL_accumulate;
    Rendering::mesh_geometry = &water_mesh_geometry;
    Rendering::mesh = &water_mesh;
    Rendering::acc = scene.acc_;
    Rendering::etas = etas;
    Rendering::n_points = n_grid_points;
    Rendering::omegas = omegas;

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

    // For closing stuff
    new_n_points += 1;
    new_n_elements += mesh_geometry->n_farfield_;
    new_n_normals += 1;

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

    // New stuff
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
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * i + 3] = mesh_geometry->n_normals_ + i + 1;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * i + 4] = mesh_geometry->n_normals_ + i;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * i + 5] = mesh_geometry->n_normals_ + i + 1;
    }
    new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_ - 1)] = mesh_geometry->n_normals_ + (mesh_geometry->n_farfield_ - 1);
    new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_ - 1) + 1] = mesh_geometry->n_normals_ + (mesh_geometry->n_farfield_ - 1);
    new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_ - 1) + 2] = mesh_geometry->n_normals_;
    new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_ - 1) + 3] = mesh_geometry->n_normals_;
    new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_ - 1) + 4] = mesh_geometry->n_normals_ + (mesh_geometry->n_farfield_ - 1);
    new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_ - 1) + 5] = mesh_geometry->n_normals_;

    // Closing
    // Making new stuff
    new_points[mesh_geometry->n_points_ + mesh_geometry->n_farfield_] = Vec3f(0.0, 0.0, height);
    new_normals[mesh_geometry->n_points_ + mesh_geometry->n_farfield_] = Vec3f(0.0, 0.0, 1.0);

    // Adds one element per boundary, created with the new points
    for (unsigned int i = 0; i < mesh_geometry->n_farfield_ - 1; ++i){
        new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_) + 3 * i] = mesh_geometry->n_points_ + i;
        new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_) + 3 * i + 1] = mesh_geometry->n_points_ + i + 1;
        new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_) + 3 * i + 2] = mesh_geometry->n_points_ + mesh_geometry->n_farfield_;
    }
    new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_) + 3 * (mesh_geometry->n_farfield_ - 1)] = mesh_geometry->n_points_ + (mesh_geometry->n_farfield_ - 1);
    new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_) + 3 * (mesh_geometry->n_farfield_ - 1) + 1] = mesh_geometry->n_points_;
    new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_) + 3 * (mesh_geometry->n_farfield_ - 1) + 2] = mesh_geometry->n_points_ + mesh_geometry->n_farfield_;

    for (unsigned int i = 0; i < mesh_geometry->n_farfield_; ++i){
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_) + 3 * i] = mesh_geometry->n_points_ + mesh_geometry->n_farfield_;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_) + 3 * i + 1] = mesh_geometry->n_points_ + mesh_geometry->n_farfield_;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_farfield_) + 3 * i + 2] = mesh_geometry->n_points_ + mesh_geometry->n_farfield_;
    }

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

void extrude_wall(MeshGeometryUnstructured_t* mesh_geometry, double height) {
    unsigned int new_n_points = mesh_geometry->n_points_ + mesh_geometry->n_wall_;
    unsigned int new_n_elements = mesh_geometry->n_elements_ + mesh_geometry->n_wall_ * 2;
    unsigned int new_n_normals = mesh_geometry->n_normals_ + mesh_geometry->n_wall_;

    // For closing stuff
    new_n_points += 1;
    new_n_elements += mesh_geometry->n_wall_;
    new_n_normals += 1;

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

    // Making new stuff
    for (unsigned int i = 0; i < mesh_geometry->n_wall_; ++i){
        new_points[i + mesh_geometry->n_points_] = mesh_geometry->points_[mesh_geometry->wall_[2 * i]] + Vec3f(0.0, 0.0, height);
    }

    for (unsigned int i = 0; i < mesh_geometry->n_wall_; ++i){
        new_normals[i + mesh_geometry->n_normals_] = mesh_geometry->points_[mesh_geometry->wall_[2 * i]].normalize();
    }

    // Adds two elements per boundary, created with the new points
    for (unsigned int i = 0; i < mesh_geometry->n_wall_ - 1; ++i){
        new_elements[3 * mesh_geometry->n_elements_ + 6 * i] = mesh_geometry->wall_[2 * i];
        new_elements[3 * mesh_geometry->n_elements_ + 6 * i + 1] = mesh_geometry->n_points_ + i;
        new_elements[3 * mesh_geometry->n_elements_ + 6 * i + 2] = mesh_geometry->n_points_ + i + 1;
        new_elements[3 * mesh_geometry->n_elements_ + 6 * i + 3] = mesh_geometry->wall_[2 * i + 1];
        new_elements[3 * mesh_geometry->n_elements_ + 6 * i + 4] = mesh_geometry->wall_[2 * i];
        new_elements[3 * mesh_geometry->n_elements_ + 6 * i + 5] = mesh_geometry->n_points_ + i + 1;
    }
    new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_ - 1)] = mesh_geometry->wall_[2 * (mesh_geometry->n_wall_ - 1)];
    new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_ - 1) + 1] = mesh_geometry->n_points_ + (mesh_geometry->n_wall_ - 1);
    new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_ - 1) + 2] = mesh_geometry->n_points_;
    new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_ - 1) + 3] = mesh_geometry->wall_[2 * (mesh_geometry->n_wall_ - 1) + 1];
    new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_ - 1) + 4] = mesh_geometry->wall_[2 * (mesh_geometry->n_wall_ - 1)];
    new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_ - 1) + 5] = mesh_geometry->n_points_ ;

    for (unsigned int i = 0; i < mesh_geometry->n_wall_ - 1; ++i){
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * i] = mesh_geometry->n_normals_ + i;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * i + 1] = mesh_geometry->n_normals_ + i;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * i + 2] = mesh_geometry->n_normals_ + i + 1;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * i + 3] = mesh_geometry->n_normals_ + i + 1;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * i + 4] = mesh_geometry->n_normals_ + i;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * i + 5] = mesh_geometry->n_normals_ + i + 1;
    }
    new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_ - 1)] = mesh_geometry->n_normals_ + (mesh_geometry->n_wall_ - 1);
    new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_ - 1) + 1] = mesh_geometry->n_normals_ + (mesh_geometry->n_wall_ - 1);
    new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_ - 1) + 2] = mesh_geometry->n_normals_;
    new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_ - 1) + 3] = mesh_geometry->n_normals_;
    new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_ - 1) + 4] = mesh_geometry->n_normals_ + (mesh_geometry->n_wall_ - 1);
    new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_ - 1) + 5] = mesh_geometry->n_normals_;

    // Closing
    // Making new stuff
    new_points[mesh_geometry->n_points_ + mesh_geometry->n_wall_] = Vec3f(0.0, 0.0, height);
    new_normals[mesh_geometry->n_points_ + mesh_geometry->n_wall_] = Vec3f(0.0, 0.0, 1.0);

    // Adds one element per boundary, created with the new points
    for (unsigned int i = 0; i < mesh_geometry->n_wall_ - 1; ++i){
        new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_) + 3 * i] = mesh_geometry->n_points_ + i;
        new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_) + 3 * i + 1] = mesh_geometry->n_points_ + i + 1;
        new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_) + 3 * i + 2] = mesh_geometry->n_points_ + mesh_geometry->n_wall_;
    }
    new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_) + 3 * (mesh_geometry->n_wall_ - 1)] = mesh_geometry->n_points_ + (mesh_geometry->n_wall_ - 1);
    new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_) + 3 * (mesh_geometry->n_wall_ - 1) + 1] = mesh_geometry->n_points_;
    new_elements[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_) + 3 * (mesh_geometry->n_wall_ - 1) + 2] = mesh_geometry->n_points_ + mesh_geometry->n_wall_;

    for (unsigned int i = 0; i < mesh_geometry->n_wall_; ++i){
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_) + 3 * i] = mesh_geometry->n_points_ + mesh_geometry->n_wall_;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_) + 3 * i + 1] = mesh_geometry->n_points_ + mesh_geometry->n_wall_;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * (mesh_geometry->n_wall_) + 3 * i + 2] = mesh_geometry->n_points_ + mesh_geometry->n_wall_;
    }

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

std::vector<std::complex<double>> get_eta(std::string filename, double &amplitude, double &omega) {
    std::string line;
    std::string token;
    double n_points;

    std::ifstream meshfile(filename);
    if (!meshfile.is_open()) {
        std::cerr << "Error: file '" << filename << "' could not be opened. Exiting." << std::endl;
        return std::vector<std::complex<double>>();
    }

    std::getline(meshfile, line);
    std::istringstream liness(line);
    liness >> token;
    if (token == "AMPLITUDE="){
        liness >> amplitude;
    }
    else {
        std::cerr << "Error: expected marker 'AMPLITUDE=', found '" << token << "'. Exiting." << std::endl;
        return std::vector<std::complex<double>>();
    }

    std::getline(meshfile, line);
    std::istringstream liness2(line);
    liness2 >> token;
    if (token == "OMEGA="){
        liness2 >> omega;
    }
    else {
        std::cerr << "Error: expected marker 'OMEGA=', found '" << token << "'. Exiting." << std::endl;
        return std::vector<std::complex<double>>();
    }

    std::getline(meshfile, line);
    std::getline(meshfile, line);
    std::istringstream liness3(line);
    liness3 >> token;
    if (token == "NPOIN="){
        liness3 >> n_points;
    }
    else {
        std::cerr << "Error: expected marker 'NPOIN=', found '" << token << "'. Exiting." << std::endl;
        return std::vector<std::complex<double>>();
    }

    std::vector<std::complex<double>> eta(n_points, std::complex<double>(0.0, 0.0));
    for (unsigned int i = 0; i < n_points; ++i) {
        std::getline(meshfile, line);
        std::istringstream liness4(line);
        double val0, val1;
        liness4 >> val0 >> val1;

        eta[i] = std::complex<double>(val0, val1);
    }

    meshfile.close();

    return eta;
}

void timestep(MeshGeometryUnstructured_t* mesh_geometry, MeshUnstructured_t* mesh, APTracer::Entities::AccelerationStructure_t* acc, std::vector<std::vector<std::complex<double>>> etas, unsigned int n_points, double time, std::vector<double> omegas) {
    for (unsigned int i = 0; i < n_points; ++i) {
        mesh_geometry->points_[i][2] = 0.0;
        for (unsigned int j = 0; j < etas.size(); ++j) {
            mesh_geometry->points_[i][2] += std::real(etas[j][i] * std::exp(std::complex<double>(0.0, -1.0) * omegas[j] * time));
        }
    }

    mesh_geometry->computeNormals(n_points);

    for (unsigned int i = 0; i < mesh_geometry->n_elements_; ++i) {
        acc->remove(mesh->triangles_[i]);
        mesh->triangles_[i]->update();
        acc->add(mesh->triangles_[i]);
    }
}

void openGL_accumulate() {
    auto t_start = std::chrono::high_resolution_clock::now();
    Rendering::renderer->accumulate();
    auto t_end = std::chrono::high_resolution_clock::now();

    std::cout << "Iteration " << Rendering::renderer->n_iter_gl_ << " done in " 
        << std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0 
        << "s." << std::endl;
    
    if (Rendering::renderer->n_iter_gl_ == Rendering::write_interval) {
        ++Rendering::n_timestep;
        std::cout << "Timestep " << Rendering::n_timestep << ", t = " << Rendering::time << std::endl;
        std::ostringstream oss;
        oss << "images/image_"<< std::setfill('0') << std::setw(4) << Rendering::n_timestep << ".png";
        Rendering::renderer->camera_->write(oss.str());
        Rendering::time += Rendering::delta_time;
        timestep(Rendering::mesh_geometry, Rendering::mesh, Rendering::acc, Rendering::etas, Rendering::n_points, Rendering::time, Rendering::omegas);
        Rendering::renderer->resetDisplay();
    }  
}
