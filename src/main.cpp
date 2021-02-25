#include <another_path_tracer.h>
#include "entities/MeshGeometryUnstructured_t.h"
#include "shapes/MeshUnstructured_t.h"
#include <entities/RandomGenerator_t.h>

#include <limits>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <complex>
#include <iomanip>
#include <chrono>
#include <array>

using APTracer::Entities::Vec3f;

double get_max_depth(MeshGeometryUnstructured_t* mesh_geometry);
void extrude_farfield(MeshGeometryUnstructured_t* mesh_geometry, double height, bool close);
void extrude_wall(MeshGeometryUnstructured_t* mesh_geometry, double height);
std::vector<std::complex<double>> get_eta(std::string filename, double &amplitude, double &omega);
void openGL_accumulate();

constexpr double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

namespace Rendering {
    APTracer::Entities::OpenGLRenderer_t* renderer = nullptr;
    unsigned int n_timestep = 0;
    double delta_time = 10.0;
    double time = n_timestep * delta_time;
    double angle_step = pi/720.0; // 0.25 deg
    MeshGeometryUnstructured_t* mesh_geometry = nullptr;
    MeshUnstructured_t* mesh = nullptr;
    APTracer::Entities::Scene_t* scene = nullptr;
    std::vector<std::vector<std::complex<double>>> etas;
    unsigned int n_points = 0;
    std::vector<double> omegas;
    unsigned int write_interval = 1000;
}


int main(int argc, char **argv){
    if (argc < 3) {
        std::cerr << "Error: input arguments with format 'mesh_file data_file'. Exiting." << std::endl;
        return(1);
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

    // Negating the sand height
    for (unsigned int i = 0; i < n_grid_points; ++i) {
        sand_mesh_geometry.points_[i][2] = -sand_mesh_geometry.points_[i][2];
    }
    sand_mesh_geometry.computeNormals(n_grid_points);

    double max_depth = get_max_depth(&sand_mesh_geometry); // Is negative

    extrude_farfield(&sand_mesh_geometry, 2 * max_depth, true);
    extrude_farfield(&water_mesh_geometry, 2 * max_depth, false);
    extrude_wall(&sand_mesh_geometry, -max_depth/4.0);

    // Setting the water height at t = 0
    for (unsigned int i = 0; i < n_grid_points; ++i) {
        water_mesh_geometry.points_[i][2] = 0.0;
        for (unsigned int j = 0; j < etas.size(); ++j) {
            water_mesh_geometry.points_[i][2] += std::real(etas[j][i] * std::exp(std::complex<double>(0.0, -1.0) * omegas[j] * Rendering::time));
        }
    }

    water_mesh_geometry.computeNormals(n_grid_points);

    // Render stuff
    APTracer::Materials::Absorber_t water_scatterer(Vec3f(0.0, 0.0, 0.0), Vec3f(0.92, 0.97, 0.99), 1000000, 2048, 1.33, 10);
    APTracer::Materials::NonAbsorber_t air(1.0, 0);

    APTracer::Materials::ReflectiveRefractive_t water(Vec3f(0.0, 0.0, 0.0), Vec3f(1.0, 1.0, 1.0), &water_scatterer);
    APTracer::Materials::Diffuse_t sand(Vec3f(0.0, 0.0, 0.0), Vec3f(1.0, 0.9217, 0.7098), 1.0);

    APTracer::Entities::TransformMatrix_t water_transform;
    APTracer::Entities::TransformMatrix_t sand_transform;

    MeshUnstructured_t water_mesh(&water, &water_transform, &water_mesh_geometry);
    MeshUnstructured_t sand_mesh(&sand, &sand_transform, &sand_mesh_geometry);

    water_mesh.transformation_->scale(Vec3f(1.0, 1.0, 8.0));
    water_mesh.update();
    sand_mesh.transformation_->scale(Vec3f(1.0, 1.0, 8.0));
    sand_mesh.update();

    APTracer::Entities::TransformMatrix_t sun_transform;

    APTracer::Entities::DirectionalLight_t sun(Vec3f(10.0, 10.0, 8.0), &sun_transform);
    sun.transformation_->scale(0.95);
    sun.transformation_->rotateZ(-0.7854);
    sun.transformation_->rotateX(-1.1781);
    sun.update();

    APTracer::Entities::Texture_t background("assets/background.jpg");
    //APTracer::Skyboxes::SkyboxFlatSun_t sky(Vec3f(0.9020, 0.9725, 1.0), &sun);
    const std::array<double, 2> sun_pos {0.62093, 0.77075};
    APTracer::Skyboxes::SkyboxTextureSun_t sky(&background, sun_pos, Vec3f(12.6373, 11.9395, 11.6477)*4, 0.035);

    APTracer::Entities::Scene_t scene;
    scene.add(water_mesh.triangles_, water_mesh.n_tris_);
    scene.add(sand_mesh.triangles_, sand_mesh.n_tris_);

    APTracer::Entities::ImgBufferOpenGL_t imgbuffer(1920, 1080);

    APTracer::Entities::TransformMatrix_t camera_transform;
    const std::array<double, 2> fov {9.0/16.0 * 80.0 * pi/180.0, 80.0 * pi/180.0};
    const std::array<unsigned int, 2> subpix {1, 1};
    std::list<Medium_t*> medium_list = {&air, &air};
    Vec3f min_sand_coord = sand_mesh.mincoord();
    APTracer::Cameras::Cam_t camera(&camera_transform, "images/output.png", Vec3f(0.0, 0.0, 1.0), fov, subpix, &imgbuffer, medium_list, &sky, 16, 1.0);
    camera.transformation_->translate(Vec3f(0.0, 2.0 * min_sand_coord[1], 0.0));
    camera.transformation_->rotateXAxis(-30.0 * pi/180);
    camera.transformation_->translate(Vec3f(0.0, 0.0, min_sand_coord[1]/2.0));
    camera.transformation_->rotateZAxis(pi + Rendering::n_timestep * Rendering::angle_step);
    camera.update();

    scene.build_acc();

    APTracer::Entities::OpenGLRenderer_t opengl_renderer(&scene, &camera, &imgbuffer);
    Rendering::renderer = &opengl_renderer;
    // We need to overload the render function with our own to apply movement etc
    opengl_renderer.render_function_ = openGL_accumulate;
    opengl_renderer.focus_point_ = Vec3f(0.0, 0.0, min_sand_coord[1]/2.0);
    Rendering::mesh_geometry = &water_mesh_geometry;
    Rendering::mesh = &water_mesh;
    Rendering::scene = &scene;
    Rendering::etas = etas;
    Rendering::n_points = n_grid_points;
    Rendering::omegas = omegas;

    opengl_renderer.initialise();
    opengl_renderer.render();

    return 0;
}

double get_max_depth(MeshGeometryUnstructured_t* mesh_geometry) {
    double depth = std::numeric_limits<double>::infinity();
    for (unsigned int i = 0; i < mesh_geometry->n_farfield_; ++i) {
        depth = std::min(depth, mesh_geometry->points_[mesh_geometry->farfield_[2 * i]][2]);
    }
    return depth;
}

void extrude_farfield(MeshGeometryUnstructured_t* mesh_geometry, double height, bool close) {
    unsigned int new_n_points = mesh_geometry->n_points_ + mesh_geometry->n_farfield_;
    unsigned int new_n_elements = mesh_geometry->n_elements_ + mesh_geometry->n_farfield_ * 2;
    unsigned int new_n_normals = mesh_geometry->n_normals_ + mesh_geometry->n_farfield_;

    // For closing stuff
    if (close) {
        new_n_points += 1;
        new_n_elements += mesh_geometry->n_farfield_;
        new_n_normals += 1;
    }

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
        Vec3f normal = mesh_geometry->points_[mesh_geometry->farfield_[2 * i]];
        normal[2] = 0;
        new_normals[i + mesh_geometry->n_normals_] = normal.normalize_inplace();
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
    if (close) {
        // Making new stuff
        new_points[mesh_geometry->n_points_ + mesh_geometry->n_farfield_] = Vec3f(0.0, 0.0, height);
        new_normals[mesh_geometry->n_points_ + mesh_geometry->n_farfield_] = Vec3f(0.0, 0.0, -1.0);

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
    unsigned int new_n_points = mesh_geometry->n_points_ + mesh_geometry->sum_n_wall_;
    unsigned int new_n_elements = mesh_geometry->n_elements_ + mesh_geometry->sum_n_wall_ * 2;
    unsigned int new_n_normals = mesh_geometry->n_normals_ + mesh_geometry->sum_n_wall_;

    // For closing stuff
    new_n_points += mesh_geometry->n_walls_;
    new_n_elements += mesh_geometry->sum_n_wall_;
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

    APTracer::Entities::Vec3f* centers = new Vec3f[mesh_geometry->n_walls_];

    unsigned int wall_index = 0;
    // Making new stuff
    for (unsigned int j = 0; j < mesh_geometry->n_walls_; ++j) {
       centers[j] = Vec3f(0.0);
        for (unsigned int i = 0; i < mesh_geometry->n_wall_[j]; ++i){
            centers[j] += mesh_geometry->points_[mesh_geometry->walls_[j][2 * i]];
        }
        centers[j] /= mesh_geometry->n_wall_[j];
        centers[j] += Vec3f(0.0, 0.0, height);

        for (unsigned int i = 0; i < mesh_geometry->n_wall_[j]; ++i){
            new_points[i + mesh_geometry->n_points_ + wall_index] = mesh_geometry->points_[mesh_geometry->walls_[j][2 * i]] + Vec3f(0.0, 0.0, height);
            Vec3f normal = new_points[i + mesh_geometry->n_points_ + wall_index] - centers[j];
            normal[2] = 0;
            new_normals[i + mesh_geometry->n_normals_ + wall_index] = normal.normalize_inplace(); 
        }
        wall_index += mesh_geometry->n_wall_[j];
    }

    wall_index = 0;
    // Adds two elements per boundary, created with the new points
    for (unsigned int j = 0; j < mesh_geometry->n_walls_; ++j) {
        for (unsigned int i = 0; i < mesh_geometry->n_wall_[j] - 1; ++i){
            new_elements[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * i] = mesh_geometry->walls_[j][2 * i];
            new_elements[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * i + 1] = mesh_geometry->n_points_ + i + wall_index;
            new_elements[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * i + 2] = mesh_geometry->n_points_ + i + 1 + wall_index;
            new_elements[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * i + 3] = mesh_geometry->walls_[j][2 * i + 1];
            new_elements[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * i + 4] = mesh_geometry->walls_[j][2 * i];
            new_elements[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * i + 5] = mesh_geometry->n_points_ + i + 1 + wall_index;
        }
        new_elements[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * (mesh_geometry->n_wall_[j] - 1)] = mesh_geometry->walls_[j][2 * (mesh_geometry->n_wall_[j] - 1)];
        new_elements[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * (mesh_geometry->n_wall_[j] - 1) + 1] = mesh_geometry->n_points_ + (mesh_geometry->n_wall_[j] - 1) + wall_index;
        new_elements[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * (mesh_geometry->n_wall_[j] - 1) + 2] = mesh_geometry->n_points_ + wall_index;
        new_elements[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * (mesh_geometry->n_wall_[j] - 1) + 3] = mesh_geometry->walls_[j][2 * (mesh_geometry->n_wall_[j] - 1) + 1];
        new_elements[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * (mesh_geometry->n_wall_[j] - 1) + 4] = mesh_geometry->walls_[j][2 * (mesh_geometry->n_wall_[j] - 1)];
        new_elements[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * (mesh_geometry->n_wall_[j] - 1) + 5] = mesh_geometry->n_points_ + wall_index;

        wall_index += mesh_geometry->n_wall_[j];
    }

    wall_index = 0;
    for (unsigned int j = 0; j < mesh_geometry->n_walls_; ++j) {
        for (unsigned int i = 0; i < mesh_geometry->n_wall_[j] - 1; ++i){
            new_element_normals[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * i] = mesh_geometry->n_normals_ + i + wall_index;
            new_element_normals[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * i + 1] = mesh_geometry->n_normals_ + i + wall_index;
            new_element_normals[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * i + 2] = mesh_geometry->n_normals_ + i + 1 + wall_index;
            new_element_normals[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * i + 3] = mesh_geometry->n_normals_ + i + 1 + wall_index;
            new_element_normals[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * i + 4] = mesh_geometry->n_normals_ + i + wall_index;
            new_element_normals[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * i + 5] = mesh_geometry->n_normals_ + i + 1 + wall_index;
        }
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * (mesh_geometry->n_wall_[j] - 1)] = mesh_geometry->n_normals_ + (mesh_geometry->n_wall_[j] - 1) + wall_index;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * (mesh_geometry->n_wall_[j] - 1) + 1] = mesh_geometry->n_normals_ + (mesh_geometry->n_wall_[j] - 1) + wall_index;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * (mesh_geometry->n_wall_[j] - 1) + 2] = mesh_geometry->n_normals_ + wall_index;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * (mesh_geometry->n_wall_[j] - 1) + 3] = mesh_geometry->n_normals_ + wall_index;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * (mesh_geometry->n_wall_[j] - 1) + 4] = mesh_geometry->n_normals_ + (mesh_geometry->n_wall_[j] - 1) + wall_index;
        new_element_normals[3 * mesh_geometry->n_elements_ + 6 * wall_index + 6 * (mesh_geometry->n_wall_[j] - 1) + 5] = mesh_geometry->n_normals_ + wall_index;

        wall_index += mesh_geometry->n_wall_[j];
    }

    // Closing
    // Making new stuff
    for (unsigned int j = 0; j < mesh_geometry->n_walls_; ++j) {
        new_points[mesh_geometry->n_points_ + mesh_geometry->sum_n_wall_ + j] = centers[j];
    }
    new_normals[mesh_geometry->n_points_ + mesh_geometry->sum_n_wall_] = Vec3f(0.0, 0.0, 1.0);

    // Adds one element per boundary, created with the new points
    wall_index = 0;
    for (unsigned int j = 0; j < mesh_geometry->n_walls_; ++j) {
        for (unsigned int i = 0; i < mesh_geometry->n_wall_[j] - 1; ++i){
            new_elements[3 * mesh_geometry->n_elements_ + 3 * wall_index + 6 * mesh_geometry->sum_n_wall_ + 3 * i] = mesh_geometry->n_points_ + i + wall_index;
            new_elements[3 * mesh_geometry->n_elements_ + 3 * wall_index + 6 * mesh_geometry->sum_n_wall_ + 3 * i + 1] = mesh_geometry->n_points_ + i + 1 + wall_index;
            new_elements[3 * mesh_geometry->n_elements_ + 3 * wall_index + 6 * mesh_geometry->sum_n_wall_ + 3 * i + 2] = mesh_geometry->n_points_ + mesh_geometry->sum_n_wall_ + j;
        }
        new_elements[3 * mesh_geometry->n_elements_ + 3 * wall_index + 6 * mesh_geometry->sum_n_wall_ + 3 * (mesh_geometry->n_wall_[j] - 1)] = mesh_geometry->n_points_ + (mesh_geometry->n_wall_[j] - 1) + wall_index;
        new_elements[3 * mesh_geometry->n_elements_ + 3 * wall_index + 6 * mesh_geometry->sum_n_wall_ + 3 * (mesh_geometry->n_wall_[j] - 1) + 1] = mesh_geometry->n_points_ + wall_index;
        new_elements[3 * mesh_geometry->n_elements_ + 3 * wall_index + 6 * mesh_geometry->sum_n_wall_ + 3 * (mesh_geometry->n_wall_[j] - 1) + 2] = mesh_geometry->n_points_ + mesh_geometry->sum_n_wall_ + j;
        wall_index += mesh_geometry->n_wall_[j];
    }

    wall_index = 0;
    for (unsigned int j = 0; j < mesh_geometry->n_walls_; ++j) {
        for (unsigned int i = 0; i < mesh_geometry->n_wall_[j]; ++i){
            new_element_normals[3 * mesh_geometry->n_elements_ + 3 * wall_index + 6 * mesh_geometry->sum_n_wall_ + 3 * i] = mesh_geometry->n_points_ + mesh_geometry->sum_n_wall_;
            new_element_normals[3 * mesh_geometry->n_elements_ + 3 * wall_index + 6 * mesh_geometry->sum_n_wall_ + 3 * i + 1] = mesh_geometry->n_points_ + mesh_geometry->sum_n_wall_;
            new_element_normals[3 * mesh_geometry->n_elements_ + 3 * wall_index + 6 * mesh_geometry->sum_n_wall_ + 3 * i + 2] = mesh_geometry->n_points_ + mesh_geometry->sum_n_wall_;
        }
        wall_index += mesh_geometry->n_wall_[j];
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
    delete [] centers;
}

std::vector<std::complex<double>> get_eta(std::string filename, double &amplitude, double &omega) {
    std::string line;
    std::string token;
    size_t n_points;

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

void timestep(MeshGeometryUnstructured_t* mesh_geometry, MeshUnstructured_t* mesh, APTracer::Entities::Scene_t* scene, std::vector<std::vector<std::complex<double>>> etas, unsigned int n_points, double time, std::vector<double> omegas) {
    for (unsigned int i = 0; i < n_points; ++i) {
        mesh_geometry->points_[i][2] = 0.0;
        for (unsigned int j = 0; j < etas.size(); ++j) {
            mesh_geometry->points_[i][2] += std::real(etas[j][i] * std::exp(std::complex<double>(0.0, -1.0) * omegas[j] * time));
        }
    }

    mesh_geometry->computeNormals(n_points);
    mesh->update();
    scene->build_acc(); // Slower than removing and adding the shapes, but will maybe fix issues with triangles dissapearing.
}

void openGL_accumulate() {
    auto t_start = std::chrono::high_resolution_clock::now();
    Rendering::renderer->accumulate();
    auto t_end = std::chrono::high_resolution_clock::now();

    std::cout << "Iteration " << Rendering::renderer->imgbuffer_->updates_ << " done in " 
        << std::chrono::duration<double, std::milli>(t_end - t_start).count()/1000.0 
        << "s." << std::endl;
    
    if (Rendering::renderer->imgbuffer_->updates_ == Rendering::write_interval) {
        ++Rendering::n_timestep;
        std::cout << "Timestep " << Rendering::n_timestep << ", t = " << Rendering::time << std::endl;
        std::ostringstream oss;
        oss << "images/image_"<< std::setfill('0') << std::setw(4) << Rendering::n_timestep << ".png";
        auto t_write_start = std::chrono::high_resolution_clock::now();
        Rendering::renderer->camera_->write(oss.str());
        auto t_write_end = std::chrono::high_resolution_clock::now();
        std::cout << "\tWriting done in " 
            << std::chrono::duration<double, std::milli>(t_write_end - t_write_start).count()/1000.0 
            << "s." << std::endl;

        auto t_step_start = std::chrono::high_resolution_clock::now();
        Rendering::time += Rendering::delta_time;
        timestep(Rendering::mesh_geometry, Rendering::mesh, Rendering::scene, Rendering::etas, Rendering::n_points, Rendering::time, Rendering::omegas);
        Rendering::renderer->camera_->transformation_->rotateZAxis(Rendering::angle_step);
        Rendering::renderer->camera_->update();
        Rendering::renderer->resetDisplay();
        auto t_step_end = std::chrono::high_resolution_clock::now();
        std::cout << "\tUpdating done in " 
            << std::chrono::duration<double, std::milli>(t_step_end - t_step_start).count()/1000.0 
            << "s." << std::endl;
    }  
}
