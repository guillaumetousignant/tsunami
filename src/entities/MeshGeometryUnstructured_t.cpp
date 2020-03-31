#include "entities/MeshGeometryUnstructured_t.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>

using APTracer::Entities::Vec3f;

MeshGeometryUnstructured_t::MeshGeometryUnstructured_t(const std::string &filename){
    std::string ext = filename.substr(filename.find_last_of(".") + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    if(ext == "su2") {
        readSU2(filename);
    }
    else {
        std::cerr << "Error: file extension '" << ext << "' not recognized. Exiting." << std::endl;
        return;
    }
}

MeshGeometryUnstructured_t::~MeshGeometryUnstructured_t(){
    if (points_ != nullptr){
        delete [] points_;
    }

    if (elements_ != nullptr){
        delete [] elements_;
    }

    if (normals_ != nullptr){
        delete [] normals_;
    }

    if (wall_ != nullptr) {
        delete [] wall_;
    }

    if (farfield_ != nullptr) {
        delete [] farfield_;
    }
}

void MeshGeometryUnstructured_t::readSU2(const std::string &filename){
    std::string line;
    std::string token;

    std::ifstream meshfile(filename);
    if (!meshfile.is_open()) {
        std::cerr << "Error: file '" << filename << "' could not be opened. Exiting." << std::endl;
        return;
    }

    std::getline(meshfile, line);
    std::getline(meshfile, line);
    std::getline(meshfile, line);
    std::istringstream liness(line);
    liness >> token;
    if (token == "NPOIN="){
        liness >> n_points_;
    }
    else {
        std::cerr << "Error: expected marker 'NPOIN=', found '" << token << "'. Exiting." << std::endl;
        return;
    }

    points_ = new Vec3f[n_points_];
    for (unsigned int i = 0; i < n_points_; ++i) {
        std::getline(meshfile, line);
        std::istringstream liness2(line);
        liness2 >> points_[i][0] >> points_[i][1] >> points_[i][2];
    }

    std::getline(meshfile, line);
    std::getline(meshfile, line);
    std::istringstream liness3(line);
    liness3 >> token;
    if (token == "NELEM="){
        liness3 >> n_elements_;
    }
    else {
        std::cerr << "Error: expected marker 'NELEM=', found '" << token << "'. Exiting." << std::endl;
        return;
    }

    elements_ = new unsigned int[3 * n_elements_];
    for (unsigned int i = 0; i < n_elements_; ++i) {
        std::getline(meshfile, line);
        std::istringstream liness2(line);
        liness2 >> token;
        if (token != "5") {
            std::cerr << "Error: expected token '5', found '" << token << "'. Exiting." << std::endl;
            return;
        }

        liness2 >> elements_[3 * i] >> elements_[3 * i + 1] >> elements_[3 * i + 2];
    }

    unsigned int n_markers;
    std::getline(meshfile, line);
    std::istringstream liness4(line);
    liness4 >> token;
    if (token == "NMARK="){
        liness4 >> n_markers;
    }
    else {
        std::cerr << "Error: expected marker 'NMARK=', found '" << token << "'. Exiting." << std::endl;
        return;
    }

    for (unsigned int i = 0; i < n_markers; ++i) {
        std::string type;
        std::getline(meshfile, line);
        std::istringstream liness5(line);
        liness5 >> token;
        if (token == "MARKER_TAG="){
            liness5 >> type;
        }
        else {
            std::cerr << "Error: expected marker 'MARKER_TAG=', found '" << token << "'. Exiting." << std::endl;
            return;
        }

        if (type == "wall") {
            std::getline(meshfile, line);
            std::istringstream liness2(line);
            liness2 >> token;
            if (token == "MARKER_ELEMS="){
                liness2 >> n_wall_;
            }
            else {
                std::cerr << "Error: expected marker 'MARKER_ELEMS=', found '" << token << "'. Exiting." << std::endl;
                return;
            }

            wall_ = new unsigned int[2 * n_wall_];
            for (unsigned int j = 0; j < n_wall_; ++j) {

                std::getline(meshfile, line);
                std::istringstream liness6(line);

                liness6 >> token;
                if (token != "3") {
                    std::cerr << "Error: expected token '3', found '" << token << "'. Exiting." << std::endl;
                    return;
                }

                liness6 >> wall_[2 * j] >> wall_[2 * j + 1];
            }
        }
        else if (type == "farfield") {
            std::getline(meshfile, line);
            std::istringstream liness2(line);
            liness2 >> token;
            if (token == "MARKER_ELEMS="){
                liness2 >> n_farfield_;
            }
            else {
                std::cerr << "Error: expected marker 'MARKER_ELEMS=', found '" << token << "'. Exiting." << std::endl;
                return;
            }

            farfield_ = new unsigned int[2 * n_farfield_];
            for (unsigned int j = 0; j < n_farfield_; ++j) {

                std::getline(meshfile, line);
                std::istringstream liness6(line);

                liness6 >> token;
                if (token != "3") {
                    std::cerr << "Error: expected token '3', found '" << token << "'. Exiting." << std::endl;
                    return;
                }

                liness6 >> farfield_[2 * j] >> farfield_[2 * j+ 1];
            }
        }
        else {
            std::cerr << "Error: expected marker tag 'wall' or 'farfield', found '" << type << "'. Exiting." << std::endl;
            return;
        }
    }
    meshfile.close();

    normals_ = new Vec3f[n_points_];

    computeNodeToFace();
    computeNormals();
}

void MeshGeometryUnstructured_t::computeNodeToFace() {
    std::vector<std::vector<unsigned int>> point_to_elements_(n_points_);
    for (unsigned int i = 0; i < n_points_; ++i){
        for (unsigned int j = 0; j < n_elements_; ++j) {
            for (unsigned int k = 0; k < 3; ++k){
                if (elements_[3 * j + k] == i) {
                    point_to_elements_[i].push_back(j);
                    break;
                }
            }
        }
    }
}

void MeshGeometryUnstructured_t::computeNormals() {
    for (unsigned int i = 0; i < n_points_; ++i) {
        Vec3f normal;
        for (unsigned int j = 0; j < point_to_elements_[i].size(); ++j) {
            normal += (points_[elements_[3 * i + 1]] - points_[elements_[3 * i]]).cross(points_[elements_[3 * i + 2]] - points_[elements_[3 * i]]).normalize_inplace();
        }
        normals_[i] = normal.normalize();
    }
}