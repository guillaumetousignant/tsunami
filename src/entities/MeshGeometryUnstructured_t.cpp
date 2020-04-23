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

    if (element_normals_ != nullptr) {
        delete [] element_normals_;
    }

    if (normals_ != nullptr){
        delete [] normals_;
    }

    if (walls_ != nullptr) {
        for (unsigned int i = 0; i < n_walls_; ++i) {
            if (walls_[i] != nullptr) {
                delete [] walls_[i];
            } 
        }
        delete [] walls_;
    }

    if (farfield_ != nullptr) {
        delete [] farfield_;
    }

    if (n_wall_ != nullptr) {
        delete [] n_wall_;
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

    n_normals_ = n_points_;
    elements_ = new unsigned int[3 * n_elements_];
    element_normals_ = new unsigned int[3 * n_elements_];
    for (unsigned int i = 0; i < n_elements_; ++i) {
        std::getline(meshfile, line);
        std::istringstream liness2(line);
        liness2 >> token;
        if (token != "5") {
            std::cerr << "Error: expected token '5', found '" << token << "'. Exiting." << std::endl;
            return;
        }
        unsigned int val0, val1, val2;
        liness2 >> val0 >> val1 >> val2;
        elements_[3 * i] = val0 - 1;
        elements_[3 * i + 1] = val1 - 1;
        elements_[3 * i + 2] = val2 - 1;
        element_normals_[3 * i] = elements_[3 * i];
        element_normals_[3 * i + 1] = elements_[3 * i + 1];
        element_normals_[3 * i + 2] = elements_[3 * i + 2];
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

    n_walls_ = n_markers - 1;
    n_wall_ = new unsigned int[n_walls_];
    unsigned int wall_index = 0;
    walls_ = new unsigned int*[n_walls_];

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
                liness2 >> n_wall_[wall_index];
            }
            else {
                std::cerr << "Error: expected marker 'MARKER_ELEMS=', found '" << token << "'. Exiting." << std::endl;
                return;
            }

            walls_[wall_index] = new unsigned int[2 * n_wall_[wall_index]];
            for (unsigned int j = 0; j < n_wall_[wall_index]; ++j) {

                std::getline(meshfile, line);
                std::istringstream liness6(line);

                liness6 >> token;
                if (token != "3") {
                    std::cerr << "Error: expected token '3', found '" << token << "'. Exiting." << std::endl;
                    return;
                }

                unsigned int val0, val1;
                liness6 >> val0 >> val1;
                walls_[wall_index][2 * j] = val0 - 1;
                walls_[wall_index][2 * j + 1] = val1 - 1;
            }
            ++wall_index;
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

                unsigned int val0, val1;
                liness6 >> val0 >> val1;
                farfield_[2 * j] = val0 - 1;
                farfield_[2 * j + 1] = val1 - 1;
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
    computeNormals(n_points_);
}

void MeshGeometryUnstructured_t::computeNodeToFace() {
    point_to_elements_ = std::vector<std::vector<unsigned int>>(n_points_, std::vector<unsigned int>());
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

void MeshGeometryUnstructured_t::computeNormals(unsigned int n_points) {
    for (unsigned int i = 0; i < n_points; ++i) {
        Vec3f normal;
        std::vector<unsigned int> point_elements = point_to_elements_[i];
        for (unsigned int j = 0; j < point_elements.size(); ++j) {
            unsigned int element_index = point_elements[j];
            unsigned int point_index0 = elements_[3 * element_index];
            unsigned int point_index1 = elements_[3 * element_index + 1];
            unsigned int point_index2 = elements_[3 * element_index + 2];
            Vec3f point0 = points_[point_index0];
            Vec3f point1 = points_[point_index1];
            Vec3f point2 = points_[point_index2];

            normal += (point2 - point0).cross(point1 - point0).normalize_inplace();
            
            //normal += (points_[elements_[3 * point_to_elements_[i][j] + 1]] - points_[elements_[3 * point_to_elements_[i][j]]]).cross(points_[elements_[3 * point_to_elements_[i][j] + 2]] - points_[elements_[3 * point_to_elements_[i][j]]]).normalize_inplace(); // wow
        }
        normals_[i] = normal.normalize();
    }
}

void MeshGeometryUnstructured_t::verify() {
    std::cout << "Tere are " << n_points_ << " points:" << std::endl;
    for (unsigned int i = 0; i < n_points_; ++i) {
        std::cout << "    " << i << "    " << points_[i] << std::endl;
    }

    std::cout << std::endl << "Tere are " << n_elements_ << " elements, with points:" << std::endl;
    for (unsigned int i = 0; i < n_elements_; ++i) {
        std::cout << "    " << i << "    " << elements_[3 * i] << "    " << elements_[3 * i + 1] << "    " << elements_[3 * i + 2] << std::endl;
    }

    std::cout << std::endl << "Tere are " << n_walls_ << " walls." << std::endl;
    for (unsigned int j = 0; j < n_walls_; ++j) {
        std::cout << std::endl << "    Tere are " << n_wall_[j] << " wall elements in wall " << j << ", with points:" << std::endl;
        for (unsigned int i = 0; i < n_wall_[j]; ++i) {
            std::cout << "        " << i << "    " << walls_[j][2 * i] << "    " << walls_[j][2 * i + 1] << std::endl;
        }
    }

    std::cout << std::endl << "Tere are " << n_farfield_ << " farfield elements, with points:" << std::endl;
    for (unsigned int i = 0; i < n_farfield_; ++i) {
        std::cout << "    " << i << "    " << farfield_[2 * i] << "    " << farfield_[2 * i + 1] << std::endl;
    }

    std::cout << "Tere are " << n_normals_ << " normals:" << std::endl;
    for (unsigned int i = 0; i < n_normals_; ++i) {
        std::cout << "    " << i << "    " << normals_[i] << std::endl;
    }

    std::cout << std::endl << "Tere are " << n_elements_ << " elements, with normals:" << std::endl;
    for (unsigned int i = 0; i < n_elements_; ++i) {
        std::cout << "    " << i << "    " << element_normals_[3 * i] << "    " << element_normals_[3 * i + 1] << "    " << element_normals_[3 * i + 2] << std::endl;
    }
}