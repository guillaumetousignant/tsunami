#include "another_path_tracer.h"
#include "entities/MeshGeometryUnstructured_t.h"

#include <string>

int main(int argc, char **argv){
    if (argc < 3) {
        std::cerr << "Error: input arguments with format 'mesh_file data_file'. Exiting." << std::endl;
    }
    std::string mesh_file = argv[1];
    std::string data_file = argv[2];

    MeshGeometryUnstructured_t mesh(mesh_file);

    return 0;
}