
#include <fstream>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"
#include <stdio.h>

#include "Mesh.hpp"


typedef Mesh<double, double, double > MeshType;

struct NodePosition {

  template<typename NODE>
  Point operator()(const NODE & n) {
    return Point(n.position().x, n.position().y, n.position().z);
  }
};

void populate_mesh(MeshType& mesh, std::vector<typename MeshType::node_type>& mesh_node, std::ifstream& nodes_file, std::ifstream& tris_file) {

  Point p;

  while(CS207::getline_parsed(nodes_file, p)) {
    mesh_node.push_back(mesh.add_node(p));
  }
  
  std::array<int,3> t;

  while(CS207::getline_parsed(tris_file, t)) {
    mesh.add_triangle(mesh_node[t[0]], mesh_node[t[1]], mesh_node[t[2]]);
  }
  std::cout << "hi";

}

int main(int argc, char* argv[])
{
 
  std::ifstream nodes_file;
  std::ifstream tris_file;
 
  int sim_file_start;

  if(std::string(argv[1]).find(".nodes") == std::string::npos) { //not in right format

    sim_file_start = 2;

    std::ifstream mesh_file(argv[1]);
    std::string s;
    for (int i = 0; i < 9; ++i) {
     getline(mesh_file, s);
    }
    
    std::ofstream nodes_output;
    nodes_output.open("generated.nodes");

    getline(mesh_file, s);

    double x, y;
    int num;
    while(s.find("ENDOFSECTION") == std::string::npos) {
      sscanf (s.c_str(), "%d %lf %lf", &num, &x, &y);
      nodes_output << x << "\t" << y << "\t" << "0" << std::endl;
      getline(mesh_file, s);
    }

    nodes_output.close();

    std::ofstream tris_output;
    tris_output.open("generated.tris");

    getline(mesh_file,s); getline(mesh_file,s);

    int con1, con2, con3;
    int node1, node2, node3;
    while(s.find("ENDOFSECTION") == std::string::npos) {
      sscanf (s.c_str(), "%d %d %d %d %d %d", &con1, &con2, &con3, &node1, &node2, &node3);
      tris_output << --node1 << "\t" << --node2 << "\t" << --node3  << std::endl;
      getline(mesh_file, s);
    }
    
    tris_output.close();

    nodes_file.open("generated.nodes");
    tris_file.open("generated.tris");
  }

  else {
   
    sim_file_start = 3;
    nodes_file.open(argv[1]);
    tris_file.open(argv[2]);
  } 

  MeshType mesh;
  std::vector<typename MeshType::node_type> mesh_node;
  
  populate_mesh(mesh, mesh_node, nodes_file, tris_file);


  CS207::SDLViewer viewer;
  viewer.launch();
  auto node_map = viewer.empty_node_map(mesh);
  viewer.add_nodes(mesh.node_begin(), mesh.node_end(), CS207::DefaultColor(), NodePosition(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
  viewer.center_view();

  for (int i = sim_file_start; i < argc; i++) {}

}
