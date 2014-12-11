
#include <fstream>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"


#include "Mesh.hpp"


typedef Mesh<double, double, double > MeshType;

struct NodePosition {

  template<typename NODE>
  Point operator()(const NODE & n) {
    return Point(n.position().x, n.position().y, n.position().z);
  }
};



int main(int argc, char* argv[])
{
  MeshType mesh; 
  std::vector<typename MeshType::node_type> mesh_node;
  std::ifstream nodes_file(argv[1]);
  Point p;

  while(CS207::getline_parsed(nodes_file, p)) {
    mesh_node.push_back(mesh.add_node(p));
  }

  std::ifstream tris_file(argv[2]);
  std::array<int,3> t;
  while(CS207::getline_parsed(tris_file, t)) {
    mesh.add_triangle(mesh_node[t[0]], mesh_node[t[1]], mesh_node[t[2]]);
  }


  CS207::SDLViewer viewer;
  viewer.launch();
  auto node_map = viewer.empty_node_map(mesh);
  viewer.add_nodes(mesh.node_begin(), mesh.node_end(), CS207::DefaultColor(), NodePosition(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
  viewer.center_view();

  //for each file, reload mesh.

}
