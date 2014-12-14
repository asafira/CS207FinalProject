
#include <fstream>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"
#include <stdio.h>

#include "Mesh.hpp"



struct NodePosition {

  template<typename NODE>
  Point operator()(const NODE & n) {
    return Point(n.position().x, n.position().y, n.position().z);
  }
};

struct node_info {

  int num_contribs;
  double sum;
};

typedef Mesh<node_info, double, double> MeshType;

void populate_mesh(MeshType& mesh, std::vector<typename MeshType::node_type>& mesh_node, std::ifstream& nodes_file, std::ifstream& tris_file) {

  Point p;

  while(CS207::getline_parsed(nodes_file, p)) {
    mesh_node.push_back(mesh.add_node(p));
  }
  
  std::array<int,3> t;

  while(CS207::getline_parsed(tris_file, t)) {
    mesh.add_triangle(mesh_node[t[0]], mesh_node[t[1]], mesh_node[t[2]]);
  }
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


  for (auto node_it = mesh.node_begin(); node_it != mesh.node_end(); ++node_it)
    (*node_it).value() = {0,0.0};

  CS207::SDLViewer viewer;
  viewer.launch();
  auto node_map = viewer.empty_node_map(mesh);
  viewer.add_nodes(mesh.node_begin(), mesh.node_end(), CS207::DefaultColor(), NodePosition(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
  viewer.center_view();

  double ux, uy, p, curl;
  std::string curr_string;
  std::size_t line_number;
  double max = 0;
  double min = 0;
  for (int i = sim_file_start; i < argc; i++) {
    std::ifstream sim(argv[i]);
    int count = 0;
    for (auto tri_it = mesh.triangle_begin(); tri_it != mesh.triangle_end(); ++tri_it) {
      for (int node_num = 0; node_num < 3; node_num++) {
        if (sim.eof()) {
          std::cout << "got to end of file before iterating over all triangles" << std::endl;
          break;
        }
        getline(sim, curr_string); 
        sscanf(curr_string.c_str(), "%lf %lf %lf %lf", &ux, &uy, &p, &curl);
        if (max < curl) {max = curl;}
        if (min > curl) {min = curl;}

        node_info current = (*tri_it).node(node_num).value();
        (*tri_it).node(node_num).value() = {1+current.num_contribs, curl+current.sum};
      }
/* 
      count++;
      if (count == 10) {std::cout << ((*tri_it).node(0).value().sum/(*tri_it).node(0).value().num_contribs) << std::endl; }

*/
    }
/*
    if (i % 2 == 0 || i % 2 == 1) {
*/
      viewer.add_nodes(mesh.node_begin(), mesh.node_end(), CS207::SimulationColor(min, max), node_map);
      viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
      viewer.set_label(i);
      viewer.center_view();
  /*  }

    else {

      viewer.add_nodes(mesh.node_begin(), mesh.node_end(), CS207::DefaultColor(), node_map);
      viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
      viewer.set_label(i);
      viewer.center_view();
  
    }
*/
    for (auto node_it = mesh.node_begin(); node_it != mesh.node_end(); ++node_it)
      (*node_it).value() = {0,0.0};

    max = 0; min = 0;

    std::cout << "-----------------------------------------" << std::endl;
    //CS207::sleep(.1);
  }

}
