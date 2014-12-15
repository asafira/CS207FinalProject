
#include <fstream>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"
#include <stdio.h>

#include "Mesh.hpp"


// Node type to keep track of the number of times a node has a solution, and the sum
// of each solution
struct node_type {

  int num_contribs;
  double sum;
};



typedef Mesh<node_type, double, double> MeshType;


/* @brief: Create a mesh and vector of nodes from properly 
           formatted node and triangle files.

   @param[in] @a mesh is a mesh to place the triangles of 
              @a tris_file 
   @param[in] mesh_node is a vector of nodes to place the
              nodes of @a nodes_file
   @param[in] @nodes_file is a nodes input file with each
              line containing three coordinates for each node
   @param[in] tris_file is a triangle input file with each
              line containing three coordinates for each 
              triangle

   @post mesh is populated with the triangles in tris_file,
         and mesh_node is populated with the nodes in
         nodes_file
*/

void populate_mesh(MeshType& mesh, 
                   std::vector<typename MeshType::node_type>& mesh_node, 
                   std::ifstream& nodes_file, std::ifstream& tris_file) 
{

  Point p;

  // Add nodes to mesh and mesh_node vector
  while(CS207::getline_parsed(nodes_file, p)) {
    mesh_node.push_back(mesh.add_node(p));
  }
  
  // Add triangles to mesh
  std::array<int,3> t;
  while(CS207::getline_parsed(tris_file, t)) {
    mesh.add_triangle(mesh_node[t[0]], mesh_node[t[1]], mesh_node[t[2]]);
  }
}


/* @brief resets the value of each node to {0, 0.0} in the mesh
   
   @param[in] mesh contains nodes of type node_type
   @post All nodes in @a mesh are set to {0, 0.0}
*/

void reset_node_values (MeshType& mesh) {
    
  for (auto node_it = mesh.node_begin(); node_it != mesh.node_end(); ++node_it) {
    (*node_it).value() = {0,0.0};
  }
  
}

/* @brief Takes in a mesh and solution files to visualize
           the results in the solution files

   @param[in] @a argc contains the number of command-line parameters
   @param[in] @a argv[] contains the commdn-line parameters; either
              1. two mesh files, file.nodes and file.tets, are the first
                 two arguments. The former contains the node coordiinates
                 on each line, the latter contains three node indices on
                 each line for each triangle in the mesh.
              2.  one mesh file used in the nudg++ simulation software,
                 (.neu extension)
        
                  Each of the above is followed by an arbitrary number of
                  solutions files, containing 3*num_triangles lines, each 
                  for nodes corresponding to each triangle.
   @post generated.nodes and generated.tris files will result if option
         2 from the above input options is used. These files are generated
         files of the format descried in option 1 above.
   @post A visualization will pop up for the curl of the velocity field
         of the fluid on the given mesh.
*/

int main(int argc, char* argv[])
{
 
  // Initialize nodes and triangles files
  std::ifstream nodes_file;
  std::ifstream tris_file;
 
 
  int sim_file_start; //keeps track of what argument begins simulation files

  // Check if input is in the correct format
  if(std::string(argv[1]).find(".nodes") == std::string::npos) 
    {

    sim_file_start = 2; //update when simulation files start

    // Load in the first argument
    std::ifstream mesh_file(argv[1]);

    // Skip the first 10 lines
    std::string s;
    for (int i = 0; i < 9; ++i) {
     getline(mesh_file, s);
    }
    
    // Open new file to output node data in correct format
    std::ofstream nodes_output;
    nodes_output.open("generated.nodes");

    // Skip first line
    getline(mesh_file, s);

    double x, y;
    int num;
  
    // Read in lines until "ENDOFSECTION" is reached
    while(s.find("ENDOFSECTION") == std::string::npos) {
      sscanf (s.c_str(), "%d %lf %lf", &num, &x, &y);
      nodes_output << x << "\t" << y << "\t" << "0" << std::endl;
      getline(mesh_file, s);
    }

    nodes_output.close();

    // Open new file to output triangle data in correct format
    std::ofstream tris_output;
    tris_output.open("generated.tris");

    // Skip 2 lines
    getline(mesh_file,s); getline(mesh_file,s);

    int con1, con2, con3;
    int node1, node2, node3;
    while(s.find("ENDOFSECTION") == std::string::npos) {
      sscanf (s.c_str(), "%d %d %d %d %d %d", &con1, &con2, &con3, &node1, &node2, &node3);
      tris_output << --node1 << "\t" << --node2 << "\t" << --node3  << std::endl;
      getline(mesh_file, s);
    }
    
    tris_output.close();

    // Open the files for later visualization
    nodes_file.open("generated.nodes");
    tris_file.open("generated.tris");
  }

  // If files were given in the correct format, just open them
  else {
   
    sim_file_start = 3;
    nodes_file.open(argv[1]);
    tris_file.open(argv[2]);
  } 

  
  MeshType mesh;
  std::vector<typename MeshType::node_type> mesh_node;
  
  populate_mesh(mesh, mesh_node, nodes_file, tris_file);

  nodes_file.close(); tris_file.close();
  
  // Set nodes to have 0 values
  reset_node_values(mesh);

  // Initialize the viewer
  CS207::SDLViewer viewer;
  viewer.launch();
  auto node_map = viewer.empty_node_map(mesh);

  
  double ux, uy, p, curl;
  std::string curr_string;
  
  double max = 0;
  double min = 0;
  
  // Go through solutions files
  for (int i = sim_file_start; i < argc; i++) {
    
    // Open next file
    std::ifstream sim(argv[i]);
    
    // Skip first line
    getline(sim, curr_string);

    // The solution for the 
    for (auto tri_it = mesh.triangle_begin(); tri_it != mesh.triangle_end(); ++tri_it) {
      
      // Read next three lines for node data for triangle (*tri_it)
      for (int node_num = 0; node_num < 3; node_num++) {
        
        // Check if end of file error (should not reach)
        if (sim.eof()) {
          std::cout << "got to end of file before iterating over all triangles" << std::endl;
          break;
        }

        // Get the values from the solution file line
        getline(sim, curr_string); 
        sscanf(curr_string.c_str(), "%lf %lf %lf %lf", &ux, &uy, &p, &curl);

        // Update max and min values
        if (max < curl) {max = curl;}
        if (min > curl) {min = curl;}

        // Update node in mesh
        node_type current = (*tri_it).node(node_num).value();
        (*tri_it).node(node_num).value() = {++current.num_contribs, curl+current.sum};
      }
    }
      
    // View solution
    viewer.add_nodes(mesh.node_begin(), mesh.node_end(), 
			CS207::SimulationColor(min, max), node_map);
    viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);
    viewer.set_label(i);
    viewer.center_view();
    
    reset_node_values(mesh);
    max = 0; min = 0;
    
    //CS207::sleep(.5);  // User-settable sleep time
  }

}
