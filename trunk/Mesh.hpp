#pragma once
/** @file Mesh.hpp
 * @brief A Mesh is composed of nodes, edges, and triangles such that:
 *  -- All triangles have three nodes and three edges.
 *  -- All edges belong to at least one triangle and at most two triangles.
 */

/** @class Mesh
 * @brief A template for 3D triangular meshes.
 *
 * Users can add triangles and retrieve nodes, edges, and triangles.
 */
template <typename N, typename E, typename T>
class Mesh {
 
  public:
  //class definitions
  class Node;
  class Edge;
  class Triangle;
  class NodeIterator;
  class EdgeIterator;
  class TriangleIterator;
  class incident_n2t_iterator;

  // define node_type for user
  typedef Node node_type;

  /** Type of indexes and sizes. Return type of Mesh::num_nodes(). */
  typedef unsigned size_type;

  //type of node_value
  typedef N node_value_type;

  //type of edge_value
  typedef E edge_value_type;

  //type of node_value
  typedef T triangle_value_type;

  //triplet of signed ints
  typedef std::array<signed, 3> int_triplet;
  
  //triplet of unsigned ints
  typedef std::array<unsigned, 3> unsigned_triplet;

  //doublet of unsigned ints
  typedef std::array<unsigned, 2> unsigned_doublet;

  /////////////////
  //  Mesh NODES //
  /////////////////
  class Node : private totally_ordered<Node> {
  public: 
    /* construct an invalid Node O(1) */
    Node() : node_id_(UINT_MAX){
    }

    /** Return this node's position. O(1) */
    const Point& position() const {
      return mesh_->nodes_[node_id_].point;
    }
    
    /** Public function that allows to modify point location O(1) */
    Point& position(){
      return mesh_->nodes_[node_id_].point;
    }

   /** Return this node's index, a number in the range [0, mesh_size). O(1) */
    size_type index() const{
      return node_id_;
    }

    /** Test whether this node and @a x are equal. O(1) */
    bool operator==(const Node& x) const {
      return (mesh_==x.mesh_ && node_id_==x.node_id_);
    }

    /** Test whether this node is less than @a x in the global order. O(1) */
    bool operator<(const Node& x) const {
      if (mesh_ < x.mesh_)    //allow comparison between different meshs
        return true;
      return (mesh_== x.mesh_ && node_id_ < x.node_id_);
    }

    /** public valid function to check if a node is valid O(1) */
    bool valid() const {
      return node_id_<mesh_->nodes_.size();
    }

    /* return the reference to user defined value of this node O(1) */
    node_value_type& value(){
      return mesh_->nodes_[node_id_].node_value;
    }

    /*return the reference to user defined value of this node O(1) */
    const node_value_type& value() const{
      return mesh_->nodes_[node_id_].node_value;
    }

    /** return degree (number of connecting triangles) of this node O(1) */
    size_type num_triangles() const{
      return mesh_->adj_n2t_[node_id_].size();
    }

    //returns the first incident_n2t_iterator of the current node O(1)
    incident_n2t_iterator incident_n2t_begin() const{
      return incident_n2t_iterator(mesh_,node_id_,0);
    }

    //returns an invalid incident_n2t_iterator O(1)
    incident_n2t_iterator incident_n2t_end() const{
      return incident_n2t_iterator(mesh_,node_id_,num_triangles());
    }  

  private:
    //Allow Mesh to access Node's private member data and functions
    friend class Mesh;
    
    //private constructor
    Node(const Mesh* mesh, size_type node_id)
      :mesh_(const_cast<Mesh*>(mesh)), node_id_(node_id){
    }

    //private variables:
    Mesh* mesh_;
    size_type node_id_;
  };
  //end NODE

  /** Return the number of nodes in the mesh. */
  size_type num_nodes() const {
    return nodes_.size();
  }

  /** add_node
   * @pre assume the node is not being previous added  
   O(1) */
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()){
    nodes_.push_back(node_data(position,node_value));
    adj_n2t_.push_back(std::vector<size_type>());
    return Node(this,nodes_.size()-1);
  }

  class Edge: private totally_ordered<Edge> {
  public: 
    /** Construct an invalid Edge. O(1) */
    Edge() : edge_id_(UINT_MAX) {
    }

    /** return length of the two connecting nodes O(1) */
    double length() const{
      //return norm_2(node1().position()-node2().position());
      return mesh_->edges_[edge_id_].edge_length;
    }
    
    /** return normal of the current edge */
    Point normal() const{
      return mesh_->edges_[edge_id_].normal;
    }

    /**return vector of the two connecting nodes O(1) */ 
    Point vector_length() const{
      return (node1().position() - node2().position());
    }

    /**return true if this edge has node n */
    bool has_node(const Node& n) const{
      return n.index()== node1().index() || n.index()==node2().index();
    }

    /**return true if this edge has node n1 and node n2 */
    bool has_node(const Node& n1, const Node& n2) const{
      return has_node(n1)  &&  has_node(n2);
    }

    /** Return a node of this Edge O(1) */
    Node node1() const{
      return Node(mesh_,mesh_->edges_[edge_id_].node_id1);
    } 

    /** Return the other node of this Edge O(1) */
    Node node2() const{
      return Node(mesh_,mesh_->edges_[edge_id_].node_id2);
    }  

    /** Return an associated triangle of this Edge O(node1().num_triangles()) */
    Triangle triangle1() const{
      return Triangle(mesh_, mesh_->adj_e2t_[edge_id_][0]);
    }

    /** Return the other triangle of this Edge 
    return an invalid triangle if this is an boundary edge O(node1().num_triangles()) */
    Triangle triangle2() const {
      return Triangle(mesh_, mesh_->adj_e2t_[edge_id_][1]);
    }
    
    /** Test whether this edge and @a x are equal. O(1) */
    bool operator==(const Edge& x) const {
      if ((mesh_==x.graph_) && ((node1()==x.node1() && node2()==x.node2()) 
	                     || (node1()==x.node2() && node2()==x.node1())))
	return true;
      return false;
    }

    /** Test whether an edge is valid O(1) */
    bool valid() {
      return edge_id_<mesh_->edges_.size();
    }

    /** Test whether this edge is less than @a x in the global order O(1) */
    bool operator<(const Edge& x) const {
      if (mesh_ < x.mesh_){  //compare pointer address, this allows comparison between different graphs
        return true;
      }
      if ( (node1().index()+node2().index()) < (x.node1().index()+x.node2().index()) )  //compare index sum
	return true;
      else if ( (node1().index()+node2().index()) == (x.node1().index()+x.node2().index()) ) 
        if    ( (node1()<x.node1() && node1()<x.node2()) ||  (node2()<x.node1() && node2()<x.node2()) ) // one of the node is smallest
          return true;
      return false;
    }

    /** return degree (number of connecting triangles) of this edge O(1) */
    size_type num_triangles() const{
      return (mesh_->adj_e2t_[edge_id_][1] == UINT_MAX) ? 1 : 2;
    }

  private:
    //Allow Mesh to access Edge's private member data and functions
    friend class Mesh;

    //private constructor
    Edge(const Mesh* mesh, size_type edge_id)
      :mesh_(const_cast<Mesh*>(mesh)), edge_id_(edge_id){
    }

    size_type index() const {
      return edge_id_;
    }

    Mesh* mesh_;
    size_type edge_id_;
  };
  //end EDGE
  
  /** compute length between two nodes */
  double edge_length(const Node& n1, const Node& n2) const{
    return norm(n1.position()-n2.position());
  }  

  /** compute normal between two nodes 
   * the normal has norm of 1
   */
  Point edge_normal(const Node& n1, const Node& n2) const{
    return cross(Point(0,0,1),n2.position()-n1.position())/
      norm(cross(Point(0,0,1),n2.position()-n1.position()));
  }

  /** Return the number of edges in the mesh.O(1) */
  size_type num_edges() const {
    return edges_.size();
  }
 
  private:
  /*return the index of the edge if edge e is already contained 
   * O(1) */
  size_type has_edge(const Edge& e) const{
    return e.edge_id_;
  }

  /* return the index of the edge if an edge between n1 and n2 
   * is already contained 
   * O(n1.num_triangles()) */
  size_type has_edge(const Node& n1, const Node& n2) const{
    for (auto n2t_it = n1.incident_n2t_begin(); n2t_it!=n1.incident_n2t_end(); ++n2t_it){
        // return edge id
        for (int i = 0; i < 3; i++) {
          if ((*n2t_it).edge(i).has_node(n1, n2))
            return (*n2t_it).edge(i).index(); 
        }
    }
    return UINT_MAX;
  }

  /*add an edge to the graph
   * @return the newly added edge or the already added edge
   * O(n1.num_triangles()) */
  Edge add_edge(const Node& n1, const Node& n2, const edge_value_type& edge_value = edge_value_type()){
    //if already has an edge defined by n1, and n2, return an empty edge
    if (has_edge(n1,n2) != UINT_MAX){
      return Edge(this, has_edge(n1,n2));
    }
    
    // add new edge
    Point normal = edge_normal(n1,n2);
    double length = edge_length(n1,n2);
    edges_.push_back(edge_data(n1.index(),n2.index(),normal,length,edge_value));    
    return Edge(this,edges_.size()-1);
  }


public:
  ////////////////////
  // MESH TRIANGLES //
  ////////////////////
  
  class Triangle : private totally_ordered<Triangle> {
  public:
    /** constructor of an invalid triangle */
    Triangle() : triangle_id_(UINT_MAX) {}

    /**return index of the current triangle*/
    size_type index() const{
      return triangle_id_;
    }

    /**return area of the current triangle O(1) */
    double area() const {
      //return 0.5 *norm(cross(node(0).position()-node(1).position()), node(0).positon()-node(2).position()));
      return mesh_->triangles_[triangle_id_].area;
    }

    /**return node(n) of this triangle O(1) */
    Node node(size_type ind) const{
      return Node(mesh_,mesh_->triangles_[triangle_id_].node_ids[ind]);
    }

    /**return edge(n)of this triangle O(1) */
    Edge edge(size_type ind) const{
      return Edge(mesh_,mesh_->triangles_[triangle_id_].edge_ids[ind]);
    }

    /**return normal(n) of this triangle O(1) */
    Point normal(size_type ind) const{
      return mesh_->triangles_[triangle_id_].normal_signs[ind] *
             mesh_->edges_[mesh_->triangles_[triangle_id_].edge_ids[ind]].normal;
    }
    /**test whether this triangle and @a x are equal O(1) */
    bool operator==(const Triangle& x) const{
      return mesh_==x.mesh_ && triangle_id_==x.triangles_id_;
    }

    /** Test whether this edge is less than @a x in the global order. O(1)*/
    bool operator<(const Triangle& x) const{
      return triangle_id_ < x.triangle_id_;
    }

    /** return the edge_id_ of an edge of this triangle that has n1 and n2 
     * @pre the current triangle must contain n1 and n2
     */
    size_type edge_index(const Node& n1, const Node& n2) const{
      assert((*this).has_node(n1) && (*this).has_node(n2));
      //size_type edge_idx = 0-1; //deliberate underflow
      if (edge(0).has_node(n1,n2))
        return edge(0).edge_id_;
      if (edge(1).has_node(n1,n2))
        return edge(1).edge_id_;
      if (edge(2).has_node(n1,n2))
        return edge(2).edge_id_;
     assert(false);
    }

    /** test whether the current triangle has a node n */
    bool has_node(const Node& n) const{
      return mesh_->triangles_[triangle_id_].node_ids[0]==n.node_id_ ||
             mesh_->triangles_[triangle_id_].node_ids[1]==n.node_id_ ||
	     mesh_->triangles_[triangle_id_].node_ids[2]==n.node_id_;
    }

    /** test whether the current triangle has an edge e */
    bool has_edge(const Edge& e) const{
      return mesh_->triangles_[triangle_id_].edge_ids[0]==e.edge_id_ ||
             mesh_->triangles_[triangle_id_].edge_ids[1]==e.edge_id_ ||
	     mesh_->triangles_[triangle_id_].edge_ids[2]==e.edge_id_;
    }


    /** return degree (number of connecting triangles) of this triangle O(1) */
    size_type num_triangles() const{

      size_type num_tri = 0;
      for (size_type count = 0; count < 3; ++count){
        if (mesh_->adj_t2t_[triangle_id_][count] != UINT_MAX)
          num_tri++;
      }

      return num_tri; 
    }

    /**test whether this triangle is valid O(1) */
    bool valid() const{
      //std::cout<<"triangle_id is: "<< triangle_id_<< " size is: " << mesh_->triangles_.size() <<
      //",  number of triangles is: "<< num_triangles() <<std::endl;
      return triangle_id_ < mesh_->triangles_.size() && num_triangles()<=3;
    }

    /*return the reference to user defined value of this triangle O(1) */
    triangle_value_type& value(){
      return mesh_->triangles_[triangle_id_].triangle_value;
    }

    /*return the reference to user defined value of this triangle O(1) */
    const triangle_value_type& value() const{
      return triangles_[triangle_id_].triangle_value;
    }

    /*returns the neighboring triangles, invalid if has less than idx neighbor O(1)*/
    Triangle neighbor(size_type idx){
      return Triangle(mesh_, mesh_->adj_t2t_[triangle_id_][idx]);
    }

  private:
    //Allow Mesh to access Triangle's private member data and functions
    friend class Mesh;
    
    //private constructor
    Triangle(const Mesh* mesh, size_type triangle_id)
      :mesh_(const_cast<Mesh*>(mesh)), triangle_id_(triangle_id){
    }

    //private constructor
    //Triangle(const Mesh* mesh, Node n1, Node n2, Node n3)
    //  : mesh_(const_cast<Mesh*>(mesh)) {
    //
    //}

    Mesh* mesh_;           // a pointer back to mesh
    size_type triangle_id_; // id of this triangle

  };
  //end Triangle 
  
  /*return the area of a triangle consisted by 3 points.
   * @pre the triangle does not need to be valid 
   * O(1) */
  double get_area(const Node& n1, const Node& n2, const Node& n3) const {
    return 0.5 * norm(cross(n1.position()-n2.position(), n1.position()-n3.position()));
  }

  /*returns total number of triangles 
   * O(1) */
  size_type num_triangles() const {
    return triangles_.size();
  }
 
  /*return true if triangle t is already contained 
   * O(1) */
  bool has_triangle(const Triangle& t) const{
    return t.triangle_id_ < triangles_.size();
  }

  /*return true if a triangle defined by the three Nodes is already contained 
   * O(n1.num_triangles()) */
  size_type has_triangle(const Node& n1, const Node& n2, const Node& n3) const{
    for (auto n2t_it = n1.incident_n2t_begin(); n2t_it!=n1.incident_n2t_end(); ++n2t_it){
      if ((*n2t_it).has_node(n2) && (*n2t_it).has_node(n3))
        return (*n2t_it).index();
    }
    return UINT_MAX;
  }

  /*add a triangle to the graph
   * overall: O(n1.num_triangles()+n2.num_triangles()+n3.num_triangles())
   * 1. check has_triangle: O(n1.num_triangles())
   * 2. update edges_ : O(n1.num_triangles())
   * 3. update triangles_: O(1)
   * 4  update adj_e2t_: O(1)
   * 5. update adj_n2t_: O(1)
   * 6. update adj_t2t_: O(n1.num_triangles()+n2.num_triangles()+n3.num_triangles())
   */
  Triangle add_triangle(const Node& n1, const Node& n2, const Node & n3,
                    const triangle_value_type& triangle_value = triangle_value_type()){
    

    
    // if this triangle is included already, return it
    if (has_triangle(n1,n2,n3) != UINT_MAX){
          return Triangle(this, has_triangle(n1,n2,n3));
    }
    
    std::vector<Edge> edges = {add_edge(n1,n2), add_edge(n2, n3), 
                               add_edge(n1, n3)};
   
    update_triangles(n1, n2, n3, edges, triangle_value); 
   
    update_adj_e2t(edges);
    
    //update adj_n2t_
    adj_n2t_[n1.index()].push_back(triangles_.size()-1);
    adj_n2t_[n2.index()].push_back(triangles_.size()-1);
    adj_n2t_[n3.index()].push_back(triangles_.size()-1);

    //update adj_t2t_
    update_adj_t2t(edges, triangles_.size() - 1);


    //return the new triangle
    return Triangle(this,triangles_.size()-1);
  }

  private:

  void update_triangles(const Node& n1, const Node& n2, const Node& n3, const std::vector<Edge> edges, const triangle_value_type& triangle_value) {

    // precompute normal values for triangle
    double area= get_area(n1,n2,n3);
    unsigned_triplet node_ids = {{n1.index(), n2.index(), n3.index()}};
    unsigned_triplet edge_ids = {{edges[0].index(), edges[1].index(), 
                                  edges[2].index()}};
    // find whether the outward normal to this triangle is +1 or -1 times
    // the normal stored in each of its edges 
    Point centroid = (n1.position()+n2.position()+n3.position())/3;
    int_triplet normal_signs;
    for (int i = 0; i < 3; i++) {
 
      Point midpoint = (edges[i].node1().position() + 
                        edges[i].node2().position() ) / 2.0;
      
      normal_signs[i] = (dot(edges[i].normal(), centroid-midpoint)<0) ? 1 : -1;
    }

    //update triangles_
    triangles_.push_back({normal_signs,node_ids,edge_ids,area,triangle_value});
   
   }

  void update_adj_e2t(const std::vector<Edge> edges) {
    for (int i = 0; i < 3; i++) {

      if (adj_e2t_.size() <= edges[i].index()) {
        unsigned_doublet new_edge_adj = {{num_triangles() - 1, UINT_MAX}};
        adj_e2t_.push_back(new_edge_adj);
      }

      else
        adj_e2t_[edges[i].index()][1] = triangles_.size() - 1;
    }


  }
  
  /*update adj_t2t_ */
  void update_adj_t2t(std::vector<Edge> edges, size_type self_id) {
    assert(is_valid_tri_id(self_id));
    assert(adj_t2t_.size() == self_id);

    unsigned_triplet neighbor_tids;
    for (int i = 0; i < 3; ++i)
      neighbor_tids[i] = (adj_e2t_[edges[i].index()][0] == triangles_.size()-1) 
                          ? UINT_MAX : adj_e2t_[edges[i].index()][0];
  
    //add in last element
    adj_t2t_.push_back(neighbor_tids);

    Triangle new_tri = Triangle(this, self_id);
    
    for (int i = 0; i < 3; ++i) {
      if (is_valid_tri_id(neighbor_tids[i])) {

        Triangle neighbor_tri = Triangle(this, neighbor_tids[i]);

        for (int j = 0; j < 3; ++j)
          if (new_tri.has_edge(neighbor_tri.edge(j)))
            adj_t2t_[neighbor_tids[i]][j] = self_id;
        }
      } 
    }
  
  /* check if the edge_id is valid */
  bool is_valid_tri_id(size_type tri_id) const{
    //std::cout<<"tri id: " << tri_id << " triangles_.size is: " << triangles_.size() << std::endl;
    return tri_id < triangles_.size();
  }
 
  
  public:
  ///////////////////
  // NODE ITERATOR //
  ///////////////////
  
  /* @brief Iterator class for nodes. a forward iterator*/
  class NodeIterator: private totally_ordered<NodeIterator>{
  public:
    // These type definitions help us use STL's iterator_traits. 
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid NodeIterator. O(1) */
    NodeIterator(){
    }

    /**node dereferencing operator O(1)*/
    Node operator*() const{
      return Node(mesh_,node_id_);
    }

    /**iterator incrementing operator O(1) */
    NodeIterator& operator++(){
      node_id_++;
      return *this;
    }

    /** iterator comparison operator O(1) */
    bool operator==(const NodeIterator& node_iter) const{
      return mesh_== node_iter.mesh_ && node_id_ == node_iter.node_id_;
    }

  private:
    //constructor O(1)
    NodeIterator(const Mesh* mesh, size_type node_id)
      : mesh_(const_cast<Mesh*>(mesh)), node_id_(node_id){
    }

    friend class Mesh;
    Mesh* mesh_;
    size_type node_id_;  
  };
  //end NODE ITERATOR  

  //returns a node_iterator of first node in graph O(1)
  NodeIterator node_begin() const{
    return NodeIterator(this,0);
  }

  //returns an invalid node_iterator O(1)
  NodeIterator node_end() const{
    return NodeIterator(this, nodes_.size());
  }


  ///////////////////
  // EDGE ITERATOR //
  ///////////////////
  /* @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
    public:
      // These type definitions help us use STL's iterator_traits.
      /** Element type. */
      typedef Edge value_type;
      /** Type of pointers to elements. */
      typedef Edge* pointer;
      /** Type of references to elements. */
      typedef Edge& reference;
      /** Iterator category. */
      typedef std::input_iterator_tag iterator_category;
      /** Difference between iterators */
      typedef std::ptrdiff_t difference_type;

      /** Construct an invalid EdgeIterator. O(1) */
      EdgeIterator() {
      }

      /* dereference operator for an edge_iterator O(1) */
      Edge operator*(){
	return Edge(mesh_,edge_id_);
      }

      /* increment operator for an edge_iterator O(1) */
      EdgeIterator& operator++(){
        edge_id_++;
        return *this;
      }

      /* comparison operator for equality O(1) */
      bool operator == (const EdgeIterator& e_it) const{
        return mesh_==e_it.mesh_ && edge_id_ == e_it.edge_id_;
      }

    private:
      //constructor O(1) 
      EdgeIterator(const Mesh* mesh,size_type edge_id)
        :mesh_(const_cast<Mesh*>(mesh)), edge_id_(edge_id){
      }  

      friend class Mesh;
      Mesh* mesh_;
      size_type edge_id_; 
  };
  //end EdgeIterator

  //returns a edge_iterator of first node in graph O(1) 
  EdgeIterator edge_begin() const{
    return EdgeIterator(this,0);
  }

  //returns an invalid edge_iterator O(1) 
  EdgeIterator edge_end() const{
    return EdgeIterator(this,edges_.size());
  }

  ///////////////////////
  // TRIANGLE ITERATOR //
  ///////////////////////
  /* @brief Iterator class for triangles. A forward iterator. */
  class TriangleIterator : private totally_ordered<TriangleIterator> {
    public:
      // These type definitions help us use STL's iterator_traits.
      /** Element type. */
      typedef Triangle value_type;
      /** Type of pointers to elements. */
      typedef Triangle* pointer;
      /** Type of references to elements. */
      typedef Triangle& reference;
      /** Iterator category. */
      typedef std::input_iterator_tag iterator_category;
      /** Difference between iterators */
      typedef std::ptrdiff_t difference_type;

      /** Construct an invalid EdgeIterator. O(1) */
      TriangleIterator(){
      }

      /* dereference operator for an triangle_iterator O(1)  */
      Triangle operator*(){
        return Triangle(mesh_,triangle_id_);
      }

      /* increment operator for an triangle_iterator O(1) */
      TriangleIterator& operator++(){
        triangle_id_++;
        return *this;
      }

      /* comparison operator for equality O(1) */
      bool operator == (const TriangleIterator& t_it) const{
        return mesh_==t_it.mesh_ && triangle_id_==t_it.triangle_id_;
      }

    private:
      //constructor O(1) 
      TriangleIterator(const Mesh* mesh,size_type triangle_id):
        mesh_(const_cast<Mesh*>(mesh)),triangle_id_(triangle_id){
      }  

      friend class Mesh;
      Mesh* mesh_;
      size_type triangle_id_; 
  };
  //end TriangleIterator

  //returns a edge_iterator of first node in graph  O(1) 
  TriangleIterator triangle_begin() const{
    return TriangleIterator(this,0);
  }

  //returns an invalid edge_iterator O(1) 
  TriangleIterator triangle_end() const{
      return TriangleIterator(this,triangles_.size());  
  }

  ///////////////////////////
  // incident_n2t_iterator //
  ///////////////////////////

  /* @brief Iterator class for triangles incident to a node. A forward iterator. */
  class incident_n2t_iterator :private totally_ordered<incident_n2t_iterator> {
  public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Triangle value_type;
    /** Type of pointers to elements. */
    typedef Triangle* pointer;
    /** Type of references to elements. */
    typedef Triangle& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid incident_n2t_iterator O(1) */
    incident_n2t_iterator() {
    }
   
    /* dereference operator of  incident_n2t_iterator O(1) */
    Triangle operator*(){
      return Triangle(mesh_,mesh_->adj_n2t_[node_id_][col_idx_]);
    }

    /* increment operator of incident_n2t_iterator O(1) */
    incident_n2t_iterator& operator++(){
      col_idx_++;
      return *this;
    }

    /* comparison operator of incident_n2t_iterator O(1) */
    bool operator==(const incident_n2t_iterator& it) const{
      return mesh_== it.mesh_ && node_id_ == it.node_id_ && col_idx_==it.col_idx_ ;
    }

    //constructor O(1) 
    incident_n2t_iterator(const Mesh* mesh,size_type nid, size_type cid)
      :mesh_(const_cast<Mesh*>(mesh)),node_id_(nid),col_idx_(cid){
    }
  private:
    friend class Graph;
    Mesh* mesh_;
    size_type node_id_;
    size_type col_idx_;
  };
  // end incident_n2t_iterator

  //helper debug functions
  void print_nodes(){
    std::cout<<"nodes information is: "<<std::endl;
    for (size_type i = 0; i<nodes_.size(); ++i){
      std::cout<< "node index: " << i <<", position: " << nodes_[i].point.x 
               << "     " << nodes_[i].point.y << "     " << nodes_[i].point.z <<std::endl;
    }
  }

  void print_edges(){
    std::cout<<"edges information is: "<<std::endl;
    for (size_type i = 0; i<edges_.size(); ++i){
      std::cout<< "edge index: " << i <<", node index: " << edges_[i].node_id1 
               << "     " << edges_[i].node_id2  << "     " <<std::endl;
    }
  }

  void print_triangles(){
    std::cout<<"triangles information is: "<<std::endl;
    for (size_type i = 0; i<triangles_.size(); ++i){
      std::cout<< "triangle index: " << i <<", node index: " << triangles_[i].node_ids[0] 
               << "     " << triangles_[i].node_ids[1]  << "     " 
               << triangles_[i].node_ids[2]  << "  area is: "<<triangles_[i].area << "  "
               << "number of neighbors is: " << Triangle(this,i).num_triangles()<<std::endl;

      std::cout<< "triangle normal: " << " ("<<Triangle(this,i).normal(0).x <<" , " << Triangle(this,i).normal(0).y 
	       << ") , ( " << Triangle(this,i).normal(1).x <<" , " << Triangle(this,i).normal(1).y
               << ") , ( " << Triangle(this,i).normal(2).x <<" , " << Triangle(this,i).normal(2).y << ") " << std::endl;
    }
  }

  void print_n2t(){
    std::cout<<"adj_n2t information is: "<<std::endl;
    for (size_type i = 0; i<nodes_.size(); ++i){
      std::cout<< "node index: " << i << ", triangle index:  ";
      for (auto it =Node(this,i).incident_n2t_begin(); it!=Node(this,i).incident_n2t_end(); ++it){
	std::cout<<(*it).triangle_id_ << "   ";
      }
      std::cout<<std::endl;
    }
  }

  void print_t2t(){
    std::cout<<"adj_t2t information is: "<<std::endl;
    for (size_type i = 0; i<adj_t2t_.size(); ++i){
      std::cout<< "triangle index: " << i << ", triangle index:  "
               << adj_t2t_[i][0] << "   " << adj_t2t_[i][1] << "   "
               << adj_t2t_[i][2] << std::endl;
    }
  }

//mesh's private data members
private:
  //friend class Node;
  //friend class Edge;
  //friend class Triangle;
  
  struct node_data{
    Point point;
    node_value_type node_value;

    //constructor
    node_data(Point p, node_value_type v):
      point(p), node_value(v){
    }
  };

  struct edge_data{
    size_type node_id1;
    size_type node_id2;
    Point normal;
    double edge_length;
    edge_value_type edge_value;

    //constructor
    edge_data(size_type n1, size_type n2, Point nor, double edge_len, edge_value_type v)
      :node_id1(n1), node_id2(n2), normal(nor), edge_length(edge_len), edge_value(v) {
    } 
  };

  struct triangle_data{
    int_triplet normal_signs;   //either +1 or -1
    unsigned_triplet node_ids;
    unsigned_triplet edge_ids; 
    double area;
    triangle_value_type triangle_value;

    //constructor 
    triangle_data(int_triplet nor_signs, unsigned_triplet n_ids, unsigned_triplet e_ids, double a,  triangle_value_type v)
      :normal_signs(nor_signs), node_ids(n_ids), edge_ids(e_ids), area(a), triangle_value(v) {
    } 
  };
  std::vector<node_data> nodes_;
  std::vector<edge_data> edges_;
  std::vector<triangle_data> triangles_;

  std::vector<unsigned_doublet> adj_e2t_;

  std::vector<std::vector<size_type>> adj_n2t_; //adjacency list of incident triangles to node
  std::vector<unsigned_triplet> adj_t2t_; //adjacency list of incident triangles to triangle

};
