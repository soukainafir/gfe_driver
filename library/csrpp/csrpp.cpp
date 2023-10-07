/**
 * Copyright (C) 2019 Dean De Leo, email: dleo[at]cwi.nl
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "csrpp.hpp"

#include <cassert>
#include <mutex>

using namespace common;
using namespace libcuckoo;
using namespace std;
#define STINGER reinterpret_cast<struct stinger*>(m_stinger_graph)
#define INT2DBL(v) ( *(reinterpret_cast<double*>(&(v))) )
#define DBL2INT(v) ( *(reinterpret_cast<int64_t*>(&(v))) )

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
extern mutex _log_mutex [[maybe_unused]];
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{::_log_mutex}; std::cout << "[CSRPP::" << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
#define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
#define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 *  Error                                                                    *
 *                                                                           *
 *****************************************************************************/
//#undef CURRENT_ERROR_TYPE

/*****************************************************************************
 *                                                                           *
 *  Initialisation                                                           *
 *                                                                           *
 *****************************************************************************/
namespace gfe::library {
    
    CSRPP::CSRPP(bool directed) : m_directed(directed) {                 
          G_MM = new graph_csrpp();
    }

    
    CSRPP::~CSRPP(){
        G_MM = nullptr;
    } 
    
    void CSRPP::preallocate_vertices(size_t num_segments) {
        G_MM->map= std::make_unique<csrpp>(num_segments);
        if(m_directed) G_MM->r_map= std::make_unique<csrpp>(num_segments); 
    }


    void CSRPP::print_lock_stat(){
#ifdef LOCK_IN 
        for (uint32_t seg_id = 0; seg_id <  G_MM->map->get_num_segments(); seg_id++) {
          auto segment = G_MM->map->get_segment(seg_id);
       
          segment.print_lock_stat();
     
    }
#endif
    }
/******************************************************************************
 *                                                                            *
 *  Conversion functions                                                      *
 *                                                                            *
 *****************************************************************************/


/******************************************************************************
 *                                                                            *
 *  Properties                                                                *
 *                                                                            *
 *****************************************************************************/
    uint64_t CSRPP::num_edges() const {
        uint64_t num_directed_edges = static_cast<uint64_t>(G_MM->map->get_total_num_edges());
        if(!m_directed) num_directed_edges /= 2; /* divided by two, because we always store the edge for both source/destination */
        return num_directed_edges;
    }

    uint64_t CSRPP::num_vertices() const {
        return G_MM->map->get_num_keys();
       // return G_MM->mapping.size();
    }

    bool CSRPP::has_vertex(uint64_t external_vertex_id) const {
        vertex_t vertex;
        return G_MM->get_node_key(uint64_to_nodekey(external_vertex_id), vertex);
    }

    double CSRPP::get_weight(uint64_t ext_source_id, uint64_t ext_destination_id) const {
        COUT_DEBUG("external: " << ext_source_id << " -> " << ext_destination_id);
        constexpr double NaN { numeric_limits<double>::signaling_NaN() };

         vertex_t src_vertex, dst_vertex;
        bool src_exists = G_MM->get_node_key(G_MM->int_to_nodekey(ext_source_id), src_vertex);

        if (!src_exists) {  return NaN; }
        bool dst_exists = G_MM->get_node_key(G_MM->int_to_nodekey(ext_destination_id), dst_vertex);

        if (!dst_exists) {  return NaN; }

        int32_t out_degree;
        vertex_t* edges = G_MM->map->get_neighbors(src_vertex, out_degree);
        if(out_degree == 1) {
          return G_MM->map->get_edge_property_for_vertex<double>(src_vertex, 0, 0);
        } else if (out_degree > 1){
          for(size_t i = 0; i < out_degree; i++) {
            if(edges[i] == dst_vertex) {
              return G_MM->map->get_edge_property_for_vertex<double>(src_vertex, i, 0);
            } 
             // 
            
          } 
          printf("error \n");
          //size_t index = upper_bound(edges, edges + out_degree, dst_vertex) - edges;
          //if(index != 1) return G_MM->map->get_edge_property_for_vertex<double>(src_vertex, index, 0);
        }
        
        return NaN;

    }

    bool CSRPP::is_directed() const {
        return m_directed;
    }

    void CSRPP::set_timeout(uint64_t seconds) {
        m_timeout = seconds;
    }

    void* CSRPP::handle(){
        return G_MM;
    }

/******************************************************************************
 *                                                                            *
 *  Updates                                                                   *
 *                                                                            *
 *****************************************************************************/
    bool CSRPP::add_vertex(uint64_t vertex_id){
        COUT_DEBUG("vertex_id: " << vertex_id);
        G_MM->add_if_not_exist(vertex_id);
        return true;
    }

bool CSRPP::remove_vertex(uint64_t external_vertex_id){
        COUT_DEBUG("vertex_id: " << external_vertex_id);
        return G_MM->delete_vertex(external_vertex_id);
}

    bool CSRPP::add_edge(graph::WeightedEdge e){
        COUT_DEBUG("edge: " << e);
        
      /*  decltype(G_MM->m_vmap)::const_accessor accessor1, accessor2; // shared lock on the dictionary
        if(!G_MM->m_vmap.find(accessor1, e.source())){ return false; }
        if(!G_MM->m_vmap.find(accessor2, e.destination())) { return false; }
*/
        vertex_t src_vertex, dst_vertex;
        bool src_exists = G_MM->get_node_key(G_MM->int_to_nodekey(e.source()), src_vertex);

        if (!src_exists) {  return false; }
        bool dst_exists = G_MM->get_node_key(G_MM->int_to_nodekey(e.destination()), dst_vertex);

        if (!dst_exists) {  return false; }

       
        // get the indices in the map
        int64_t weight = DBL2INT(e.m_weight);
       /* #pragma omp parallel for 
            for (node_t t = 0; t < 100000; t++) {
                G_MM->map->add(v, v2);
            }
        */
        G_MM->map->add_with_weight<double>(src_vertex, dst_vertex, e.m_weight);
        
        if(m_directed) {
            G_MM->r_map->add(dst_vertex, src_vertex);
        } else{
            G_MM->map->add_with_weight<double>(dst_vertex, src_vertex, e.m_weight);
        }
        
       //G_MM->r_map->add(dst_vertex, src_vertex);

        //add edge prop insertion in direct map like stinger


        return true;
    }

    bool CSRPP::add_edge_v2(gfe::graph::WeightedEdge e){
        COUT_DEBUG("edge: " << e);
        vertex_t src_vertex, dst_vertex;
        
        
        
        // get the indices in the adjacency lists
        
       /* decltype(G_MM->m_vmap)::const_accessor accessor1, accessor2; // shared lock on the dictionary
        if(!G_MM->m_vmap.find(accessor1, e.source())){
        //    global_lock.lock();
            src_vertex = G_MM->add_vertex(e.source());
        //    global_lock.unlock();
        } else src_vertex = accessor1->second;
 
        if(!G_MM->m_vmap.find(accessor2, e.destination())) {
         //   global_lock.lock();
            dst_vertex = G_MM->add_vertex(e.destination());
        //    global_lock.unlock();
        } else dst_vertex = accessor2->second;
        */
     /*   bool exists = get_node_key(e.source(), src_vertex);       
        if (!exists) {
          src_vertex = map->add_vertex_rr(0);
          r_map->add_vertex_rr(0);              
        }
        bool exists = get_node_key(e.destination(), dst_vertex);       
        if (!exists) {
          dst_vertex = map->add_vertex_rr(0);
          r_map->add_vertex_rr(0);              
        }*/
        src_vertex = G_MM->add_if_not_exist(e.source());
        dst_vertex = G_MM->add_if_not_exist(e.destination());
        
        /*  bool src_exists = G_MM->get_node_key(uint32_to_nodekey(e.source()), src_vertex);
          bool dst_exists = G_MM->get_node_key(uint32_to_nodekey(e.destination()), dst_vertex);
          {
          scoped_lock<SpinLock> lock(m_mutex_vtx);
            
          if (!src_exists ) {
            src_vertex = G_MM->add_vertex(uint32_to_nodekey(e.source()));
          }
          if (!dst_exists) {
            
            dst_vertex = G_MM->add_vertex(uint32_to_nodekey(e.destination()));
          }
          }*/
        
        
       // printf("num vertices = %d\n", num_vertices());
        
        int64_t weight = DBL2INT(e.m_weight);
        
         G_MM->map->add_with_weight<double>(src_vertex, dst_vertex, e.m_weight);
         if(m_directed) {
            G_MM->r_map->add(dst_vertex, src_vertex);
        } else{
            G_MM->map->add_with_weight<double>(dst_vertex, src_vertex, e.m_weight);
        }
        return true;
    }

    bool CSRPP::remove_edge(graph::Edge e){
         COUT_DEBUG("edge: " << e);
        vertex_t src_vertex, dst_vertex;
        bool src_exists = G_MM->get_node_key(G_MM->int_to_nodekey(e.source()), src_vertex);

        if (!src_exists) {  return false; }
        bool dst_exists = G_MM->get_node_key(G_MM->int_to_nodekey(e.destination()), dst_vertex);

        if (!dst_exists) {  return false; }
         G_MM->map->remove(src_vertex, dst_vertex, false);
         G_MM->map->remove(dst_vertex, src_vertex, false);
        return true;  
    }
    
    //GraphAlytics
    
      /* Perform a BFS from source_vertex_id to all the other vertices in the graph.
         * @param source_vertex_id the vertex where to start the search
         * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
         */
         void CSRPP::bfs(uint64_t source_vertex_id, const char* dump2file){}

        /**
         * Execute the PageRank algorithm for the specified number of iterations.
         *
         * @param num_iterations the number of iterations to execute the algorithm
         * @param damping_factor weight for the PageRank algorithm, it affects the score associated to the sink nodes in the graphs
         * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
         */
        //  void CSRPP::pagerank(uint64_t num_iterations, double damping_factor, const char* dump2file = nullptr){}

        /**
         * Weakly connected components (WCC), associate each node to a connected component of the graph
         * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
         */
         void CSRPP::wcc(const char* dump2file){}

        /**
         * Community Detection using Label-Propagation. Associate a label to each vertex of the graph, according to its neighbours.
         * @param max_iterations max number of iterations to perform
         * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
         */
         void CSRPP::cdlp(uint64_t max_iterations, const char* dump2file){}

        /**
         * Local clustering coefficient. Associate to each vertex the ratio between the number of its outgoing edges and the number of
         * possible remaining edges.
         * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
         */
         void CSRPP::lcc(const char* dump2file){}

        /**
         * Single-source shortest paths. Compute the weight related to the shortest path from the source to any other vertex in the graph.
         * @param source_vertex_id the vertex where to start the search
         * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
         */
         void CSRPP::sssp(uint64_t source_vertex_id, const char* dump2file){}

        /**
         * Retrieve the internal handle to the library implementation
         */
        
         void CSRPP::load(const std::string& path){}

/******************************************************************************
 *                                                                            *
 *  Dump                                                                      *
 *                                                                            *
 *****************************************************************************/

    void CSRPP::dump_ostream(std::ostream& out) const {
       // struct CSRPP* graph = STINGER;

        out << "[CSRPP] Vertices: " << num_vertices() << ", edges: " << num_edges() << ", directed: " << std::boolalpha << is_directed() << ", size: NOT YET IMPLEMENTED bytes" << "\n";
        for (uint32_t seg_id = 0; seg_id < G_MM->map->get_num_segments();seg_id++) {
            for (uint32_t v_id = 0; v_id < G_MM->map->NUM_VERTICES_PER_SEGMENT; v_id++) {
                auto &segment = G_MM->map->get_segment(seg_id);
                vertex_info_t info = segment.at(v_id);
                int32_t &out_degree = info.length;
                if (out_degree == -1) continue;
                vertex_t vertex(seg_id, v_id);
                //print external IDs
                uint64_t external_vertex_id = G_MM->map->get_vertex_property<long>(vertex, 0);
                if(external_vertex_id != numeric_limits<uint64_t>::max()){
                    out << "[" << external_vertex_id << " (internal ID: " << G_MM->map->get_absolute_vertex_id(vertex) << ")] ";
                } else {
                    out << "[" << external_vertex_id << "] ";
                }
                //print the outgoing edges
                out << " weight: Not Yet Implemented " << "degree out: " << out_degree;
                if(out_degree == 0){
                    out << "\n";
                } else {
                    out << ", outgoing edges: \n";
                    vertex_t* neighbors = info.get_neighbors();
                    for (edge_t w_idx = 0; w_idx < out_degree; w_idx ++)
                    {
                        vertex_t &w =  neighbors[w_idx];
                        uint64_t vertex_id_name = G_MM->map->get_vertex_property<long>(w, 0);;
                        if(vertex_id_name != numeric_limits<uint64_t>::max()){
                            out << vertex_id_name << ", internal ID: " << G_MM->map->get_absolute_vertex_id(w);
                        } else {
                            out << G_MM->map->get_absolute_vertex_id(w);
                        }
                        out << ", ";

                        /*out << "type: " << dump_edge_type(graph, STINGER_EDGE_TYPE) << ", " <<
                            "weight: Not Yet Implemented" << "\n";*/
                    }
                }
            }
        }
        std::flush(out);
    }


} // namespace
