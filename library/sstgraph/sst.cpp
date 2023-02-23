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

#include "sst_internal.hpp"
#include "sst.hpp"
#if defined(LLAMA_HASHMAP_WITH_TBB)
#include "tbb/concurrent_hash_map.h"
#else
#include "third-party/libcuckoo/cuckoohash_map.hh"
#endif

#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <fstream>
#include <limits>
#include <mutex>
#include <queue>
#include <thread>
#include <map>
#include "common/time.hpp"
#include <cassert>
#include <mutex>
#include "third-party/gapbs/gapbs.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>
#include <cmath>
#include <omp.h>


#include "common/system.hpp"
#include "third-party/libcuckoo/cuckoohash_map.hh"
#include "utility/timeout_service.hpp"

using namespace common;
using namespace gfe::utility;
using namespace libcuckoo;
using namespace std;

namespace gapbs { class Bitmap; }
namespace gapbs { template <typename T> class SlidingQueue; }
namespace gapbs { template <typename T> class pvector; }

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

// Macros
#define TIMER_INIT auto time_start = chrono::steady_clock::now();
#define CHECK_TIMEOUT if(m_timeout > 0 && (chrono::steady_clock::now() - time_start) > chrono::seconds(m_timeout)) { \
        RAISE_EXCEPTION(TimeoutError, "Timeout occurred after: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - time_start).count() << " seconds") };


/******************************************************************************
 *                                                                            *
 *  Utility functions                                                         *
 *                                                                            *
 ******************************************************************************/
// dump the content to the given file
static void save(cuckoohash_map<uint64_t, double>& result, const char* dump2file){
    if(dump2file == nullptr) return; // nop
    COUT_DEBUG("save the results to: " << dump2file)

    fstream handle(dump2file, ios_base::out);
    if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

    auto list_entries = result.lock_table();

    for(const auto& p : list_entries){
        handle << p.first << " " << p.second << "\n";
    }

    list_entries.unlock();
    handle.close();
}
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

    SSTGraph::SSTGraph(bool directed) : m_directed(directed) {
        uint32_t num_nodes = 0;
        uint64_t num_edges = 0;
        SparseMatrixV<true, uint32_t> g(num_nodes, num_nodes);

    }

    SSTGraph::SSTGraph(bool directed, uint32_t num_vertices) : m_directed(directed), m_num_vertices(num_vertices){
        g = new SparseMatrixV<true, uint32_t>(num_vertices, num_vertices);
    }


    SSTGraph::~SSTGraph(){
        g = nullptr;
    }

/******************************************************************************
 *                                                                            *
 *  Conversion functions                                                      *
 *                                                                            *
 *****************************************************************************/

    uint64_t SSTGraph::get_internal_vertex_id(uint64_t external_vertex_id) const {
#if defined(SST_HASHMAP_WITH_TBB)
        decltype(m_vmap)::const_accessor accessor;
    if ( m_vmap.find(accessor, external_vertex_id ) ){
        return accessor->second;
    } else {
        ERROR("The given vertex does not exist: " << external_vertex_id);
    }
#else
    return 0;
#endif
    }
/******************************************************************************
 *                                                                            *
 *  Properties                                                                *
 *                                                                            *
 *****************************************************************************/
    uint64_t SSTGraph::num_edges() const {
       // uint64_t num_directed_edges = static_cast<uint64_t>(G_MM->map->get_total_num_edges());
       // if(!m_directed) num_directed_edges /= 2; /* divided by two, because we always store the edge for both source/destination */
       // return num_directed_edges;
       return 0;
    }

    uint64_t SSTGraph::num_vertices() const {
        return g->get_rows();
    }

    bool SSTGraph::has_vertex(uint64_t external_vertex_id) const {
        uint32_t sst_source_vertex_id = get_internal_vertex_id(external_vertex_id);
        return g->has(sst_source_vertex_id, sst_source_vertex_id);
    }

    double SSTGraph::get_weight(uint64_t ext_source_id, uint64_t ext_destination_id) const {
       /* COUT_DEBUG("external: " << ext_source_id << " -> " << ext_destination_id);
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
            for(size_t i = 0; i < out_degree; i) {
                if(edges[i] == dst_vertex) {
                    return G_MM->map->get_edge_property_for_vertex<double>(src_vertex, i, 0);
                }
                //

            }
            printf("error \n");
            //size_t index = upper_bound(edges, edges  out_degree, dst_vertex) - edges;
            //if(index != 1) return G_MM->map->get_edge_property_for_vertex<double>(src_vertex, index, 0);
        }
*/
        return 0;

    }

    bool SSTGraph::is_directed() const {
        return m_directed;
    }

    void SSTGraph::set_timeout(uint64_t seconds) {
        m_timeout = seconds;
    }

    void* SSTGraph::handle(){
        return g;
    }

/******************************************************************************
 *                                                                            *
 *  Updates                                                                   *
 *                                                                            *
 *****************************************************************************/
    bool SSTGraph::add_vertex(uint64_t vertex_id) {
        COUT_DEBUG("vertex_id: " << vertex_id);
        uint32_t sst_source_vertex_id = get_internal_vertex_id(vertex_id);
        uint32_t value = 0;
        if (g->has(sst_source_vertex_id, sst_source_vertex_id)) return false;
       else {
           /* if (value > UINT64_MAX) {
                if constexpr (std::is
           _integral<int64_t>::value) {
                    value = value % UINT64_MAX;
                } else {
                    value = UINT64_MAX;
                }
            }*/
            g->insert(sst_source_vertex_id, sst_source_vertex_id, 0);
            return true;
        }

    }

    bool SSTGraph::remove_vertex(uint64_t external_vertex_id){
        COUT_DEBUG("vertex_id: " << external_vertex_id);
        auto sst_source_vertex_id = get_internal_vertex_id(external_vertex_id);
        g->remove(sst_source_vertex_id, sst_source_vertex_id);

            if (g->has(sst_source_vertex_id, sst_source_vertex_id)) {
                COUT_DEBUG("have something we removed while removing elements");
                return false;
            }

        return true;
    }

    //
    bool SSTGraph::add_edge(graph::WeightedEdge e){
        COUT_DEBUG("edge: " << e);
        auto sst_source_vertex_id = get_internal_vertex_id(e.source());
        auto sst_destination_vertex_id = get_internal_vertex_id(e.destination());

        // get the indices in the map
        int64_t weight = DBL2INT(e.m_weight);

        g->insert(sst_source_vertex_id, sst_destination_vertex_id, weight);

        return true;
    }

    bool SSTGraph::add_edge_v2(gfe::graph::WeightedEdge e){
        COUT_DEBUG("edge: " << e);
        int64_t sst_source_vertex_id;
        int64_t sst_destination_vertex_id;
        try {
             sst_source_vertex_id = get_internal_vertex_id(e.source());

        } catch(...) { // the vertex does not exist
           g->insert(sst_source_vertex_id, sst_source_vertex_id, 0);
        }

        try {
            sst_destination_vertex_id = get_internal_vertex_id(e.destination());

        } catch(...) { // the vertex does not exist
            g->insert(sst_destination_vertex_id, sst_destination_vertex_id, 0);
        }

        int64_t weight = DBL2INT(e.m_weight);
        g->insert(sst_source_vertex_id, sst_destination_vertex_id, weight);
        return true;
    }

    bool SSTGraph::remove_edge(graph::Edge e){
        /*COUT_DEBUG("edge: " << e);
        vertex_t src_vertex, dst_vertex;
        bool src_exists = G_MM->get_node_key(G_MM->int_to_nodekey(e.source()), src_vertex);

        if (!src_exists) {  return false; }
        bool dst_exists = G_MM->get_node_key(G_MM->int_to_nodekey(e.destination()), dst_vertex);

        if (!dst_exists) {  return false; }
        G_MM->map->remove(src_vertex, dst_vertex, false);
        G_MM->map->remove(dst_vertex, src_vertex, false);*/
        return true;
    }

    //GraphAlytics

    /* Perform a BFS from source_vertex_id to all the other vertices in the graph.
       * @param source_vertex_id the vertex where to start the search
       * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
       */
    void SSTGraph::bfs(uint64_t source_vertex_id, const char* dump2file){}

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
    void SSTGraph::wcc(const char* dump2file){}

    /**
     * Community Detection using Label-Propagation. Associate a label to each vertex of the graph, according to its neighbours.
     * @param max_iterations max number of iterations to perform
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    void SSTGraph::cdlp(uint64_t max_iterations, const char* dump2file){}

    /**
     * Local clustering coefficient. Associate to each vertex the ratio between the number of its outgoing edges and the number of
     * possible remaining edges.
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    void SSTGraph::lcc(const char* dump2file){}

    /**
     * Single-source shortest paths. Compute the weight related to the shortest path from the source to any other vertex in the graph.
     * @param source_vertex_id the vertex where to start the search
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    void SSTGraph::sssp(uint64_t source_vertex_id, const char* dump2file){}

    /**
     * Retrieve the internal handle to the library implementation
     */

    void SSTGraph::load(const std::string& path){}

    struct PR_Vertex_R {
        double* &scores;
        gapbs::pvector<double>& outgoing_contrib;
        SparseMatrixV<true, uint32_t>* &G;
        double &dangling_sum;
        PR_Vertex_R(double* &_scores, SparseMatrixV<true, unsigned int>*& _G, gapbs::pvector<double>& _outgoing_contrib, double& _dangling_sum) : dangling_sum(_dangling_sum), scores(_scores), G(_G), outgoing_contrib(_outgoing_contrib) {}
        inline bool operator()(uint32_t i) {
            if(G->getDegree(i) == 0){ // this is a sink
                dangling_sum = scores[i];
            } else {
                outgoing_contrib[i] = scores[i] / G->getDegree(i) == 0;
            }
            return true;
        }
    };


    struct PR_Vertex_F {
        double* &scores;
        SparseMatrixV<true, uint32_t>* &G;
        double &dangling_sum;
        const double &base_score;
        double &damping_factor;
        double* &incoming_total;
        PR_Vertex_F(double* &_scores, SparseMatrixV<true, unsigned int>*& _G,
                    double& _dangling_sum, double& _damping_factor, const double& _base_score, double* &_incoming_total) : dangling_sum(_dangling_sum), scores(_scores), G(_G),
                                                                                                                           base_score(_base_score), incoming_total(_incoming_total), damping_factor(_damping_factor) {}
        inline bool operator()(uint32_t i) {
            scores[i] = base_score +  damping_factor * (incoming_total[i] +  dangling_sum);
            return true;
        }
    };


// template <class vertex>
    /*template <typename T, typename SM> struct PR_F {
        double* &scores;
        SparseMatrixV<true, uint32_t>* &G;
        gapbs::pvector<double>& outgoing_contrib;
        double base_score;
        const double damping_factor;
        double *incoming_total;
        double dangling_sum;
        // vertex* V;
        // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
        PR_F(T *_scores, gapbs::pvector<double>& _outgoing_contrib, const SM &_G)
                : scores(_scores), outgoing_contrib(_outgoing_contrib), G(_G) {}
        inline bool update(uint32_t s, uint32_t d) {
            incoming_total[s] += outgoing_contrib[d];
            return true;
        }

        inline bool updateAtomic([[maybe_unused]] el_t s,
                                 [[maybe_unused]] el_t d) { // atomic Update
            printf("should never be called for now since its always dense\n");

            return true;
        }
        inline bool cond([[maybe_unused]] el_t d) { return true; }
    };
*/
    /*
GAP Benchmark Suite
Kernel: PageRank (PR)
Author: Scott Beamer

Will return pagerank scores for all vertices once total change < epsilon

This PR implementation uses the traditional iterative approach. This is done
to ease comparisons to other implementations (often use same algorithm), but
it is not necesarily the fastest way to implement it. It does perform the
updates in the pull direction to remove the need for atomics.
*/


    static unique_ptr<double[]> do_pagerank(SparseMatrixV<true, uint32_t>* g, uint64_t num_vertices, uint64_t num_iterations, double damping_factor, utility::TimeoutService& timer) {
        // init
        const double init_score = 1.0 / num_vertices;
        const double base_score = (1.0 - damping_factor) / num_vertices;
        unique_ptr<double[]> ptr_scores{ new double[num_vertices]() }; // avoid memory leaks
        unique_ptr<double[]> ptr_incoming{ new double[num_vertices]() }; // avoid memory leaks
        double* incoming_total = ptr_incoming.get();
        double* scores = ptr_scores.get();
#pragma omp parallel for
        for(uint64_t v = 0; v < num_vertices; v++){
            scores[v] = init_score;
            incoming_total = 0;
        }
        gapbs::pvector<double> outgoing_contrib(num_vertices, 0.0);

        VertexSubset Frontier = VertexSubset(0, num_vertices, true);

        // pagerank iterations
        for(uint64_t iteration = 0; iteration < num_iterations && !timer.is_timeout(); iteration++){
            double dangling_sum = 0.0;

            // for each node, precompute its contribution to all of its outgoing neighbours and, if it's a sink,
            // add its rank to the `dangling sum' (to be added to all nodes).
            g->vertexMap(Frontier, PR_Vertex_R(scores, g, outgoing_contrib, dangling_sum), false);
            //printf("pass 1\n");
            dangling_sum /= num_vertices;
//        COUT_DEBUG("[" << iteration << "] base score: " << base_score << ", dangling sum: " << dangling_sum);
            g->vertexMap(Frontier, PR_Vertex_F(scores, g, dangling_sum, damping_factor, base_score, incoming_total), false);

        }
        return ptr_scores;
    }

    void SSTGraph::pagerank(uint64_t num_iterations, double damping_factor, const char* dump2file){
        utility::TimeoutService tcheck ( chrono::seconds{ m_timeout } );
        // common::Timer timer; timer.start();
        size_t num_nodes = g->get_rows();
        printf("num keys = %d\n", num_nodes);
        // execute the algorithm
        unique_ptr<double[]> ptr_rank = do_pagerank(g, num_nodes, num_iterations, damping_factor, tcheck);
        printf("pass 0.0\n");
        // if(tcheck.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

        // translate the vertex ids and set the distance of the non reached vertices to max()
        double* rank = ptr_rank.get();
        map<uint32_t, double> external_ids;
        // #pragma omp parallel for
        /*for (uint32_t seg_id = 0; seg_id < G_MM->map->get_num_segments();seg_id) {
            for (uint32_t v_id = 0; v_id < G_MM->map->NUM_VERTICES_PER_SEGMENT; v_id) {
                //vertex_info_t info = segment->at(v_id); // throws invalid_vertex?
                //if (UNLIKELY(!info.is_valid())) { continue; }
                vertex_t vertex(seg_id, v_id);
                // 1.what is the real node ID, in the external domain (e.g. user id)
                int32_t key;
                G_MM->nodeid_to_nodekey(vertex, key);
                // if(key == nullptr) continue; // this means that a mapping does not exist. It should never occur as atm we don't support vertex deletions
                // 2. retrieve the distance / weight
                size_t internal_id = G_MM->map->get_absolute_vertex_id(vertex);
                double score = rank[internal_id];
                // 3. make the association vertex name - score
                external_ids.insert({key, score});
            }
        }*/
        //if(tcheck.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }
        //if(dump2file != nullptr)
        //  save0(external_ids, dump2file);

        // store the results in the given file
        // if(dump2file != nullptr){
        COUT_DEBUG("save the results to: " << "rank_values.log");
        //fstream handle("rank_values.log", ios_base::out);
        //  if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

        //   auto hashtable = external_ids.lock_table();
        //for (auto itr = external_ids.begin(); itr != external_ids.end(); itr) {
        //  handle << itr->first << " " << itr->second << "\n";
        //}
        /*for(const auto& keyvaluepair : external_ids){
            handle << keyvaluepair.first << " " << keyvaluepair.second << "\n";
        }*/


        // handle.close();
        // }
    }
/******************************************************************************
 *                                                                            *
 *  Dump                                                                      *
 *                                                                            *
 *****************************************************************************/

    void SSTGraph::dump_ostream(std::ostream& out) const {
        // struct CSRPP* graph = STINGER;

    /*    out << "[CSRPP] Vertices: " << num_vertices() << ", edges: " << num_edges() << ", directed: " << std::boolalpha << is_directed() << ", size: NOT YET IMPLEMENTED bytes" << "\n";
        for (uint32_t seg_id = 0; seg_id < G_MM->map->get_num_segments();seg_id) {
            for (uint32_t v_id = 0; v_id < G_MM->map->NUM_VERTICES_PER_SEGMENT; v_id) {
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
                    for (edge_t w_idx = 0; w_idx < out_degree; w_idx )
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
    /*                }
                }
            }
        }
        std::flush(out);*/
    }


} // namespace