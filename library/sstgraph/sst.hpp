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

#pragma once


#include "common/spinlock.hpp"
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

#include "library/interface.hpp"
#include "third-party/libcuckoo/cuckoohash_map.hh"
#include "tbb/concurrent_hash_map.h"
#include <vector>
/*
class Spa; // forward declaration
class ll_mlcsr_ro_graph; // forward declaration
namespace gfe::utility { class TimeoutService; } // forward declaration
[[maybe_unused]] static void _bm_run_llama(); // bm tool
[[maybe_unused]] static void _test_perf_run_llama(); // Performance test

*/
template <bool is_csr_ , typename value_type> class SparseMatrixV;

#define SST_HASHMAP_WITH_TBB


namespace gfe::library {

// CSR has its own mapping impl~ to translate uint64_t to int64_t and viceversa. Problem is does not support deletions.
// If the macro is defined, use the stinger mapping to perform the translations external vertex id -> internal vertex id
// Otherwise, use libcuckoo


    class SSTGraph : public virtual UpdateInterface, public virtual LoaderInterface, public virtual GraphalyticsInterface {
    public:
       // graph_csrpp* G_MM { nullptr };

    protected:
        SSTGraph(const SSTGraph&) = delete;
        SSTGraph& operator=(const SSTGraph&) = delete;
        SparseMatrixV<true, uint32_t>* g { nullptr };
        //graph_csrpp* G_MM { nullptr };
        common::SpinLock m_mutex_vtx;
        //assert(G_MM != nullptr && "CSRPP allocation");
        // if(G_MM == nullptr) ERROR("Cannot allocate CSRPP graph");
        //graph_csrpp &G_MM; // opaque object, container of the handle to the stinger graph
        const bool m_directed; // is the graph directed
        std::atomic<uint64_t> m_num_vertices = 0; // number of vertices
        uint64_t m_timeout = 0; // available time, in seconds, to complete the computation

#if defined(SST_HASHMAP_WITH_TBB)
        tbb::concurrent_hash_map<uint64_t, /* node_t */ uint64_t> m_vmap; // vertex dictionary, from external vertex ID to internal vertex ID
#endif

    public:

        /**
         * Initialise the graph instance
         * @param directed is the graph directed?
         */
        SSTGraph(bool directed);

        SSTGraph(bool directed, uint32_t num_vertices);
        /**
         * Initialise the graph instance
         * @param num_vertices number of vertices
         */
        //  CSRPP(size_t num_vertices);

        /**
         * Destructor
         */
        virtual ~SSTGraph();

        uint64_t get_internal_vertex_id(uint64_t external_vertex_id) const;

        /**
         * Get the number of edges contained in the graph
         */
        virtual uint64_t num_edges() const;

        /**
         * Get the number of nodes stored in the graph
         */
        virtual uint64_t num_vertices() const;

        /**
         * Returns true if the given vertex is present, false otherwise
         */
        virtual bool has_vertex(uint64_t vertex_id) const;

        /**
         * Returns the weight of the given edge is the edge is present, or NaN otherwise
         */
        virtual double get_weight(uint64_t source, uint64_t destination) const;

        /**
         * Check whether the graph is directed
         */
        virtual bool is_directed() const;

        /**
         * Dump the content of the graph to stdout
         */
        virtual void dump_ostream(std::ostream& out) const;

        /**
         * Impose a timeout on each graph computation. A computation that does not terminate by the given seconds will raise a TimeoutError.
         */
        virtual void set_timeout(uint64_t seconds);

        /**
         * Add the given vertex to the graph
         * @return true if the vertex has been inserted, false otherwise
         */
        virtual bool add_vertex(uint64_t vertex_id);

        virtual bool insert_batch_vertices(std::vector<uint64_t> &vertices);
        /**
         * Remove the given vertex and all edges attached to it.
         * @return true in case of success, false otherwise
         */
        virtual bool remove_vertex(uint64_t vertex_id);

        /**
         * Add the given edge in the graph
         * @return true if the edge has been inserted or updated, false in case of error
         */
        virtual bool add_edge(graph::WeightedEdge e);

        /**
         * Add the given edge in the graph. Implicitly create the referred vertices if they do not already exist
         * @return true if the edge has been inserted, false otherwise (e.g. this edge already exists)
         */
        virtual bool add_edge_v2(gfe::graph::WeightedEdge e);

        /**
         * Remove the given edge from the graph
         * @return true if the given edge has been removed, false otherwise (e.g. this edge does not exist)
         */
        virtual bool remove_edge(graph::Edge e);

        /**
         * Perform a BFS from source_vertex_id to all the other vertices in the graph.
         * @param source_vertex_id the vertex where to start the search
         * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
         */
        virtual void bfs(uint64_t source_vertex_id, const char* dump2file = nullptr);

        /**
         * Execute the PageRank algorithm for the specified number of iterations.
         *
         * @param num_iterations the number of iterations to execute the algorithm
         * @param damping_factor weight for the PageRank algorithm, it affects the score associated to the sink nodes in the graphs
         * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
         */
        virtual void pagerank(uint64_t num_iterations, double damping_factor = 0.85, const char* dump2file = nullptr);

        /**
         * Weakly connected components (WCC), associate each node to a connected component of the graph
         * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
         */
        virtual void wcc(const char* dump2file = nullptr);

        /**
         * Community Detection using Label-Propagation. Associate a label to each vertex of the graph, according to its neighbours.
         * @param max_iterations max number of iterations to perform
         * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
         */
        virtual void cdlp(uint64_t max_iterations, const char* dump2file = nullptr);

        /**
         * Local clustering coefficient. Associate to each vertex the ratio between the number of its outgoing edges and the number of
         * possible remaining edges.
         * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
         */
        virtual void lcc(const char* dump2file = nullptr);

        /**
         * Single-source shortest paths. Compute the weight related to the shortest path from the source to any other vertex in the graph.
         * @param source_vertex_id the vertex where to start the search
         * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
         */
        virtual void sssp(uint64_t source_vertex_id, const char* dump2file = nullptr);

        /**
         * Retrieve the internal handle to the library implementation
         */
        virtual void* handle();
        virtual void load(const std::string& path);

    };

/**
 * Implementation of the Graphalytics algorithms using the GAP BS source code
 * https://github.com/sbeamer/gapbs
 */
    /*  class StingerRef : public virtual Stinger {
      public:

          /**
           * Initialise the graph instance
           * @param directed is the graph directed?
           */
    //     StingerRef(bool directed);

    /**
     * Destructor
     */
    //    ~StingerRef();

    /**
     * Perform a BFS from source_vertex_id to all the other vertices in the graph.
     * @param source_vertex_id the vertex where to start the search
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    //  virtual void bfs(uint64_t source_vertex_id, const char* dump2file = nullptr);

    /**
     * Execute the PageRank algorithm for the specified number of iterations.
     *
     * @param num_iterations the number of iterations to execute the algorithm
     * @param damping_factor weight for the PageRank algorithm, it affects the score associated to the sink nodes in the graphs
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    //  virtual void pagerank(uint64_t num_iterations, double damping_factor = 0.85, const char* dump2file = nullptr);

    /**
     * Weakly connected components (WCC), associate each node to a connected component of the graph
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    //  virtual void wcc(const char* dump2file = nullptr);

    /**
     * Single-source shortest paths. Compute the weight related to the shortest path from the source to any other vertex in the graph.
     * @param source_vertex_id the vertex where to start the search
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    //  virtual void sssp(uint64_t source_vertex_id, const char* dump2file = nullptr);
    // };

} // namespace