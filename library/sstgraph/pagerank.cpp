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
#include "sst.hpp"

#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <fstream>
#include <limits>
#include <mutex>
#include <queue>
#include <thread>
#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>
#include <cmath>
#include <omp.h>



#include "common/system.hpp"
#include "common/timer.hpp"
#include "third-party/gapbs/gapbs.hpp"
#include "third-party/libcuckoo/cuckoohash_map.hh"
#include "utility/timeout_service.hpp"

using namespace common;
using namespace gfe::utility;
using namespace libcuckoo;
using namespace std;

// Macros
#define TIMER_INIT auto time_start = chrono::steady_clock::now();
#define CHECK_TIMEOUT if(m_timeout > 0 && (chrono::steady_clock::now() - time_start) > chrono::seconds(m_timeout)) { \
        RAISE_EXCEPTION(TimeoutError, "Timeout occurred after: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - time_start).count() << " seconds") };

/******************************************************************************
 *                                                                            *
 *  Debug                                                                     *
 *                                                                            *
 *****************************************************************************/
//#define DEBUG
extern mutex _log_mutex [[maybe_unused]];
#define COUT_DEBUG_FORCE(msg) { scoped_lock<mutex> lock(::_log_mutex); std::cout << "[Stinger::" << __FUNCTION__ << "] [" << common::concurrency::get_thread_id() << "] " << msg << std::endl; }
#if defined(DEBUG)
#define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
#define COUT_DEBUG(msg)
#endif


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

template <typename T, typename SM> struct PR_Vertex {
    T *scores;
    T *outgoing_contrib;
    const SM &G;
    double dangling_sum;
    PR_Vertex(T *_scores, const SM &_G, T *_outgoing_contrib, double _dangling_sum) : dangling_sum(_dangling_sum), scores(_scores), G(_G), outgoing_contrib(_outgoing_contrib) {}
    inline bool operator()(uint32_t i) {
        if(G.getDegree(i) == 0){ // this is a sink
            dangling_sum = scores[i];
        } else {
            outgoing_contrib[i] = scores[i] / G.getDegree(i) == 0;
        }
        return true;
    }
};


template <typename T, typename SM> struct PR_Vertex_F {
    T *scores;
    const SM &G;
    double dangling_sum;
    const double base_score;
    const double damping_factor;
    const double *incoming_total;
    PR_Vertex_F(T *_scores, const SM &_G,
                double _dangling_sum, double _damping_factor, double _base_score, double *_incoming_total) : dangling_sum(_dangling_sum), scores(_scores), G(_G),
                base_score(_base_score), incoming_total(_incoming_total), damping_factor(_damping_factor) {}
    inline bool operator()(uint32_t i) {
        scores[i] = base_score +  damping_factor * (incoming_total[i] +  dangling_sum);
        return true;
    }
};


// template <class vertex>
template <typename T, typename SM> struct PR_F {
    T *scores, *outgoing_contrib;
    const double base_score;
    const double damping_factor;
    double *incoming_total;
    double dangling_sum;
    // vertex* V;
    // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
    const SM &G;
    PR_F(T *_scores, T *_outgoing_contrib, const SM &_G)
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
}; // from ligra readme: for cond which always ret true, ret cond_true// return
// cond_true(d); }};
/******************************************************************************
 *                                                                            *
 *  Pagerank                                                                  *
 *                                                                            *
 ******************************************************************************/
namespace gfe::library {


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



    static
    unique_ptr<double[]> do_pagerank(gfe::library::SSTGraph* g, uint64_t num_vertices, uint64_t num_iterations, double damping_factor, utility::TimeoutService& timer) {
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
        for(uint64_t iteration = 0; iteration < num_iterations && !timer.is_timeout(); iteration){
            double dangling_sum = 0.0;

            // for each node, precompute its contribution to all of its outgoing neighbours and, if it's a sink,
            // add its rank to the `dangling sum' (to be added to all nodes).
            g->vertexMap(Frontier, PR_Vertex(scores, outgoing_contrib, g, dangling_sum), false);
            //printf("pass 1\n");
            dangling_sum /= num_vertices;
//        COUT_DEBUG("[" << iteration << "] base score: " << base_score << ", dangling sum: " << dangling_sum);
            auto frontier2 = g->edgeMap(Frontier, PR_F(p_curr, p_next, G), false, 20);
            g->vertexMap(Frontier, PR_Vertex_F(scores, g, dangling_sum, damping_factor, base_score, incoming_total), false);

        }
        return ptr_scores;
    }

    void SSTGraph::pagerank(uint64_t num_iterations, double damping_factor, const char* dump2file){
        utility::TimeoutService tcheck ( chrono::seconds{ m_timeout } );
        common::Timer timer; timer.start();
        size_t num_nodes = g->get_rows();
        printf("num keys = %d\n", num_nodes);
        // execute the algorithm
        unique_ptr<double[]> ptr_rank = do_pagerank(g, num_nodes, num_iterations, damping_factor, tcheck);
        printf("pass 0.0\n");
        if(tcheck.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

        // translate the vertex ids and set the distance of the non reached vertices to max()
        double* rank = ptr_rank.get();
        map<uint32_t, double> external_ids;
        // #pragma omp parallel for
        for (uint32_t seg_id = 0; seg_id < G_MM->map->get_num_segments();seg_id) {
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
        }
        if(tcheck.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }
        //if(dump2file != nullptr)
        //  save0(external_ids, dump2file);

        // store the results in the given file
        // if(dump2file != nullptr){
        COUT_DEBUG("save the results to: " << "rank_values.log");
        fstream handle("rank_values.log", ios_base::out);
        if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

        //   auto hashtable = external_ids.lock_table();
        for (auto itr = external_ids.begin(); itr != external_ids.end(); itr) {
            handle << itr->first << " " << itr->second << "\n";
        }
        /*for(const auto& keyvaluepair : external_ids){
            handle << keyvaluepair.first << " " << keyvaluepair.second << "\n";
        }*/

        handle.close();
        // }
    }


} // namespace