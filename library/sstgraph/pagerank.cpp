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
    unique_ptr<double[]> do_pagerank(gfe::library::SSTGraph* g, uint64_t num_nodes, uint64_t num_iterations, double damping_factor, utility::TimeoutService& timer) {
        // init
        const double init_score = 1.0 / num_nodes;
        const double base_score = (1.0 - damping_factor) / num_nodes;
        unique_ptr<double[]> ptr_scores{ new double[num_nodes]() }; // avoid memory leaks
        double* scores = ptr_scores.get();
#pragma omp parallel for
        for(uint64_t v = 0; v < num_nodes; v){
            scores[v] = init_score;
        }
        gapbs::pvector<double> outgoing_contrib(num_nodes, 0.0);

        // pagerank iterations
        for(uint64_t iteration = 0; iteration < num_iterations && !timer.is_timeout(); iteration){
            double dangling_sum = 0.0;

            // for each node, precompute its contribution to all of its outgoing neighbours and, if it's a sink,
            // add its rank to the `dangling sum' (to be added to all nodes).
#pragma omp parallel for collapse(2) schedule(dynamic,4096) reduction(:dangling_sum)
            for (uint32_t seg_id = 0; seg_id < G_MM->map->get_num_segments();seg_id) {
                for (uint32_t v_id = 0; v_id < G_MM->map->NUM_VERTICES_PER_SEGMENT; v_id) {
                    //vertex_info_t info = segment->at(v_id); // throws invalid_vertex?
                    //if (UNLIKELY(!info.is_valid())) { continue; }
                    auto &segment = G_MM->map->get_segment(seg_id);
                    vertex_info_t info = segment.at(v_id);
                    if(info.length == -1) continue;
                    vertex_t vertex(seg_id, v_id);
                    size_t v = G_MM->map->get_absolute_vertex_id(vertex);
                    if (v >= num_nodes) {
                        continue;
                    }
                    if(info.length == 0){ // this is a sink
                        dangling_sum = scores[v];
                    } else {
                        outgoing_contrib[v] = scores[v] / info.length;
                    }

                }
            }
            //printf("pass 1\n");
            dangling_sum /= num_nodes;
//        COUT_DEBUG("[" << iteration << "] base score: " << base_score << ", dangling sum: " << dangling_sum);

            // compute the new score for each node in the graph
#pragma omp parallel for schedule(dynamic,4096) collapse(2)
            for (uint32_t seg_id = 0; seg_id < G_MM->r_map->get_num_segments(); seg_id) {
                for (uint32_t v_id = 0; v_id < G_MM->r_map->NUM_VERTICES_PER_SEGMENT; v_id) {
                    auto &segment = G_MM->r_map->get_segment(seg_id);
                    vertex_info_t info = segment.at(v_id);

                    if(info.length == -1) continue;
                    double incoming_total = 0;

                    vertex_t vertex(seg_id, v_id);
                    size_t rank = G_MM->r_map->get_absolute_vertex_id(vertex);
                    if (unlikely(rank >= num_nodes)) {
                        continue;
                    }
                    // printf("pass 2\n");
                    vertex_t* in_neighbors = info.get_neighbors();
                    int32_t &length = info.length;
                    for (edge_t w_idx = 0; w_idx < length; w_idx ) {
                        vertex_t &w =  in_neighbors[w_idx];
                        /*#ifdef DELETION
                                                if(w.deleted == 1) {
                                                    length;
                                                    continue;
                                                }
                        #endif*/		//printf("pass 3\n");
                        size_t rank_id = G_MM->map->get_absolute_vertex_id(w);
                        incoming_total = outgoing_contrib[rank_id];
                    }
                    scores[rank] = base_score  damping_factor * (incoming_total  dangling_sum);
                }
            }


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

// Implementation based on stinger_alg/src/page_rank.c
    /*   void gm_pagerank(uint64_t max_iter, double damp, const char* dump2file){
           TIMER_INIT
           double tol = 0.001;
           bool norm = false;
           double diff = 0.0;
           int32_t cnt = 0;
           double N = 0.0;
           size_t num_nodes = G_MM->map->get_num_keys();
           double* rank = new double[num_nodes]();
           uint64_t prop_id = G_MM->map->register_property(sizeof(double));
           double* G_rank_nxt = gm_rt_allocate_double(num_nodes,gm_rt_thread_id());
           combined_rank_out__degree_0_t* G_combined_rank_out__degree_ = gm_rt_allocate_struct<combined_rank_out__degree_0_t>(num_nodes, gm_rt_thread_id());

           cnt = 0 ;
           N = (double)(num_nodes) ;

   #pragma omp parallel for collapse(2) schedule(dynamic,4096)
           for (uint32_t seg_id = 0; seg_id < G_MM->map->get_num_segments();seg_id) {
               for (uint32_t v_id = 0; v_id < G_MM->map->NUM_VERTICES_PER_SEGMENT; v_id) {
                   //vertex_info_t info = segment->at(v_id); // throws invalid_vertex?
                   //if (UNLIKELY(!info.is_valid())) { continue; }
                   auto &segment = G_MM->map->get_segment(seg_id);
                   vertex_info_t info = segment.at(v_id);
                   if(info.length == -1) continue;
                   vertex_t vertex(seg_id, v_id);
                   size_t rank_id = G_MM->map->get_absolute_vertex_id(vertex);
                   if (rank_id >= num_nodes) {
                       continue;
                   }

                   assert(rank_id < num_nodes);
                   assert(info.length >= 0);
                   G_combined_rank_out__degree_[rank_id].out__degree = (int32_t) (info.length);
                   G_combined_rank_out__degree_[rank_id].G_rank = 1 / N;
               }
           }

           do
           {
               double dangling_factor = 0.0;

               diff = 0.0 ;
               dangling_factor = (double)0 ;
               if (norm)
               {
                   double __S2 = 0.0;

                   __S2 = 0.0 ;
   #pragma omp parallel
                   {
                       double __S2_prv = 0.0;

                       __S2_prv = 0.0 ;

   #pragma omp for nowait
                       for (node_t v = 0; v < num_nodes; v )
                       {
                           if ((G_combined_rank_out__degree_[v].out__degree == 0))
                           {
                               __S2_prv = __S2_prv  G_combined_rank_out__degree_[v].G_rank ;
                           }
                       }


                       GM_ATOMIC_ADD<double>(&__S2, __S2_prv);
                   }
                   dangling_factor = damp / N * __S2 ;
               }
   #pragma omp parallel
               {
                   double diff_prv = 0.0;

                   diff_prv = 0.0 ;
                   //const size_t num_nodes = num_nodes;

                   // figure out if this is well balance
               #pragma omp for nowait schedule(dynamic,4096) collapse(2)
               for (uint32_t seg_id = 0; seg_id < G_MM->map->get_num_segments(); seg_id) {
                   for (uint32_t v_id = 0; v_id < G_MM->map->NUM_VERTICES_PER_SEGMENT; v_id) {
                       auto &segment = G_MM->r_map->get_segment(seg_id);
                       vertex_info_t info = segment.at(v_id);

                       if(info.length == -1) continue;
                       double val = 0.0;
                       double __S3 = 0.0;

                       __S3 = 0.0 ;
                       vertex_t vertex(seg_id, v_id);
                       size_t rank = G_MM->map->get_absolute_vertex_id(vertex);
                       if (unlikely(rank >= num_nodes)) {
                           continue;
                       }
                       vertex_t* neighbors = info.get_neighbors();
                       int32_t &length = info.length;
                       for (edge_t w_idx = 0; w_idx < length; w_idx )
                       {
                           vertex_t &w =  neighbors[w_idx];
   #ifdef DELETION
                           if(w.deleted == 1) {
                               length;
                               continue;
                           }
   #endif
                           size_t rank_id = G_MM->map->get_absolute_vertex_id(w);
                           __S3 = __S3  G_combined_rank_out__degree_[rank_id].G_rank / ((double)G_combined_rank_out__degree_[rank_id].out__degree) ;
                       }

                       val = (1 - damp) / N  damp * __S3  dangling_factor ;
                       diff_prv = diff_prv   std::abs((val - G_combined_rank_out__degree_[rank].G_rank))  ;
                       G_rank_nxt[rank] = val ;
                   }
               }

                   GM_ATOMIC_ADD<double>(&diff, diff_prv);
               }

   #pragma omp parallel for
               for (node_t i5 = 0; i5 < num_nodes; i5 )
               {
                   G_combined_rank_out__degree_[i5].G_rank = G_rank_nxt[i5] ;
               }
               cnt = cnt  1 ;
           }
           while ((diff > tol) && (cnt < max_iter));

   #pragma omp parallel for collapse(2)
       for (uint32_t seg_id = 0; seg_id < G_MM->map->get_num_segments(); seg_id) {
           for (uint32_t v_id = 0; v_id < G_MM->map->NUM_VERTICES_PER_SEGMENT; v_id) {

               auto &segment = G_MM->map->get_segment(seg_id);
               vertex_info_t info = segment.at(v_id);
               if(info.length == -1) continue;

               vertex_t vertex(seg_id, v_id);
               size_t rank = G_MM->map->get_absolute_vertex_id(vertex);
               if (unlikely(rank >= num_nodes)) {
                   continue;
               }

               G_MM->map->set_vertex_property<double>(vertex, prop_id, G_combined_rank_out__degree_[rank].G_rank);
           }
       }
       gm_rt_cleanup();
           //
           // store the final results (if required)
        /*   cuckoohash_map<uint64_t, double> result;
           to_external_ids(rank, num_registered_vertices, &result); // convert the internal logical IDs into the external IDs
           save(result, dump2file);*/


    // auto list_entries = result.lock_table();
    /*   for (uint32_t seg_id = 0; seg_id < G_MM->map->get_num_segments();seg_id) {
           for (uint32_t v_id = 0; v_id < G_MM->map->NUM_VERTICES_PER_SEGMENT; v_id) {
               vertex_t vertex(seg_id, v_id);
               size_t rank = G_MM->map->get_absolute_vertex_id(vertex);
             //  std::cout  << G_MM->_numeric32_reverse_key[rank] << " " << G_MM->map->get_vertex_property<double>(vertex, prop_id) << std::endl;
           }
       }

       };*/

} // namespace