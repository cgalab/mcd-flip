#include "geom.h"

#include <numeric>
#include <cmath>
#include <iomanip>

#ifndef NDEBUG
/** Local sanity check. */
void
Edge::
assert_valid() const {
  assert(opposite);
  assert(this == opposite->opposite);

  assert(this == prev->next);
  assert(this == next->prev);
  assert(v == opposite->prev->v);
}
#endif



/** Comperator function for sorting edges of a vertex.
 *
 * Sorts counterclockwise, starting in the direction of the positive
 * x-axis. */
int
Edge::
edge_direction_comperator(const Edge* ea, const Edge* eb, const Vertex* v) {
  const Vertex* const va = ea->v;
  const Vertex* const vb = eb->v;

  assert(va != vb);
  assert(va != v);
  assert(vb != v);
  int sign;

  /* sign of determinant does not provide a total order,
   * since it's quite possible for a<b, b<c and c<a.
   *
   * So we also split by quadrants.
   *
   * quadrants:
   *  (is_right << 1) + is_top          quadrant #
   *           1 3                         1 0
   *           0 2                         2 3
   */
  static const int quadrant_map[4] = {2, 1, 3, 0};
  unsigned quadrant_a, quadrant_b, is_right, is_top;

  is_right = (va->x >= v->x);
  is_top   = (va->y >= v->y);
  quadrant_a = quadrant_map[(is_right << 1) + is_top];

  is_right = (vb->x >= v->x);
  is_top   = (vb->y >= v->y);
  quadrant_b = quadrant_map[(is_right << 1) + is_top];

  if (quadrant_a == quadrant_b) {
    return Vertex::orientation(*v, *vb, *va);
  } else {
    return quadrant_a - quadrant_b;
  }
}


DCEL::
DCEL(VertexList&& vertices,
     const InputEdgeSet& edges)
  : all_vertices(std::move(vertices))
{
  all_edges.reserve(edges.size()*2);
  std::vector< std::vector<Edge*> > vertex_edgelist;  // Edges pointing away from each vertex.

  vertex_edgelist.resize(all_vertices.size());

  for (const auto& e : edges) {
    Edge* e1 = &all_edges.back() + 1;
    Edge* e2 = e1 + 1;

    all_edges.emplace_back(Edge(e2, &all_vertices[e.first], all_edges.size()));
    all_edges.emplace_back(Edge(e1, &all_vertices[e.second], all_edges.size()));

    vertex_edgelist[e.first ].push_back(e2);
    vertex_edgelist[e.second].push_back(e1);
  }

  for (unsigned i=0; i < all_vertices.size(); ++i) {
    std::vector<Edge*>& vel = vertex_edgelist[i];
    Vertex* v = &all_vertices[i];

    std::sort(vel.begin(), vel.end(), [v](const Edge* a, const Edge *b){ return Edge::edge_direction_comperator(a, b, v); });

    Edge* prev = vel.back();
    for (auto outgoing_e_ptr : vel) {
      prev->prev = outgoing_e_ptr->opposite;
      outgoing_e_ptr->opposite->next = prev;

      prev = outgoing_e_ptr;
    }

    v->degree = vel.size();
    v->incident_edge = vel.front()->opposite;
    v->update_vertex_is_of_higher_degree();
  }

  assert_valid();

  //num_faces = ...
}

void
DCEL::
improve_convex_decomposition() {
}

#if 0
void
DCEL::
find_convex_decomposition_many(unsigned num_iterations) {
  DBG_FUNC_BEGIN(DBG_DECOMPOSITION_LOOP);

  unsigned initial_num_faces_to_beat = num_faces;
  unsigned num_faces_to_beat = initial_num_faces_to_beat;
  bool have_solution = false;
  int solution_from_iter = -1;
  bool current_is_best = false;
  SavedDecomposition best = SavedDecomposition(working_set, num_faces);
  for (unsigned iter = 0; iter < num_iterations; ++iter) {
    DBG(DBG_DECOMPOSITION_LOOP) << "Resetting constraints";

    assert_valid();
    flip_random_edges_and_reset_constraints();
    unconstrain_random_edges();

    current_is_best = (num_faces < num_faces_to_beat);
    if (current_is_best) {
      DBG(DBG_DECOMPOSITION_LOOP) << "Iteration " << iter << "/" << num_iterations << ": This solution: " << num_faces << "; NEW BEST; previous best: " << num_faces_to_beat;
      have_solution = true;
      solution_from_iter = iter;
      num_faces_to_beat = num_faces;
      best = SavedDecomposition(working_set, num_faces);
    } else {
      DBG(DBG_DECOMPOSITION_LOOP) << "Iteration " << iter << "/" << num_iterations << ": This solution: " << num_faces << "; current best: " << num_faces_to_beat;
    }
  }

  if (have_solution) {
    DBG(DBG_GENERIC | DBG_DECOMPOSITION_LOOP) << "Done " << num_iterations << " iterations.  Best now is  " << num_faces_to_beat << " from " << initial_num_faces_to_beat
      << " (" << std::setw(3) << (initial_num_faces_to_beat-num_faces_to_beat) << ")"
      << "; found in iteration (#/max/#triangles): "
      << solution_from_iter
      << ", " << num_iterations
      << ", " << working_set.shuffled_edges.size()
      ;
  } else {
    DBG(DBG_GENERIC | DBG_DECOMPOSITION_LOOP) << "Done " << num_iterations << " iterations.  No new best: " << num_faces_to_beat << " from " << initial_num_faces_to_beat
      << " (  -)"
      << "; found in iteration (#/max/#triangles): "
      << solution_from_iter
      << ", " << num_iterations
      << ", " << working_set.shuffled_edges.size()
      ;
  };
  /** We keep the current decomposition and triangulation if it is at least as good as the one before.  We do not require it be better. */
  if (num_faces > num_faces_to_beat) {
    DBG(DBG_DECOMPOSITION_LOOP) << "Re-injecting saved decomposition with " << best.saved_num_faces << " faces because we have " << num_faces << " right now.";
    reinject_saved_decomposition(std::move(best));
  } else {
    DBG(DBG_DECOMPOSITION_LOOP) << "Keeping currently best known decomposition with " << num_faces << " faces.";
    #if 0
    LOG(INFO) << "Keeping currently best known decomposition with " << num_faces << " faces."; /* We like this during poor man's timing tests. */
    #endif
  }
  assert_valid();
  assert_hole_shooting_reset();
  assert(num_faces_to_beat == num_faces);

  DBG_FUNC_END(DBG_DECOMPOSITION_LOOP);
}

/** Finds an initial, or improves an existing convex decomposition.
 *
 * First we unconstrain all possible in a random way, then we try to improve locally.
 *
 * To restart the process, run reset_constraints().
 */
void
DCEL::
find_convex_decomposition() {
  DBG_FUNC_BEGIN(DBG_DECOMPOSITION_LOOP);

  assert_valid();
  if (!initial_constrained && working_set.num_my_triangles == num_faces) {
    flip_random_edges_and_reset_constraints();
    unconstrain_random_edges();
  };
  assert_valid();
  shoot_holes();

  assert_valid();
  assert_hole_shooting_reset();

  DBG_FUNC_END(DBG_DECOMPOSITION_LOOP);
}
#endif


/** Runs validity checks for all the edges */
#ifndef NDEBUG
void
DCEL::
assert_valid() const {
  for (const auto &v : all_vertices) {
    v.assert_valid();
  }
  for (const auto &e : all_edges) {
    if (!e.is_alive) continue;
    e.assert_valid();
  }
}
#endif

/** Write the (constraint) segments to the output stream in obj format
 *
 * If we have a vertex list, also print those.
 */
void
DCEL::
write_obj_segments(bool dump_vertices, std::ostream &o) {
  if (dump_vertices) {
    for (const auto &v : all_vertices) {
      o << "v " << std::setprecision(15) << v.x << " " << v.y << " 0" << std::endl;
    }
  }

  for (const auto &e : all_edges) {
    if (!e.is_alive) continue;
    if (e.opposite && e.opposite < &e) continue;

    int tail_idx = (e.get_tail() - all_vertices.data())+1;
    int head_idx = (e.v - all_vertices.data())+1;

    o << "l "
      << tail_idx
      << " "
      << head_idx
      << std::endl;
  }
}

std::ostream& operator<<(std::ostream& os, const Vertex& v) {
  os << "v"
     << v.idx()
     << "(" << v.x << "; " << v.y << ")";
  return os;
}
std::ostream& operator<<(std::ostream& os, const Edge& e) {
  os << "e"
     << e.idx()
     << "(" << *e.get_tail()
     << "->" << *e.v;
  os << "; n/p/o";
  os << ": ";

  os << (e.next->idx())
     << "/" << (e.prev->idx())
     << "/" << (e.opposite ? e.opposite->idx() : -1);
  os << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os, const DCEL& d) {
  os << "DCEL" << std::endl;
  os << " #edges: " << d.all_edges.size() << std::endl;
  for (unsigned i=0; i<d.all_edges.size(); ++i) {
    const Edge &e = d.all_edges[i];

    os << "   edge #" << e
       << std::endl;
  }
  return os;
}
/* vim: set fdm=marker: */
