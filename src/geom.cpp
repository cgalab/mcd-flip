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
  assert(prev);
  assert(next);

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

  assert(ea->opposite->v == v);
  assert(eb->opposite->v == v);

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
  int quadrant_a, quadrant_b;
  unsigned is_right, is_top;

  is_right = (va->x >= v->x);
  is_top   = (va->y >= v->y);
  quadrant_a = quadrant_map[(is_right << 1) + is_top];

  is_right = (vb->x >= v->x);
  is_top   = (vb->y >= v->y);
  quadrant_b = quadrant_map[(is_right << 1) + is_top];

  #if 0
  if (v->idx_ == 10) {
    DBG(DBG_SETUP) << "from " << *v;
    DBG(DBG_SETUP) << " ea goes to " << *va;
    DBG(DBG_SETUP) << "  qa " << quadrant_a;
    DBG(DBG_SETUP) << " eb goes to " << *vb;
    DBG(DBG_SETUP) << "  qb " << quadrant_b;
    if (quadrant_a == quadrant_b) {
      DBG(DBG_SETUP) << " orient:" << Vertex::orientation(*v, *vb, *va);
    } else {
      DBG(DBG_SETUP) << " diff  :" << (quadrant_a - quadrant_b);
    }
  }
  #endif

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
  DBG_FUNC_BEGIN(DBG_SETUP);

  num_faces = 1 - all_vertices.size() + edges.size();
  DBG(DBG_SETUP) << "# faces: " << num_faces;

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

    std::sort(vel.begin(), vel.end(), [v](const Edge* a, const Edge *b){ return Edge::edge_direction_comperator(a, b, v) < 0; });
    // DBG(DBG_SETUP) << "at " << *v;

    Edge* prev = vel.back();
    for (auto outgoing_e_ptr : vel) {
      // DBG(DBG_SETUP) << " e to " << *outgoing_e_ptr->v;
      prev->prev = outgoing_e_ptr->opposite;
      outgoing_e_ptr->opposite->next = prev;

      prev = outgoing_e_ptr;
    }
    v->degree = vel.size();
    v->incident_edge = vel.front()->opposite;
  }

  /* Find CH vertices */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
  auto it_v_min = std::min_element(all_vertices.begin(), all_vertices.end(), [](const Vertex& a, const Vertex& b){ return (a.x < b.x) || ((a.x == b.x) && (a.y < b.y)); });
#pragma GCC diagnostic pop
  Vertex *v_min = &*it_v_min;
  DBG(DBG_SETUP) << "v_min:" << *v_min;
  Edge *ch_edge = NULL;
  /* Walk around v in clockwise order.  Once the det(v, prev, cur) is postive, cur is a CH vertex. */

  const Vertex *v_prev = v_min->incident_edge->opposite->prev->opposite->v;
  DCEL::AroundVertexFacesIterator it(v_min->incident_edge);
  for (; *it; ++it) {
    const Vertex *v_cur = it->opposite->v;
    if (Vertex::orientation(*v_min, *v_prev, *v_cur) > 0) {
      ch_edge = it->opposite;
      break;
    }
    v_prev = v_cur;
  };
  assert(ch_edge);
  unsigned num_ch = 0;
  Edge *e = ch_edge;
  do {
    e->is_on_ch = true;
    DBG(DBG_SETUP) << "Noting that " << *e << " is on the CH";
    assert( Vertex::orientation(*e->opposite->v, *e->v, *e->next->v) <= 0);
    e = e->next;
    ++num_ch;
  } while (e != ch_edge);
  DBG(DBG_SETUP) << "# ch vertices: " << num_ch;

  DBG(DBG_SETUP) << "# vertices: " << all_vertices.size();
  DBG(DBG_SETUP) << "# edges: " << edges.size();

  for (auto& v : all_vertices) {
    v.update_vertex_is_of_higher_degree();
    if (v.is_of_higher_degree) {
      DBG(DBG_SETUP) << " of higher degree: " << v;
      v.idx_in_higher_degree_vertices = higher_degree_vertices.size();
      higher_degree_vertices.push_back(&v);
    }
  }

  assert_valid();

  DBG_FUNC_END(DBG_SETUP);
}

void
DCEL::
improve_convex_decomposition() {
  DBG_FUNC_BEGIN(DBG_IMPROVE);

  do {
    /* Pick a vertex */
    Vertex* v;
    int v_idx_in_higher_degree_vertices;
    {
      std::uniform_int_distribution<unsigned> vertex_picker(0, higher_degree_vertices.size() - 1);
      v_idx_in_higher_degree_vertices = vertex_picker(random_engine);

      v = higher_degree_vertices[v_idx_in_higher_degree_vertices];
      assert(v->idx_in_higher_degree_vertices == v_idx_in_higher_degree_vertices);
      DBG(DBG_IMPROVE) << "picked:" << *v;
    }

    /* Pick an edge */
    Edge* e;
    {
      std::uniform_int_distribution<unsigned> edge_picker(0, v->degree - 1);
      unsigned e_idx = edge_picker(random_engine);

      DCEL::AroundVertexFacesIterator it(v->incident_edge);
      for (unsigned i=e_idx; i>0; ++it, --i) {};
      assert(*it);

      e = *it;
      DBG(DBG_IMPROVE) << "first picked:" << *e;
      DCEL::AroundVertexFacesIterator it2(e);
      for (; *it2 && !(it2->can_remove_at_tip()); ++it2) {};
      assert(*it2);
      assert(it2->can_remove_at_tip());
      e = *it2;
      assert(! e->is_on_ch && !e->opposite->is_on_ch);
      DBG(DBG_IMPROVE) << "picked:" << *e;
    };

    /* Move left or right */
    std::uniform_int_distribution<unsigned> direction_picker(0, 1);
    bool move_right = direction_picker(random_engine);

    /* Check if legal to move.  It must not coincide with any edges once moved,
     * and the angle at the tail end must not get > pi.
     *
     * It coincides either when there only is a triangle on the side, or when
     * wheter are only pi-degree vertices between the tail and the new tip.
     */
    if (move_right) {
      if (e->v == e->opposite->next->next->v) {
        DBG(DBG_IMPROVE) << " move would collapse a triangle";
        break;
      }
      Vertex* tail = e->opposite->v;
      const Vertex* new_tip = e->opposite->prev->opposite->v;
      const Vertex* in_between = e->opposite->next->v;
      assert(new_tip != in_between);
      DBG(DBG_IMPROVE) << "tail       :" << *tail;
      DBG(DBG_IMPROVE) << "in_between :" << *in_between;
      DBG(DBG_IMPROVE) << "new_tip    :" << *new_tip;
      int o = Vertex::orientation(*tail, *in_between, *new_tip);
      DBG(DBG_IMPROVE) << "o:" << o;
      if (o == 0) {
        DBG(DBG_IMPROVE) << " moving onto collinear edges";
        break;
      }
      assert(o > 0);

      const Vertex* other_behind_tail = e->prev->opposite->v;
      if (Vertex::orientation(*other_behind_tail, *tail, *new_tip) < 0) {
        DBG(DBG_IMPROVE) << " would create reflex tail";
        break;
      }

      /* We can move the tip */
      Vertex* old_v = e->v;
      Edge* old_next = e->next;
      Edge* moved_over = e->opposite->prev;
      Edge* new_opposite_prev = moved_over->prev;

      old_v->incident_edge = moved_over;
      assert(old_v->incident_edge->v == old_v);

      moved_over->set_next(old_next);
      e->set_next(moved_over);
      new_opposite_prev->set_next(e->opposite);

      bool new_v_old_is_of_higher_degree = new_opposite_prev->v->is_of_higher_degree;
      e->v = new_opposite_prev->v;
      assert(e->v == moved_over->opposite->v);

      old_v->dec_degree();
      e->v->inc_degree();

      bool tail_old_is_of_higher_degree = tail->is_of_higher_degree;
      tail->update_vertex_is_of_higher_degree();

      assert(e->v->is_of_higher_degree);
      assert(old_v == higher_degree_vertices[v_idx_in_higher_degree_vertices]);

      if (! old_v->is_of_higher_degree) {
        DBG(DBG_IMPROVE) << " Old_v is no longer a higher degree vertex; Removing old vertex from higher_degree_vertex";
        higher_degree_vertices_remove(v_idx_in_higher_degree_vertices);
      }
      if (new_v_old_is_of_higher_degree) {
        DBG(DBG_IMPROVE) << "new_v previously was a higher degree vertex";
      } else {
        DBG(DBG_IMPROVE) << "Adding new higher_degree_vertex";
        higher_degree_vertices_append(e->v);
      }
      if (tail_old_is_of_higher_degree) {
        DBG(DBG_IMPROVE) << "Tail previously was of higher degree with degree " << tail->degree;
        if (!tail->is_of_higher_degree) {
          DBG(DBG_IMPROVE) << "But is no more!";
          higher_degree_vertices_remove(tail);
        }
      } else if (tail->is_of_higher_degree) {
        DBG(DBG_IMPROVE) << "Adding tail to higher_degree_vertex";
        higher_degree_vertices_append(tail);
      };

      assert_valid();

      if (moved_over->can_remove()) {
        DBG(DBG_IMPROVE) << "(1) We can remove an edge!";
        remove_edge(moved_over);
        assert_valid();
      };
      if (e->opposite->prev->can_remove()) {
        DBG(DBG_IMPROVE) << "(2) We can remove an edge!";
        remove_edge(e->opposite->prev);
        assert_valid();
      };
      if (e->opposite->next->can_remove()) {
        DBG(DBG_IMPROVE) << "(3) We can remove an edge!";
        remove_edge(e->opposite->next);
        assert_valid();
      };
    } else {
      if (e == e->next->next->next) {
        DBG(DBG_IMPROVE) << " move would collapse a triangle";
        break;
      }
      Vertex* tail = e->opposite->v;
      const Vertex* new_tip = e->next->v;
      const Vertex* in_between = e->prev->opposite->v;
      assert(new_tip != in_between);
      DBG(DBG_IMPROVE) << "tail       :" << *tail;
      DBG(DBG_IMPROVE) << "in_between :" << *in_between;
      DBG(DBG_IMPROVE) << "new_tip    :" << *new_tip;
      int o = Vertex::orientation(*tail, *in_between, *new_tip);
      DBG(DBG_IMPROVE) << "o:" << o;
      if (o == 0) {
        DBG(DBG_IMPROVE) << " moving onto collinear edges";
        break;
      }
      assert(o < 0);

      const Vertex* other_behind_tail = e->opposite->next->v;
      if (Vertex::orientation(*other_behind_tail, *tail, *new_tip) > 0) {
        DBG(DBG_IMPROVE) << " would create reflex tail";
        break;
      }

      /* We can move the tip */
      Vertex* old_v = e->v;
      Edge* old_op_prev = e->opposite->prev;
      Edge* moved_over = e->next;
      Edge* new_next = moved_over->next;

      old_v->incident_edge = old_op_prev;
      assert(old_v->incident_edge->v == old_v);

      moved_over->set_next(e->opposite);
      old_op_prev->set_next(moved_over);
      e->set_next(new_next);

      bool new_v_old_is_of_higher_degree = moved_over->v->is_of_higher_degree;
      e->v = moved_over->v;
      assert(e->v == new_next->opposite->v);

      old_v->dec_degree();
      e->v->inc_degree();

      bool tail_old_is_of_higher_degree = tail->is_of_higher_degree;
      tail->update_vertex_is_of_higher_degree();

      assert(e->v->is_of_higher_degree);
      assert(old_v == higher_degree_vertices[v_idx_in_higher_degree_vertices]);

      if (! old_v->is_of_higher_degree) {
        DBG(DBG_IMPROVE) << " Old_v is no longer a higher degree vertex; Removing old vertex from higher_degree_vertex";
        higher_degree_vertices_remove(v_idx_in_higher_degree_vertices);
      }
      if (new_v_old_is_of_higher_degree) {
        DBG(DBG_IMPROVE) << "new_v previously was a higher degree vertex";
      } else {
        DBG(DBG_IMPROVE) << "Adding new higher_degree_vertex";
        higher_degree_vertices_append(e->v);
      }
      if (tail_old_is_of_higher_degree) {
        DBG(DBG_IMPROVE) << "Tail previously was of higher degree with degree " << tail->degree;
        if (!tail->is_of_higher_degree) {
          DBG(DBG_IMPROVE) << "But is no more!";
          higher_degree_vertices_remove(tail);
        }
      } else if (tail->is_of_higher_degree) {
        DBG(DBG_IMPROVE) << "Adding tail to higher_degree_vertex";
        higher_degree_vertices_append(tail);
      };

      assert_valid();

      if (moved_over->can_remove()) {
        DBG(DBG_IMPROVE) << "(1) We can remove an edge!";
        remove_edge(moved_over);
        assert_valid();
      };
      if (new_next->can_remove()) {
        DBG(DBG_IMPROVE) << "(2) We can remove an edge!";
        remove_edge(new_next);
        assert_valid();
      };
      if (e->prev->can_remove()) {
        DBG(DBG_IMPROVE) << "(3) We can remove an edge!";
        remove_edge(e->prev);
        assert_valid();
      };
    }
  } while (0);

  DBG_FUNC_END(DBG_IMPROVE);
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
  unsigned cnt_high_degree = 0;
  for (const auto &v : all_vertices) {
    v.assert_valid();
    if (v.is_of_higher_degree) {
      ++cnt_high_degree;
      assert(v.idx_in_higher_degree_vertices >= 0);
      assert(higher_degree_vertices[ v.idx_in_higher_degree_vertices ] == &v);
    } else {
      assert(v.idx_in_higher_degree_vertices == -1);
    }
  }
  for (const auto &e : all_edges) {
    if (!e.is_alive) continue;
    e.assert_valid();
  }

  assert(cnt_high_degree == higher_degree_vertices.size());
  for (unsigned i=0; i<higher_degree_vertices.size(); ++i) {
    assert(higher_degree_vertices[i]->idx_in_higher_degree_vertices >= 0);
    assert((unsigned)higher_degree_vertices[i]->idx_in_higher_degree_vertices == i);
    assert(higher_degree_vertices[i]->is_of_higher_degree);
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
