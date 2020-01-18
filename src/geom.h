#pragma once

#include "mcd.h"
#include "tools.h"

#include <vector>
#include <unordered_set>

class Edge;
class DCEL;

// {{{ Vertex
class Vertex {
public:
  const double x;
  const double y;

private:
  unsigned degree = 0;
  bool is_of_higher_degree = false;
  int idx_in_higher_degree_vertices = -1;
  Edge* incident_edge = NULL; /* Edge pointing towards v */

  bool get_vertex_is_of_higher_degree() const;
  void update_vertex_is_of_higher_degree() {
    is_of_higher_degree = get_vertex_is_of_higher_degree();
  }

  void dec_degree() {
    --degree;
    update_vertex_is_of_higher_degree();
  }
  void inc_degree() {
    ++degree;
    assert(get_vertex_is_of_higher_degree());
    is_of_higher_degree = true;
  }

#ifndef NDEBUG
  unsigned find_degree() const;
#endif

public:
#ifndef NDEBUG
  const int idx_;
  Vertex(double x_, double y_, int idx)
    : x(x_)
    , y(y_)
    , idx_(idx)
  {}
  int idx() const { return idx_; };
  void assert_valid() const;
#else
  Vertex(double x_, double y_, int)
    : x(x_)
    , y(y_)
  {}
  int idx() const { return 0; };
  void assert_valid() const {};
#endif

  /** determine Vertices' a, b, c's relative orientation.
   *
   * returns 1 if they are counterclockwise, 0 if collinear, -1 if clockwise.
   */
  static int orientation(const Vertex& a, const Vertex &b, const Vertex &c) {
    double det1 = (a.x - c.x) * (b.y - c.y);
    double det2 = (a.y - c.y) * (b.x - c.x);
    double det = det1 - det2;
    return signum(det);
  }

  friend class Edge;
  friend class DCEL;
  friend std::ostream& operator<<(std::ostream&, const Vertex&);
};
using VertexList = std::vector<Vertex>;
// }}}

struct PairHasher {
    std::size_t operator()(const std::pair<unsigned, unsigned> &p) const { return std::hash<unsigned>()(p.first) ^ std::hash<unsigned>()(p.second); }
};
using InputEdgeSet = std::unordered_set<std::pair<unsigned,unsigned>, PairHasher>;

// {{{ Edge
/** A multi-layer half-edge data structure.
 *
 * Edges come in two forms:  constraints and triangulation edges.
 *
 * The triangulation edges themselves make up a classic DCEL,
 * with buddy-pointer, next and prev pointer.
 *
 * Some of the triangulation edges are also constraints.  These will,
 * additionally have next/prev links to the next constrained edge.
 */
class Edge {
  private:
  bool is_alive = true;
  bool is_on_ch = false;   /** Whether this is the outside buddy of an edge on the CH.  This means the left side is outside and the angle spanned to the next is >= pi */
  Edge *opposite;          /** Pointer to the buddy of this edge. */
  Edge *next = NULL;       /** Pointer to the next edge of this triangle. This edge will start at Vertex v. */
  Edge *prev = NULL;       /** Pointer to the previous face of this triangle. */
  Vertex *v;               /** Vertex this edge points to.  Tail of prev. */

#ifndef NDEBUG
  const int idx_;
#endif

  Edge(Edge *opposite_, Vertex *v_, [[maybe_unused]] int idx)
    : opposite(opposite_)
    , v(v_)
#ifndef NDEBUG
    , idx_(idx)
#endif
  {}
  friend class DCEL;
  friend class Vertex;
  friend std::ostream& operator<<(std::ostream&, const DCEL&);

#ifndef NDEBUG
  int idx() const { return idx_; };
#else
  int idx() const { return 0; };
#endif

public:
  void print_tip_and_neighbors() const {
    assert(opposite);

    const Vertex &tip = *v;
    const Vertex &left = *next->v;
    const Vertex &right = *opposite->prev->prev->v;

    DBG(DBG_GENERIC) << "right (" << right.x << ", " << right.y << ")";
    DBG(DBG_GENERIC) << "tip   (" << tip.x << ", " << tip.y << "); tail (" << get_tail()->x << ", " << get_tail()->y << ")";
    DBG(DBG_GENERIC) << "left  (" << left.x << ", " << left.y << ")";
    DBG(DBG_GENERIC) << "det (r,t,l)" << Vertex::orientation(right, tip, left);
  }

protected:
  bool can_remove_at_tip() const {
    assert(v->is_of_higher_degree); /* otherwise we shouldn't call you in the first place */
    assert(opposite);

    const Vertex &tip = *v;
    const Vertex &left = *next->v;
    const Vertex &right = *opposite->prev->prev->v;

    return ((&left != &right) && Vertex::orientation(right, tip, left) >= 0);
  }
private:
  static int edge_direction_comperator(const Edge* ea, const Edge* eb, const Vertex* v);
public:
  /** Check whether this edge can be removed
   *
   * Edges can be removed if they are not on the CH,
   * and if removing them will not result in a reflex vertex
   * at their tip or tail.
   */
  bool can_remove() const {
    assert(is_alive);
    return (!is_on_ch &&
            !opposite->is_on_ch &&
            v->is_of_higher_degree &&
            opposite->v->is_of_higher_degree &&
            can_remove_at_tip() &&
            opposite->can_remove_at_tip());
  }

  /** remove this edge
   */
  void remove() {
    assert(can_remove());
    assert_valid();

    v->incident_edge = opposite->prev;
    opposite->v->incident_edge = prev;


    next->prev = opposite->prev;
    prev->next = opposite->next;

    opposite->next->prev = prev;
    opposite->prev->next = next;

    prev = NULL;
    next = NULL;
    is_alive = false;
    opposite->prev = NULL;
    opposite->next = NULL;
    opposite->is_alive = false;

    v->dec_degree();
    opposite->v->dec_degree();

    assert_valid();
  }

#if 0
  /** Check whether this edge can be flipped.
   *
   * It assumes the graph is triangulated and everything is still constrained. */
  bool can_flip() const {
    if (opposite == NULL) return false;

    assert(working_set_depth == opposite->working_set_depth);

    const Vertex &tip = *v;
    const Vertex &tail = *opposite->v;
    const Vertex &left = *next->v;
    const Vertex &right = *opposite->next->v;

    assert(&left != &right);
    assert(Vertex::orientation(tail, tip, left) > 0);
    assert(Vertex::orientation(tip, tail, right) > 0);

    return ((Vertex::orientation(right, left, tail) > 0) &&
            (Vertex::orientation(left, right, tip) > 0));
  }

  /** Flip this triangulation edge.
   *
   * Note that this only works on triangulation edges and completely
   * ignores next/prev constrained pointers.  Those are assumed to
   * be fixed with reset_all_constraints() soon.
   *
   * As such, assert_valid() will fail until reset_all_constraints() is run.
   */
  void flip() {
    assert(can_flip());
    const auto triangle = [](Edge* a, Edge *b, Edge *c) {
      const auto link_edges = [](Edge* first, Edge *second) {
        first->next = second;
        second->prev = first;
      };
      link_edges(a, b);
      link_edges(b, c);
      link_edges(c, a);
    };

    Edge * l = this;     /* This edge points upwards */
    Edge * r = opposite; /* Our buddy is to our right, pointing downwards */
    Edge * tl = l->next; /* top left */
    Edge * bl = l->prev; /* bottom left */
    Edge * tr = r->prev; /* top right */
    Edge * br = r->next; /* bottom right */

    /* After the flip, this edge points from right to left and is "below" our buddy which points from left to right */
    l->v = tl->v;
    r->v = br->v;
    triangle(l, bl, br);
    triangle(r, tr, tl);
  }
#endif

  Vertex* get_tail() const { return prev->v; }

#ifndef NDEBUG
  void assert_valid() const;
#else
  void assert_valid() const {}
#endif

  friend std::ostream& operator<<(std::ostream&, const Edge&);
};
// }}}

class DCEL {
  /** Iterate around v once, returning a handle for each face.
   *
   * We walk around a vertex in clockwise order.
   */
  // {{{ Iterators
  class AroundVertexFacesIterator {
  protected:
    Edge * const e_start;
    Edge * e_cur;

  public:
    AroundVertexFacesIterator(Edge* e_vertex)
      : e_start(e_vertex)
      , e_cur(e_vertex)
    {
      assert(e_vertex);
      assert(e_vertex->is_alive);
    }

    Edge* operator*() const { return e_cur; }
    Edge* operator->() const { return e_cur; }

    AroundVertexFacesIterator& operator++() {
      assert(e_cur != NULL);
      assert(e_cur->v == e_start->v);
      assert(e_cur->is_alive);

      e_cur = e_cur->next->opposite;
      if (e_cur == e_start) {
        e_cur = NULL;
        return *this;
      };
      return *this;
    }
  };
  // }}}}}}

  /* Helper functions */
  private:
    std::vector<Vertex*> higher_degree_vertices;

  /* The state of the DCEL */
  private:
    /* These stay fixed over all iterations.
     *
     * (not necessarily their content, but * at least the set and order) */
    VertexList all_vertices; /** The list of all vertices. */
    FixedVector<Edge> all_edges; /** The list of all edges */
    unsigned num_faces; /** The number of faces of the entire graph right now. */

  /* public interface */
  public:
    /** Initialize the DCEL with the vertices and a triangulation of their CH */
    DCEL(VertexList&& vertices, const InputEdgeSet& edges);

    void write_obj_segments(bool dump_vertices, std::ostream &o);
    unsigned get_num_faces() const { return num_faces; }
    void improve_convex_decomposition();

    friend std::ostream& operator<<(std::ostream&, const DCEL&);
    friend class Vertex;

  // Debugging things:
public:
#ifndef NDEBUG
  void assert_valid() const;
#else
  void assert_valid() const {}
#endif
};

std::ostream& operator<<(std::ostream&, const Vertex&);
std::ostream& operator<<(std::ostream&, const Edge&);
std::ostream& operator<<(std::ostream&, const DCEL&);


/** Updates if this vertex is of degree 4 or higher.
 *
 * Maybe also, if it's a collinear point and of degree >=3,
 * and also if this is a vertex on the convex hull with degree > 2;
 *
 * Basically. this is to know if we should bother trying to improve things here.
 */
inline bool
Vertex::
get_vertex_is_of_higher_degree() const {
  assert(degree >= 2);

  if (degree == 2) {
    return false;
  } else if (degree >= 4) {
    return true;
  } else {
    assert(degree == 3);
    /* Degree 3 vertex */
    /* Could still have an incident angle of exactly pi; check for that */
    const Edge* e = incident_edge;
    const Vertex &v1 = *e->opposite->v;
    const Vertex &v2 = *e->next->opposite->next->v;
    const Vertex &v3 = *e->next->v;
    if (Vertex::orientation(*this, v1, v2) == 0) {
      return true;
    };
    if (Vertex::orientation(*this, v2, v3) == 0) {
      return true;
    };
    if (Vertex::orientation(*this, v3, v1) == 0) {
      return true;
    };

    /* Could be on the CH */
    if (e->is_on_ch || e->next->opposite->is_on_ch || e->opposite->prev->is_on_ch) {
      return true;
    }

    return false;
  }
}

inline void
Vertex::
assert_valid() const {
  assert(incident_edge->is_alive);
  assert(incident_edge->v == this);
  assert(get_vertex_is_of_higher_degree() == is_of_higher_degree);
  assert(degree == find_degree());
}

#ifndef NDEBUG
inline unsigned
Vertex::
find_degree() const {
  DCEL::AroundVertexFacesIterator it(incident_edge);
  unsigned d = 0;
  for (; *it; ++it) { ++d; };
  return d;
}
#endif

/* vim: set fdm=marker: */
