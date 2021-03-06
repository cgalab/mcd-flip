/* CG:SHOP 2020: Minimum Convex Decomposition -- Flip Tool
*
*  Copyright 2019, 2020 Peter Palfraader
*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "io.h"
#include "geom.h"

#include "gitversion.h"

#include <fstream>
#include <iostream>
#include <chrono>
#include <getopt.h>
#include <csignal>

INITIALIZE_EASYLOGGINGPP
unsigned DBG_INDENT_CTR = 0;
std::default_random_engine random_engine;
bool main_loop_interrupted = false;
const double DCEL::default_move_freedom_in_direction_probability = 0.60;
const unsigned DCEL::default_move_freedom_in_direction_new_pick_ctr = 250000;
const unsigned DCEL::default_move_distance_prob_bound = 40;

/*seconds*/

static void
setup_logging(int argc, char* argv[]) {
  START_EASYLOGGINGPP(argc, argv);

  el::Configurations defaultConf;
  defaultConf.setGlobally(el::ConfigurationType::Format, "%datetime{%H%m%s.%g} %levshort %msg");
  el::Loggers::reconfigureAllLoggers(defaultConf);
  el::Loggers::addFlag( el::LoggingFlag::ColoredTerminalOutput );
}

[[noreturn]]
static void
usage(const char *progname, int err) {
  std::ostream &f = err ? std::cerr : std::cout;

  f << "Usage: " << progname << "[options] <INPUT> <OUTPUT>" << std::endl
    << "  Options" << std::endl
    << "    --seed NUM             seed of the RNG" << std::endl
    << "    --full-obj             also print vertex coordinates to .obj file" << std::endl
    << "    --lower-bound NUM      return immediately when this is reached" << std::endl
    << "    --improve-time SECONDS After improving, keep working for at least this long" << std::endl
    << "    --max-time NUM         Do not start a new run after NUM seconds (overrides improve-* bounds)" << std::endl
    << "    --log-interval SECONDS Report on state regularly." << std::endl
    << "    --move_freedom_in_direction_probability"  " (default: " << DCEL::default_move_freedom_in_direction_probability << ")" << std::endl
    << "    --move_freedom_in_direction_new_pick_ctr" " (default: " << DCEL::default_move_freedom_in_direction_new_pick_ctr << ")" << std::endl
    << "    --move_distance_prob_bound"               " (default: " << DCEL::default_move_distance_prob_bound << ")" << std::endl
    << std::endl
  ;
  exit(err);
}


static void
signalHandler( int signum ) {
   LOG(INFO) << "Interrupt signal (" << signum << ") received.\n";
   main_loop_interrupted = true;
}

int main(int argc, char *argv[]) {
  const char * const short_options = "hS:fFB:i:T:L:";
  const option long_options[] = {
    { "help"        , no_argument      , 0, 'h'},
    { "seed"        , required_argument, 0, 'S'},
    { "full-obj"    , no_argument      , 0, 'f'},
    { "lower-bound" , required_argument, 0, 'B'},
    { "improve-time", required_argument, 0, 'i'},
    { "max-time"    , required_argument, 0, 'T'},
    { "log-interval", required_argument, 0, 'L'},
    { "move_freedom_in_direction_probability", required_argument, 0, '1'},
    { "move_freedom_in_direction_new_pick_ctr", required_argument, 0, '2'},
    { "move_distance_prob_bound", required_argument, 0, '3'},
    { 0, 0, 0, 0}
  };

  setup_logging(argc, argv);

  long requested_seed = 0;
  bool full_obj = false;
  unsigned lower_bound = 0;
  unsigned improvement_time = 1;
  int max_time = 3;
  unsigned log_interval = 60;
  double move_freedom_in_direction_probability = DCEL::default_move_freedom_in_direction_probability;
  unsigned move_freedom_in_direction_new_pick_ctr = DCEL::default_move_freedom_in_direction_new_pick_ctr;
  unsigned move_distance_prob_bound = DCEL::default_move_distance_prob_bound;

  while (1) {
    int option_index = 0;
    int r = getopt_long(argc, argv, short_options, long_options, &option_index);

    if (r == -1) break;
    switch (r) {
      case 'h':
        usage(argv[0], 0);
        break;

      case 'S':
        requested_seed = atol(optarg);
        break;

      case 'f':
        full_obj = true;
        break;

      case 'B':
        lower_bound = atol(optarg);
        break;

      case 'i':
        improvement_time = atol(optarg);
        break;

      case 'T':
        max_time = atol(optarg);
        break;

      case 'L':
        log_interval = atol(optarg);
        break;

      case '1':
        move_freedom_in_direction_probability = std::stod(optarg);
        break;

      case '2':
        move_freedom_in_direction_new_pick_ctr = atol(optarg);
        break;

      case '3':
        move_distance_prob_bound = atol(optarg);
        break;

      default:
        std::cerr << "Invalid option " << (char)r << std::endl;
        exit(1);
    }
  }

  if (argc - optind > 2) {
    usage(argv[0], 1);
  }

  std::istream *in = &std::cin;
  std::ostream *out = &std::cout;
  std::ifstream filestreamin;
  std::ofstream filestreamout;

  if (argc - optind >= 1) {
    std::string fn(argv[optind]);
    if (fn != "-") {
      filestreamin.open(fn);
      in = &filestreamin;
    }
  }
  if (argc - optind >= 2) {
    std::string fn(argv[optind + 1]);
    if (fn != "-") {
      filestreamout.open(fn);
      out = &filestreamout;
    }
  }

  signal(SIGINT, signalHandler);
  signal(SIGTERM, signalHandler);
  signal(SIGUSR1, signalHandler);
  signal(SIGUSR2, signalHandler);

  auto start_time = std::chrono::system_clock::now();
  auto end_time = start_time + std::chrono::seconds(max_time);
  auto last_info_time = std::chrono::system_clock::now();
  auto info_interval = std::chrono::seconds(log_interval);
  std::chrono::time_point<std::chrono::system_clock> solution_found_at;

  unsigned num_iters = 0;
  unsigned num_iters_since_improved = 0;
  int best_seed = 0;

  std::random_device real_rng("/dev/urandom");
  int seed = requested_seed != 0 ? requested_seed : real_rng();
  std::cout << "version: " << GITVERSION << std::endl;
  std::cout << "random_seed: " << seed << std::endl << std::flush;
  random_engine.seed(seed);

  std::unique_ptr<DCEL> decl;
  std::pair<VertexList, InputEdgeSet> p = load_obj(*in);
  decl = std::make_unique<DCEL>(
    std::move(p.first),
    p.second,
    move_freedom_in_direction_probability,
    move_freedom_in_direction_new_pick_ctr,
    move_distance_prob_bound);

  unsigned initial_to_beat = decl->get_num_faces();
  unsigned current_num_faces = initial_to_beat;
  while (1) {
    decl->assert_valid();
    decl->improve_convex_decomposition();
    decl->assert_valid();

    ++num_iters;
    ++num_iters_since_improved;

    auto now = std::chrono::system_clock::now();

    if (decl->get_num_faces() < current_num_faces) {
      num_iters_since_improved = 0;
      solution_found_at = now;
      current_num_faces = decl->get_num_faces();
    }

    if (now > last_info_time + info_interval) {
      if (current_num_faces < initial_to_beat) {
        LOG(INFO) << "Iter " << num_iters << " overall and " << num_iters_since_improved << " since improved; Current num faces " << current_num_faces << "; initial to beat: " << initial_to_beat;
      } else {
        LOG(INFO) << "Iter " << num_iters << " overall; to beat: " << initial_to_beat;
      }
      last_info_time = now;
    }

    if (UNLIKELY(current_num_faces <= lower_bound)) {
      LOG(INFO) << "We hit the lower bound of " << lower_bound;
      std::cout << "exit_reason: lower_bound" << std::endl;
      break;
    } else if (UNLIKELY((current_num_faces < initial_to_beat) && (now > solution_found_at + std::chrono::seconds(improvement_time)))) {
      LOG(INFO) << "We ran for " << num_iters << " overall and " << num_iters_since_improved << " since improved.  We did improve on " << initial_to_beat << " faces by " << (initial_to_beat - current_num_faces);
      std::cout << "exit_reason: found-after-improvement_time" << std::endl;
      std::cout << "improvement_time: " << improvement_time << std::endl;
      break;
    } else if (UNLIKELY(max_time != 0 && now > end_time)) {
      LOG(INFO) << "We ran for max-time of " << max_time << " seconds.  (We did " << num_iters << " iterations total.)";
      std::cout << "exit_reason: timeout" << std::endl;
      break;
    } else if (UNLIKELY(main_loop_interrupted)) {
      std::cout << "exit_reason: interrupt" << std::endl;
      break;
    }
  }

  DBG(DBG_GENERIC) << "Random seed was " << seed;
  std::cout << "num_cvx_areas: " << decl->get_num_faces() << std::endl;
  std::cout << "num_iters: " << num_iters << std::endl;
  if (current_num_faces < initial_to_beat) {
    auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - solution_found_at);
    std::cout << "time_since_found: " << milliseconds.count()/1000. << std::endl;
  }
  auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
  std::cout << "run_time: " << milliseconds.count()/1000. << std::endl;
  decl->write_obj_segments(full_obj, *out);
  return 0;
}
