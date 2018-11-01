// C++ standard library includes
#include <string>
#include <vector>
#include <set>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <random>


// CAF includes
#include "caf/all.hpp"
#include "caf/io/all.hpp"

// Boost includes
CAF_PUSH_WARNINGS
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/random.hpp>
#include <boost/random.hpp>
CAF_POP_WARNINGS

// Own includes
#include "is_probable_prime.hpp"
#include "int512_serialization.hpp"

using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::unordered_map;
using namespace std::chrono;

using boost::multiprecision::int512_t;
using boost::multiprecision::gcd;
using boost::random::seed_seq;
typedef boost::random::independent_bits_engine<boost::random::mt19937, 512,
                                               boost::multiprecision::uint512_t>
                                               generator512_type;

using namespace caf;

using calc_atom = atom_constant<atom("calc")>;
using start_calc_atom = atom_constant<atom("start")>;
using result_atom = atom_constant<atom("result")>;
using quit_atom = atom_constant<atom("quit")>;
using spawn_actor_atom = atom_constant<atom("spawn")>;
using new_actor_atom = atom_constant<atom("new_actor")>;

struct config : actor_system_config {
  string host = "localhost";
  uint16_t port = 3334;
  size_t num_workers = 0;
  string mode = "coordinator";
  config() {
    add_message_type<std::pair<int512_t, int512_t>>("std::pair<int512_t, int512_t>");
    add_message_type<int512_t>("int512_t");
    add_message_type<vector<int512_t>>("vector<int512_t>");
    opt_group{custom_options_, "global"}
    .add(host, "host,H", "server host (ignored in server mode)")
    .add(port, "coordinator-port,p", "coordinator port")
    .add(num_workers, "num-workers,w", "number of workers (in manager mode)")
    .add(mode, "mode,m", "one of 'coordinator', 'manager' or 'client'");
  }
};

namespace worker {
  int512_t modular_pow(int512_t base, int exponent,
                       int512_t modulus) {
    int512_t result = 1;
    while (exponent > 0) {
      if (exponent & 1)
        result = (result * base) % modulus;
      exponent = exponent >> 1;
      base = (base * base) % modulus;
    }
    return result;
  }

  struct state {
    generator512_type gen;
    int512_t n = 0;
    int512_t x = 0;
    int512_t y = 0;
    int512_t c = 0;
    int512_t d = 0;
  };

  void rho(stateful_actor<state>* self, int512_t n, actor ret) {
    if(!self->mailbox().empty()) return;
    auto& state = self->state;
    cout << "got task " << n << endl;

    if (n == 1) {
      if(!self->mailbox().empty()) return;
      self->send(ret, result_atom::value, std::make_pair(n, int512_t("1")), self);
      aout(self) << "sent result n == 1" << endl;
      return;
    }
    if((n%2) == 0) {
      if(!self->mailbox().empty()) return;
      self->send(ret, result_atom::value, std::make_pair(n, int512_t("2")), self);
      aout(self) << "sent result " << n << " n%2" << endl;
      return;
    }

    state.n = n;
    state.x = (state.gen() % (state.n - 2)) + 2;
    state.y = state.x;
    state.c = (state.gen() % (state.n - 1)) + 1;
    state.d = 1;

    if(!self->mailbox().empty()) return;

    while (true) {
      state.x = (modular_pow(state.x, 2, state.n) + state.c + state.n) % state.n;
      state.y = (modular_pow(state.y, 2, state.n) + state.c + state.n) % state.n;
      state.y = (modular_pow(state.y, 2, state.n) + state.c + state.n) % state.n;
      state.d = gcd(abs(state.x - state.y), state.n);

      if (state.d == state.n) {
        rho(self, state.n, ret);
        return ;
      }

      if(!self->mailbox().empty()) return;

      if(state.d != 1) {
        // only send result if it is prime
        self->send(ret, result_atom::value, std::make_pair(n, state.d), self);
        cout << "sent result " << state.d << endl;
        return;
      }
    }
  }

  /**
   * actor calculates a prime factor for the given number n
   * @return 0 if calculation was postponed!!!
   */
  behavior rho_actor(stateful_actor<state>* self, unsigned int seed) {
    std::default_random_engine gen;
    gen.seed(seed);
    std::uniform_int_distribution<unsigned int> dist(0, std::numeric_limits<unsigned int>::max());
    seed_seq ss = {dist(gen), dist(gen), dist(gen), dist(gen), dist(gen), dist(gen), dist(gen)};
    self->state.gen.seed(ss);

    return {
      [=](calc_atom, int512_t n, actor ret) {
        rho(self, n, ret);
      },

      // ignore handlers
      [=](quit_atom) {}
    };
  }
}


namespace coordinator {
  struct state {
    std::vector<actor> workers;
    int512_t number;
    std::set<int512_t> problem_parts;
    std::vector<int512_t> prime_factors;
    actor requester;
    bool running = false;
  };

  bool isContained(std::vector<int512_t>& vec, int512_t& value) {
    return (std::find(vec.begin(), vec.end(), value) != vec.end());
  }

  int divide(int512_t& number, int512_t& factor) {
    int i = 0;
    while ((number % factor) == 0) {
      number /= factor;
      ++i;
    }
    return i;
  }

  void insert(stateful_actor<state>* self, int512_t& div) {
    if (div == 1) {
      // nop
    } else if(is_probable_prime(div)) {
      // needs to be like this
      if (!isContained(self->state.prime_factors, div)) {
        auto amount = divide(self->state.number, div);
        std::fill_n(std::back_inserter(self->state.prime_factors), amount, div);
        cout << "filled in " << amount << " * " << div << endl;
      }
    } else {
      self->state.problem_parts.insert(div);
      cout << "new problem part " << div << endl;
    }
  }

  void checkProblems(stateful_actor<state>* self, int512_t factor) {
    std::vector<int512_t> del;
    std::vector<int512_t> ins;

    for (auto p : self->state.problem_parts) {
      while ((p != factor) && ((p % factor) == 0)) {
        del.push_back(p);
        p /= factor;
      }
      ins.push_back(p);
    }

    for(auto& d : del) {
      self->state.problem_parts.erase(d);
    }
    for(auto& i : ins) {
      insert(self, i);
    }
  }

  behavior coordinator_actor(stateful_actor<state>* self) {
    self->set_down_handler([&](const down_msg& msg) {
      self->state.workers.erase(std::remove(self->state.workers.begin(),
                                self->state.workers.end(), msg.source),
                                self->state.workers.end());
    });

    return {
      [=](start_calc_atom, actor& requester, int512_t n) {
        self->state.running = true;
        self->state.prime_factors.clear();
        self->state.problem_parts.clear();
        self->state.number = n;
        self->state.requester = requester;
        self->state.problem_parts.insert(n);
        for(auto& w : self->state.workers) {
          self->send(w, calc_atom::value, n, self);
        }
      },
      [=](result_atom res_atom, std::pair<int512_t, int512_t>& res, actor sender){
        if(!self->state.running) return;
        auto& task = res.first;
        auto& factor = res.second;
        auto& problem_parts = self->state.problem_parts;
        auto& prime_factors = self->state.prime_factors;

        // task is "solved"
        problem_parts.erase(task);

        // check factor and remainder and insert it if possible
        insert(self, factor);
        int512_t remainder = task / factor;
        insert(self, remainder);

        // now try to divide open problems by factor
        checkProblems(self, factor);

        if (problem_parts.empty()) {
          self->state.running = false;
          self->send(self->state.requester, res_atom, prime_factors);
          for (auto& w : self->state.workers) {
            self->send(w, quit_atom::value);
          }
        }
        else {
          self->send(sender, calc_atom::value, *problem_parts.begin(), self);
        }
      },
      [=](new_actor_atom, std::vector<actor>& workers) {
        self->state.workers.insert(self->state.workers.end(), workers.begin(), workers.end());
        for(auto& w : workers) {
          // TODO server segfaults on downhandler // self->monitor(w);
          // start new actors with task
          if (!self->state.problem_parts.empty()) {
            self->send(w, calc_atom::value, *self->state.problem_parts.begin(), self);
          }
        }
      },
      [](calc_atom, int512_t, actor) {/*ignore*/},
      [](quit_atom) {/*ignore*/}
    };
  }

  void run(actor_system& sys, const config& cfg) {
    cout << "coordinator started" << endl;
    auto coord = sys.spawn(coordinator::coordinator_actor);

    // publish this coordinator
    io::publish(coord, cfg.port, nullptr, true);

    cout << "press [enter] to quit" << endl;
    std::string dummy;
    getline(cin, dummy);

    anon_send_exit(coord, exit_reason::user_shutdown);
  }
}


namespace manager {
  struct state {
    std::vector<actor> workers;
  };

  behavior manager_actor(stateful_actor<state>* self, const actor& coordinator) {
    self->monitor(coordinator);
    srand(time(nullptr));

    self->set_exit_handler([=](const exit_msg& msg){
      for (auto& w : self->state.workers) {
        self->send_exit(w, msg.reason);
      }
    });
    self->set_down_handler([=](const down_msg&) {
      aout(self) << "lost connection to server" << std::endl;
    });

    return {
      [=](spawn_actor_atom, size_t amount) {
        std::vector<actor> workers;
        for(size_t i = 0; i < amount; ++i) {
          auto w = self->spawn(worker::rho_actor, rand());
          workers.push_back(w);
        }
        self->send(coordinator, new_actor_atom::value, workers);
        self->state.workers.insert(self->state.workers.end(), workers.begin(), workers.end());
      }
    };
  }

  void run(actor_system& sys, const config& cfg) {
    auto coordinator = sys.middleman().remote_actor(cfg.host, cfg.port);
    if (!coordinator) {
      cerr << "unable to connect to coordinator: "
           << sys.render(coordinator.error()) << std::endl;
      return;
    }
    cout << "successfully connected to coordinator" << endl;

    auto man = sys.spawn(manager_actor, *coordinator);
    anon_send(man, spawn_atom::value, cfg.num_workers);

    cout << "press [enter] to quit" << endl;
    std::string dummy;
    getline(cin, dummy);

    anon_send_exit(man, exit_reason::user_shutdown);
  }
}


namespace client {
  struct state {
    milliseconds start;
    int512_t number;
  };

  behavior client_actor(stateful_actor<state>* self, const actor& coordinator) {
    self->monitor(coordinator);
    self->set_down_handler([=](const down_msg&) {
      aout(self) << "lost connection to server" << endl;
    });

    return {
      // blocks this client actor until the result is received.
      [=](start_calc_atom calc, int512_t n) {
        self->state.start = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
        self->send(coordinator, calc, self, n);
        self->state.number = n;
      },
      [=](result_atom, std::vector<int512_t>& factors) {
        milliseconds duration = duration_cast<milliseconds>(system_clock::now().time_since_epoch()) - self->state.start;
        std::sort(factors.begin(), factors.end());
        std::stringstream ss;
        ss << std::setprecision(std::numeric_limits<int512_t>::max_digits10)
             << self->state.number << " = ";
        for(auto& fac : factors) {
          ss << fac << " * ";
        }
        string out = ss.str();
        out.resize(out.size()-3);
        aout(self) << out << endl;
        aout(self) << "in " << duration << endl;
      }
    };
  }

  void run(actor_system& sys, const config& cfg) {
    auto coordinator = sys.middleman().remote_actor(cfg.host, cfg.port);
    if (!coordinator) {
      cerr << "unable to connect to coordinator: "
           << sys.render(coordinator.error()) << std::endl;
      return;
    }
    cout << "successfully connected to coordinator" << endl;

    auto client = sys.spawn(client_actor, *coordinator);
    cout << "client started\n"
            "/quit - quits the client program\n"
            "/connect <host>:<port> - connect to given address\n"
            "/calc <int512_t> - calculate primefactors for given number"
            << endl;

    string line;
    while (std::getline(cin, line)) {
      std::vector<string> parts;
      split(parts, line, is_any_of(" "), token_compress_on);

      if (!parts.empty()) {
        if (parts[0] == "/quit") {
          break;
        } else if (parts[0] == "/calc" && parts.size() == 2) {
          anon_send(client, start_calc_atom::value, int512_t(parts[1]));
        }
        parts.clear();
      }
    }

    anon_send(client, exit_reason::user_shutdown);
  }
}

// dispatches to run_* function depending on selected mode
void caf_main(actor_system &sys, const config& cfg) {
  using map_t = unordered_map<string, void (*)(actor_system &, const config &)>;
  map_t modes{
      {"coordinator", coordinator::run},
      {"manager", manager::run},
      {"client", client::run}
  };
  auto i = modes.find(cfg.mode);
  if (i != modes.end())
    (i->second)(sys, cfg);
  else
    cerr << "*** invalid mode specified" << endl;
}

CAF_MAIN(io::middleman)
