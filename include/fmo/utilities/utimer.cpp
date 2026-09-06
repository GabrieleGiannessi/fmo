#include <chrono>
#include <iostream>

#include "fmo/utilities/utimer.hpp"

struct utimer::Impl {
  std::chrono::steady_clock::time_point start;
  std::string message;
  long *elapsed_ms;
  bool silent;
};

utimer::utimer(const std::string &message)
    : utimer(message, nullptr, false) {}

utimer::utimer(const std::string &message, long *elapsed_ms)
    : utimer(message, elapsed_ms, false) {}

utimer::utimer(const std::string &message, long *elapsed_ms, bool silent)
    : impl(new Impl{std::chrono::steady_clock::now(), message, elapsed_ms,
                    silent}) {}

utimer::~utimer() {
  const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - impl->start).count();

  if (impl->elapsed_ms != nullptr) {
    *impl->elapsed_ms = elapsed;
  }
  if (!impl->silent) {
    std::cout << impl->message << " computed in " << elapsed << " msec"
              << std::endl;
  }
  delete impl;
};
