#pragma once

#include <string>

class utimer {
public:
  explicit utimer(const std::string &message);
  utimer(const std::string &message, long *elapsed_ms);
  utimer(const std::string &message, long *elapsed_ms, bool silent);
  ~utimer();

  utimer(const utimer &) = delete;
  utimer &operator=(const utimer &) = delete;

private:
  struct Impl;
  Impl *impl;
};
