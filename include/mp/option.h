/*
 Command-line option parser.

 Copyright (C) 2014 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#ifndef MP_OPTION_H_
#define MP_OPTION_H_

#include <vector>

#include "mp/error.h"

namespace mp {

// An option error.
class OptionError : public Error {
public:
  explicit OptionError(fmt::StringRef message) : Error(message) {}
};

struct Option {
  char name;
  const char *description;
  void *handler;
  bool (*on_option)(void *);
};

// A list of options.
class OptionList {
 private:
  typedef std::vector<Option> OptionContainer;
  OptionContainer options_;
  bool sorted_;

  template <typename Handler, bool (Handler::*on_option)()>
  static bool OnOption(void *handler) {
    return (static_cast<Handler*>(handler)->*on_option)();
  }

  template <typename Handler, bool (Handler::*on_option)()>
  void Add(char name, const char *description, Handler &h) {
    Option option = {name, description, &h, OnOption<Handler, on_option>};
    options_.push_back(option);
    sorted_ = false;
  }

 public:
  template <typename Handler>
  class Builder {
   private:
    OptionList &options_;
    Handler &handler_;

   public:
    Builder(OptionList &options, Handler &h)
      : options_(options), handler_(h) {}

    // Adds an option.
    // on_option: method called when option has been parsed; returns true
    //            to continue parsing, false to stop
    template <bool (Handler::*on_option)()>
    void Add(char name, const char *description) {
      options_.Add<Handler, on_option>(name, description, handler_);
    }
  };

  explicit OptionList() : sorted_(true) {}

  typedef OptionContainer::const_iterator iterator;

  iterator begin() const { return options_.begin(); }
  iterator end() const { return options_.end(); }

  bool sorted() const { return sorted_; }

  // Sorts the option list.
  void Sort();

  // Finds an option with the specified name in the list.
  // Requires list to be sorted.
  const Option *Find(char name) const;
};

// Parses command-line options until the first argument that doesn't
// start with '-'. Returns the option that terminated parsing or 0 if
// parsing continued till the end.
char ParseOptions(char **&args, OptionList &options);
}  // namespace mp

#endif  // MP_OPTION_H_
