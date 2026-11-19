/*
 String utils

 Copyright (C) 2021 AMPL Optimization Inc

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

 Authors: Gleb Belov
 */

#include <sstream>

#include <algorithm>

#include "mp/utils-string.h"

namespace mp {

std::vector<std::string> split_string(
    const char* str) {
  std::vector<std::string> result;
  std::string word;
  std::istringstream iss(str);
  while (iss >> word) {
    if (!word.empty())
      result.push_back(word);
  }
  return result;
}

bool begins_with(std::string_view fullString, std::string_view start) {
  if (fullString.length() >= start.length())
    return (0 == fullString.compare(0, start.length(), start));
  return false;
}

bool ends_with(std::string_view fullString, std::string_view ending) {
  if (fullString.length() >= ending.length())
    return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
  return false;
}

namespace {

/// Whether a line is laid out by its author (a table): two spaces
/// after the indentation
bool IsLaidOut(std::string_view line) {
  const auto start = line.find_first_not_of(' ');
  return start != line.npos && line.find("  ", start) != line.npos;
}

/// Append \a line, wrapped at \a width, to \a out.
/// Continuation lines get the line's indentation.
void WrapLine(std::string_view line, std::size_t width, std::string& out) {
  const auto ind = std::min(line.find_first_not_of(' '), line.size());
  const std::string_view indent = line.substr(0, ind);
  std::string_view rest = line.substr(ind);
  const std::size_t room = width > ind ? width - ind : 1;
  out.append(indent);
  while (rest.size() > room) {
    auto cut = rest.rfind(' ', room);            // last space that fits
    if (cut == rest.npos || 0 == cut)            // a longer word: after it
      cut = rest.find(' ');
    if (cut == rest.npos)
      break;
    out.append(rest.substr(0, cut));
    out += '\n';
    out.append(indent);
    rest.remove_prefix(cut);
    rest.remove_prefix(std::min(rest.find_first_not_of(' '), rest.size()));
  }
  out.append(rest);
}

}  // namespace

std::string WrapLines(std::string_view text, std::size_t width) {
  width = std::max<std::size_t>(width, 1);
  std::string out;
  out.reserve(text.size() + text.size() / width + 1);
  while (true) {
    const auto eol = text.find('\n');
    const auto line = text.substr(0, eol);
    if (IsLaidOut(line))
      out.append(line);
    else
      WrapLine(line, width, out);
    if (eol == text.npos)
      break;
    out += '\n';
    text.remove_prefix(eol + 1);
  }
  return out;
}

}  // namespace mp
