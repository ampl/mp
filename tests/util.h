#ifndef TESTS_UTIL_H_
#define TESTS_UTIL_H_

#include <string>

std::string ReadFile(const char *name);
void WriteFile(const char *name, const char *data);

#endif  // TESTS_UTIL_H_
