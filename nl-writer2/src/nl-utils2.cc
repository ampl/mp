#include "mp/nl-utils2.h"

/////////////////////// TECHNICAL /////////////////////////
namespace mp {

void File::Printf(const char *format, ...) {
  va_list args;
  va_start (args, format);
  std::vfprintf (f_, format, args);
  va_end (args);
}

File NLUtils::
openf(const std::string& fname, int Close, const char *mode)
{
  File f;
  if (Close) {
    std::remove(fname.c_str());
    return f;
  }
  f.Open(fname.c_str(), mode);
  if (!f) {
    log_warning("can't open %s", fname.c_str());
    return f;
  }
  if (if_show_filenames())
    log_message("File %s\n", fname.c_str());
  return f;
}

void NLUtils::log_message(const char *format, ...) {
  va_list args;
  va_start (args, format);
  std::vprintf (format, args);
  va_end (args);
}

void NLUtils::log_warning(const char *format, ...) {
  std::fprintf(stderr, "WARNING: ");
  va_list args;
  va_start (args, format);
  std::vfprintf (stderr, format, args);
  va_end (args);
  std::fprintf(stderr, "\n");
}

}  // namespace mp
