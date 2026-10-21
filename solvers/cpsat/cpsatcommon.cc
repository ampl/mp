#include <cmath>
#include <limits>
#include <stdexcept>

#include "google/protobuf/descriptor.h"
#include "google/protobuf/message.h"

#include "mp/format.h"
#include "cpsatcommon.h"

namespace mp {

int CpsatCommon::NumVars() const {
  return (int)data()->vars.size();
}

int64_t CpsatCommon::LowerBoundToInt64(double lb) const {
  const double vb = data()->var_bound;
  if (lb <= -vb)
    return (int64_t)-vb;
  return (int64_t)std::ceil(lb - 1e-9);
}

int64_t CpsatCommon::UpperBoundToInt64(double ub) const {
  const double vb = data()->var_bound;
  if (ub >= vb)
    return (int64_t)vb;
  return (int64_t)std::floor(ub + 1e-9);
}


//////////////////// SatParameters access via reflection ///////////////////

namespace {

namespace pb = google::protobuf;

const pb::FieldDescriptor& FindField(
    const operations_research::sat::SatParameters& p, const char* key) {
  auto* fd = p.GetDescriptor()->FindFieldByName(key);
  if (!fd)
    throw std::runtime_error(
        fmt::format("Unknown CP-SAT parameter '{}'", key));
  return *fd;
}

}  // namespace

void CpsatCommon::GetSolverOption(const char* key, int& value) const {
  const auto& p = data()->params;
  auto& fd = FindField(p, key);
  auto* r = p.GetReflection();
  switch (fd.cpp_type()) {
  case pb::FieldDescriptor::CPPTYPE_INT32: value = r->GetInt32(p, &fd); break;
  case pb::FieldDescriptor::CPPTYPE_INT64: value = (int)r->GetInt64(p, &fd); break;
  case pb::FieldDescriptor::CPPTYPE_BOOL:  value = r->GetBool(p, &fd); break;
  case pb::FieldDescriptor::CPPTYPE_ENUM:
    value = r->GetEnumValue(p, &fd); break;
  default:
    throw std::runtime_error(fmt::format(
        "CP-SAT parameter '{}' is not an integer", key));
  }
}

void CpsatCommon::SetSolverOption(const char* key, int value) {
  auto& p = data()->params;
  auto& fd = FindField(p, key);
  auto* r = p.GetReflection();
  switch (fd.cpp_type()) {
  case pb::FieldDescriptor::CPPTYPE_INT32: r->SetInt32(&p, &fd, value); break;
  case pb::FieldDescriptor::CPPTYPE_INT64: r->SetInt64(&p, &fd, value); break;
  case pb::FieldDescriptor::CPPTYPE_BOOL:  r->SetBool(&p, &fd, value != 0); break;
  case pb::FieldDescriptor::CPPTYPE_ENUM:  r->SetEnumValue(&p, &fd, value); break;
  case pb::FieldDescriptor::CPPTYPE_DOUBLE: r->SetDouble(&p, &fd, value); break;
  default:
    throw std::runtime_error(fmt::format(
        "CP-SAT parameter '{}' is not an integer", key));
  }
}

void CpsatCommon::GetSolverOption(const char* key, double& value) const {
  const auto& p = data()->params;
  auto& fd = FindField(p, key);
  auto* r = p.GetReflection();
  switch (fd.cpp_type()) {
  case pb::FieldDescriptor::CPPTYPE_DOUBLE: value = r->GetDouble(p, &fd); break;
  case pb::FieldDescriptor::CPPTYPE_FLOAT:  value = r->GetFloat(p, &fd); break;
  default:
    throw std::runtime_error(fmt::format(
        "CP-SAT parameter '{}' is not a double", key));
  }
}

void CpsatCommon::SetSolverOption(const char* key, double value) {
  auto& p = data()->params;
  auto& fd = FindField(p, key);
  auto* r = p.GetReflection();
  switch (fd.cpp_type()) {
  case pb::FieldDescriptor::CPPTYPE_DOUBLE: r->SetDouble(&p, &fd, value); break;
  case pb::FieldDescriptor::CPPTYPE_FLOAT:  r->SetFloat(&p, &fd, (float)value); break;
  default:
    throw std::runtime_error(fmt::format(
        "CP-SAT parameter '{}' is not a double", key));
  }
}

void CpsatCommon::GetSolverOption(const char* key, std::string& value) const {
  const auto& p = data()->params;
  auto& fd = FindField(p, key);
  if (fd.cpp_type() != pb::FieldDescriptor::CPPTYPE_STRING)
    throw std::runtime_error(fmt::format(
        "CP-SAT parameter '{}' is not a string", key));
  value = p.GetReflection()->GetString(p, &fd);
}

void CpsatCommon::SetSolverOption(const char* key, const std::string& value) {
  auto& p = data()->params;
  auto& fd = FindField(p, key);
  if (fd.cpp_type() != pb::FieldDescriptor::CPPTYPE_STRING)
    throw std::runtime_error(fmt::format(
        "CP-SAT parameter '{}' is not a string", key));
  p.GetReflection()->SetString(&p, &fd, value);
}

} // namespace mp
