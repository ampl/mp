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

std::pair<int64_t, int64_t>
CpsatCommon::VarDomainToInt64(double lb, double ub) const {
  constexpr double kMax = kCpsatMaxVarBound;
  if (lb > ub)                  // contradictory: an empty domain
    return { 1, 0 };
  if (lb > kMax || ub < -kMax)
    throw std::runtime_error(fmt::format(
        "CP-SAT: variable {} has bounds [{:.17g}, {:.17g}] beyond the "
        "range CP-SAT can represent (+-{:g})", NumVars(), lb, ub, kMax));
  const double vb = data()->var_bound;
  const bool lb_inf = !std::isfinite(lb), ub_inf = !std::isfinite(ub);
  const double t = data()->inteps;
  double l = lb_inf ? -vb : RoundBoundInward(lb, true, t);
  double u = ub_inf ? vb : RoundBoundInward(ub, false, t);
  if (lb_inf && !ub_inf && u <= -vb)   // keep the width var_bound
    l = u - vb;
  if (ub_inf && !lb_inf && l >= vb)
    u = l + vb;
  l = std::max(l, -kMax);
  u = std::min(u, kMax);
  return { static_cast<int64_t>(l), static_cast<int64_t>(u) };
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
