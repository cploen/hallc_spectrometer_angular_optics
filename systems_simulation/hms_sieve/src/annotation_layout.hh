#ifndef HMS_ANNOTATION_LAYOUT_HH
#define HMS_ANNOTATION_LAYOUT_HH

// Display metadata only. No coordinates, transforms, or ray equations here.
#include <array>
#include <string>

namespace hms_display {
enum class Layer { base, labels, lab, hms, sieve, code, branches, rays };
struct LayerName { Layer layer; const char* name; };
inline constexpr std::array<LayerName, 8> layers{{
  {Layer::base, "HMS_base"}, {Layer::labels, "HMS_labels"},
  {Layer::lab, "HMS_lab_frame"}, {Layer::hms, "HMS_transport_frame"},
  {Layer::sieve, "HMS_sieve_frame"}, {Layer::code, "HMS_code_names"},
  {Layer::branches, "HMS_branch_names"}, {Layer::rays, "HMS_rays"}
}};

inline bool starts(const std::string& name, const std::string& prefix) {
  return name.compare(0, prefix.size(), prefix) == 0;
}

// Every existing primitive belongs to exactly one drawing layer. Annotations
// are separate actions; changing visibility never moves or replaces a point.
inline Layer primitiveLayer(const std::string& name) {
  if (starts(name, "LAB_")) return Layer::lab;
  if (starts(name, "HMS_TRANSPORT_")) return Layer::hms;
  if (starts(name, "SIEVE_LOCAL_")) return Layer::sieve;
  if (starts(name, "A_endpoint") || starts(name, "B_fit_target") ||
      starts(name, "C_HCANA_coordinate_step") || name == "synthetic_ray_vertex")
    return Layer::rays;
  return Layer::base;
}

struct QuantityName {
  const char* conceptTag;
  const char* code;
  const char* branch;
  const char* meaning;
};
// One row = one stored quantity. Both naming layers use these SAME rows.
// Source: GEOMETRY_AUDIT.md, saved-versus-constructed branch provenance table.
// These are a naming key, NEVER labels claiming synthetic A/B/C are event data.
inline constexpr std::array<QuantityName, 11> quantities{{
  {"T", "xtar", "H.gtr.x", "target x intercept"},
  {"T", "ytar", "H.gtr.y", "target y intercept"},
  {"T", "xptar", "H.gtr.th", "target x slope"},
  {"T", "yptar", "H.gtr.ph", "target y slope"},
  {"S", "xsieve", "H.extcor.xsieve", "projected x"},
  {"S", "ysieve", "H.extcor.ysieve", "projected y"},
  {"V", "reactx", "H.react.x", "reaction lab X"},
  {"V", "reacty", "H.react.y", "reaction lab Y"},
  {"V", "reactz = ztar", "H.react.z", "reaction lab Z"},
  {"R", "xbpm_tar", "H.rb.raster.fr_xbpm_tar", "raster/BPM X; reference OPEN"},
  {"R", "ybpm_tar", "H.rb.raster.fr_ybpm_tar", "raster/BPM Y; reference OPEN"}
}};

inline std::string codeLabel(const QuantityName& q) {
  return std::string("[") + q.conceptTag + "] " + q.code + "  (" + q.meaning + ")";
}
inline std::string aliasLabel(const QuantityName& q) {
  return std::string("[") + q.conceptTag + "] " + q.code + " = " + q.branch;
}
}  // namespace hms_display
#endif
