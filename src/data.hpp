#include "tinyjson.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

#include <string>
#include <vector>

using std::cerr;
using std::endl;
using std::ifstream;
using std::ostringstream;
using std::string;
using std::vector;

// 定义相对分子质量
#define M_Fe 55.85
#define M_O 16.
#define M_P 31.
#define M_P2O5 142.
#define M_Si 28.
#define M_SiO2 60.
#define M_FeO 71.85
#define M_Fe2O3 159.7
#define M_FeS 87.85
#define M_Ti 48.
#define M_TiO2 80.
#define M_Mn 55.
#define M_MnO 71.
#define M_MnO2 87.
#define M_V 56.
#define M_V2O5 182.
#define M_C 12.
#define M_CH4 16.
#define M_S 32.
#define M_CO2 44.
#define M_H2O 18.
#define M_CO 28.
#define M_N2 28.
#define M_H 1.
#define M_H2 2.
#define IdealR 22.4
#define M_Fe2SiO4 232.

// 焦炭与喷煤热值
#define coalHeatp 1172
#define blowHeat 1005

#define CpO2 1.510
#define CpN2 1.417
#define CpH2O 1.784
#define CpCO 1.31
#define CpCO2 1.807
#define CpH2 1.297
 	 	



inline string readFileIntoString(const string &path) {
  auto ss = ostringstream{};
  ifstream input_file(path);
  if (!input_file.is_open()) {
    cerr << "Could not open the file - '" << path << "'" << endl;
    exit(EXIT_FAILURE);
  }
  ss << input_file.rdbuf();
  return ss.str();
}

inline tiny::TinyJson get_json(tiny::TinyJson org, string key) {
  string ret_str = org.Get<string>(key);
  tiny::TinyJson ret;
  ret.ReadJson(ret_str);
  return ret;
}

inline tiny::TinyJson get_json_from_xobj(tiny::xobject xobj, string key) {
  tiny::TinyJson ret;
  for (int i = 0; i < xobj.Count(); i++) {
    xobj.Enter(i);
    auto ret_str = xobj.Get<string>(key);
    ret.ReadJson(ret_str);
  }
  return ret;
}

inline double get_double_from_xobj(tiny::xobject xobj, string key) {
  double ret;
  for (int i = 0; i < xobj.Count(); i++) {
    xobj.Enter(i);
    ret = xobj.Get<double>(key);
  }
  return ret;
}

struct Burden_Element {
  double OMEGA;
  tiny::TinyJson ELEMENT;
  tiny::TinyJson CONTENT;
  double MOIST;
  Burden_Element() = default;
  ~Burden_Element() = default;
  double get_element_volumn(string target) {
    return 0.01 * OMEGA * ELEMENT.Getd(target);
  }
  double get_content_volumn(string target) {
    return 0.01 * OMEGA * CONTENT.Getd(target);
  }
  static Burden_Element build_burden_element(tiny::xobject xobj) {
    return Burden_Element{get_double_from_xobj(xobj, "OMEGA"),
                          get_json_from_xobj(xobj, "ELEMENT"),
                          get_json_from_xobj(xobj, "CONTENT"),
                          get_double_from_xobj(xobj, "MOIST")};
  }
};

struct Fuel_Element {
  double OMEGA;
  double CF;
  double S;
  tiny::TinyJson CONTENT;
  double MOIST;
  Fuel_Element() = default;
  ~Fuel_Element() = default;
  double get_content_volumn(string target) {
    return 0.01 * OMEGA * CONTENT.Getd(target);
  }
  static Fuel_Element build_fuel_element(tiny::xobject xobj) {
    return Fuel_Element{
        get_double_from_xobj(xobj, "OMEGA"), get_double_from_xobj(xobj, "CF"),
        get_double_from_xobj(xobj, "S"), get_json_from_xobj(xobj, "CONTENT"),
        get_double_from_xobj(xobj, "MOIST")};
  }
};

struct Burdens {
  vector<Burden_Element> ores;
  Burden_Element partical_Fe;
  Burden_Element lime;
  Burden_Element dust;
  std::vector<Burden_Element> others;
  Burdens() = default;
  ~Burdens() = default;
  double get_ores_element_volumn(string elementName) {
    double oreContent = {};
    for (auto ore : ores) {
      oreContent += ore.get_content_volumn(elementName);
    }
    return oreContent;
  }
  double get_ores_content_volumn(string elementName) {
    double oreContent = {};
    for (auto ore : ores) {
      oreContent += 0.01 * ore.OMEGA * ore.CONTENT.Getd(elementName);
    }
    return oreContent;
  }
  double get_total_input_content(string contentName) {
    return get_ores_content_volumn(contentName) +
           partical_Fe.get_content_volumn(contentName) +
           lime.get_content_volumn(contentName) -
           dust.get_content_volumn(contentName);
  }
};

struct Diffusion {
  tiny::TinyJson HOT_METAL;
  tiny::TinyJson SLAG;
  tiny::TinyJson GAS;
  Diffusion() = default;
  ~Diffusion() = default;
  static Diffusion build_diffusion(tiny::xobject xobj) {
    return Diffusion{
        get_json_from_xobj(xobj, "HOT_METAL"),
        get_json_from_xobj(xobj, "SLAG"),
        get_json_from_xobj(xobj, "GAS"),
    };
  }
};

struct Condition {
  Diffusion diffusion;
  tiny::TinyJson hot_metal;
  tiny::TinyJson org;
  Condition() = default;
  ~Condition() = default;
};
