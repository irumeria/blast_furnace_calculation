#include "data.hpp"
#include "matrix.hpp"
#include <array>
#include <cmath>
#include <cstddef>
#include <map>

using namespace std;

struct BF_Info {
  Burdens burdens;
  Fuel_Element coke;
  Fuel_Element coal;
  Condition condition;
  vector<string> oresLabel;
  map<string, double> slag;
  map<string, double> blast;
  map<string, double> topgas;
  map<string, double> total_heat;
  map<string, double> area_heat;
  double RO_min;
  double RO;
  BF_Info() = default;
  ~BF_Info() = default;

  BF_Info &processData(tiny::TinyJson &root, vector<string> oreNames) {
    auto burden_xo = root.Get<tiny::xobject>("BURDEN");
    oresLabel = oreNames;
    for (auto i = 0; i < burden_xo.Count(); i++) {
      burden_xo.Enter(i);
      for (auto ore : oreNames) {
        burdens.ores.push_back(Burden_Element::build_burden_element(
            burden_xo.Get<tiny::xobject>(ore)));
      }

      auto aust_xobj = burden_xo.Get<tiny::xobject>("aust");
      auto partical_Fe_xobj = burden_xo.Get<tiny::xobject>("PARTICAL_FE");
      auto lime_xobj = burden_xo.Get<tiny::xobject>("LIME");
      auto dust_xobj = burden_xo.Get<tiny::xobject>("DUST");
      burdens.partical_Fe =
          Burden_Element::build_burden_element(partical_Fe_xobj);
      burdens.lime = Burden_Element::build_burden_element(lime_xobj);
      burdens.dust = Burden_Element::build_burden_element(dust_xobj);
    }

    auto fuel_xo = root.Get<tiny::xobject>("FUEL");
    for (auto i = 0; i < fuel_xo.Count(); i++) {
      fuel_xo.Enter(i);
      auto coke_xobj = fuel_xo.Get<tiny::xobject>("COKE");
      auto coal_xobj = fuel_xo.Get<tiny::xobject>("COAL");
      coke = Fuel_Element::build_fuel_element(coke_xobj);
      coal = Fuel_Element::build_fuel_element(coal_xobj);
    }

    auto condition_json = get_json(root, "CONDITION");
    condition.diffusion = Diffusion::build_diffusion(
        condition_json.Get<tiny::xobject>("DIFUSSION"));
    condition.hot_metal = get_json(condition_json, "HOT_METAL");
    condition.org = condition_json;
    return *this;
  }
  BF_Info &get_burden_ratio(vector<string> equVariables) {
    map<string, double> balance_material = {};
    // Fe????????????
    auto tFe =
        10.0 * condition.hot_metal.Getd("Fe") +
        10.0 * condition.hot_metal.Getd("Fe") *
            condition.diffusion.SLAG.Getd("Fe") /
            condition.diffusion.HOT_METAL.Getd("Fe") +
        burdens.dust.get_element_volumn("TFe") - 0. -
        burdens.partical_Fe.get_element_volumn("TFe") -
        M_Fe / M_FeO *
            (coke.get_content_volumn("FeO") + coal.get_content_volumn("FeO")) -
        M_Fe / M_FeS *
            (coke.get_content_volumn("FeS") + coal.get_content_volumn("FeS"));

    // P ????????????
    auto tP = 10. * condition.hot_metal.Getd("P") +
              10. * condition.hot_metal.Getd("P") *
                  condition.diffusion.SLAG.Getd("P") /
                  condition.diffusion.HOT_METAL.Getd("P") +
              burdens.dust.get_element_volumn("P") - 0. -
              burdens.partical_Fe.get_element_volumn("P") -
              2. * M_P / M_P2O5 * coke.get_content_volumn("P2O5") -
              2. * M_P / M_P2O5 * coal.get_content_volumn("P2O5");

    // Mn ????????????
    auto tMn = 10. * condition.hot_metal.Getd("Mn") +
               10. * condition.hot_metal.Getd("Mn") *
                   condition.diffusion.SLAG.Getd("Mn") /
                   condition.diffusion.HOT_METAL.Getd("Mn") +
               burdens.dust.get_element_volumn("Mn") - 0. -
               burdens.partical_Fe.get_element_volumn("Mn") -
               M_Mn / M_MnO * coke.get_content_volumn("MnO") -
               M_Mn / M_MnO * coal.get_content_volumn("MnO");

    // S ????????????
    auto tS = 10. * condition.hot_metal.Getd("S") +
              10. * condition.hot_metal.Getd("S") *
                  condition.diffusion.SLAG.Getd("S") /
                  condition.diffusion.GAS.Getd("S") +
              burdens.dust.get_element_volumn("S") - 0. -
              burdens.partical_Fe.get_element_volumn("S") -
              M_S / M_FeS * coke.get_content_volumn("FeS") -
              M_S / M_FeS * coal.get_content_volumn("FeS");

    balance_material.insert({"TFe", tFe});
    balance_material.insert({"P", tP});
    balance_material.insert({"Mn", tMn});
    balance_material.insert({"S", tS});

    // ????????????????????????, ????????????
    {
      vector<double> oreRatio;

      vector<vector<double>> b = {};
      vector<vector<double>> vect = {};
      for (auto x : equVariables) {
        auto _vect = vector<double>();
        auto _b = vector<double>();
        for (auto ore : burdens.ores) {
          _vect.push_back(0.01 * ore.ELEMENT.Getd(x));
          _b.push_back(0);
        }
        _b[0] = balance_material[x];
        vect.push_back(_vect);
        b.push_back(_b);
      }
      cout << "????????????: " << endl;
      printMatrix(vect, 5);
      vector<vector<double>> vect_inv =
          vector<vector<double>>(vect.size(), vect[0]);

      // ???????????????????????????
      get_M_Inverse(vect, vect_inv);

      //????????????
      auto ret = mCross(vect_inv, b);

      cout << "????????????: " << endl;
      printMatrix(ret, 5);

      // ??????????????????
      for (size_t i = 0; i < ret.size(); i++) {
        burdens.ores[i].OMEGA = ret[i][0];
      }
    }

    cout << tP << "\n"
         << tFe << "\n"
         << burdens.ores[0].OMEGA << "\n"
         << burdens.ores[1].OMEGA << "\n"
         << endl;
    return *this;
  }
  BF_Info &get_basic_volumn() {
    // ??????????????????
    auto R0 = condition.org.Getd("SLAG_BASIC"); // ????????????
    auto Q = R0 * (burdens.get_total_input_content("SiO2") +
                   coke.get_content_volumn("SiO2") +
                   coal.get_content_volumn("SiO2") -
                   M_SiO2 / M_Si * 10. * condition.hot_metal.Getd("Si")) -
             (burdens.get_total_input_content("CaO") +
              coke.get_content_volumn("CaO") + coal.get_content_volumn("CaO"));
    Q = Q / (0.01 * (burdens.lime.CONTENT.Getd("CaO") -
                     R0 * burdens.lime.CONTENT.Getd("SiO2")));

    cout << "Q:" << Q << "\n";
    if (Q < 0) {
      Q = 0.0; // ??????????????? ?????????
    }

    burdens.lime.OMEGA = Q;
    cout << "Q:" << Q << "\n";

    return *this;
  }
  BF_Info &get_slag_volumn() {
    vector<string> slag_labels = {};
    vector<double> slag_volumn = {};
    slag_labels =
        vector<string>{"CaO", "SiO2", "MgO", "Al2O3", "(K+Na)2O", "MnO",
                       "S/2", "V2O5", "FeO", "TiO2",  "total"};
    slag_volumn = vector<double>(slag_labels.size(), 0.0);

    // ?????? "CaO", "SiO2", "MgO", "Al2O3","(K+Na)2O"
    for (size_t i = 0; i < 5; i++) {
      slag_volumn[i] = burdens.get_total_input_content(slag_labels[i]) +
                       coke.get_content_volumn(slag_labels[i]) +
                       coal.get_content_volumn(slag_labels[i]);
    }
    // SiO2 ???????????????
    slag_volumn[1] -= M_SiO2 / M_Si * 10. * condition.hot_metal.Getd("Si");
    // MnO
    slag_volumn[5] = M_MnO / M_Mn * 10. * condition.hot_metal.Getd("Mn") *
                     (1 - condition.diffusion.HOT_METAL.Getd("Mn")) /
                     condition.diffusion.HOT_METAL.Getd("Mn");
    // S/2
    auto oreS_2 = burdens.get_ores_content_volumn("S");
    slag_volumn[6] =
        1. / 2. *
        ((1 - condition.diffusion.GAS.Getd("S")) *
             (oreS_2 + burdens.partical_Fe.get_content_volumn("S") +
              burdens.lime.get_content_volumn("S") +
              0.01 * coke.OMEGA * coke.S + 0.01 * coal.OMEGA * coal.S) -
         burdens.dust.get_content_volumn("S") -
         10 * condition.hot_metal.Getd("S"));
    // V2O5
    slag_volumn[7] = M_V2O5 / (2 * M_V) * 10. * condition.hot_metal.Getd("V") *
                     (1 - condition.diffusion.HOT_METAL.Getd("V")) /
                     condition.diffusion.HOT_METAL.Getd("V");

    // FeO
    slag_volumn[8] = M_FeO / M_Fe * 10. * condition.hot_metal.Getd("Fe") *
                     (1 - condition.diffusion.HOT_METAL.Getd("Fe")) /
                     condition.diffusion.HOT_METAL.Getd("Fe");
    // TiO2
    slag_volumn[9] = M_TiO2 / M_Ti * 10. * condition.hot_metal.Getd("Ti") *
                     (1 - condition.diffusion.HOT_METAL.Getd("Ti")) /
                     condition.diffusion.HOT_METAL.Getd("Ti");
    // total
    slag = {};
    for (size_t i = 0; i < slag_labels.size() - 1; i++) {
      slag_volumn[slag_labels.size() - 1] += slag_volumn[i];
      slag.insert({slag_labels[i], slag_volumn[i]});
    }
    slag.insert({slag_labels[slag_labels.size() - 1],
                 slag_volumn[slag_labels.size() - 1]});
    cout << "\n";
    for (auto x : slag) {
      cout << x.first << ": " << x.second << "\t"
           << 100 * x.second / slag_volumn[slag_labels.size() - 1] << "%"
           << endl;
    }

    return *this;
  }
  BF_Info &deSulfur() {
    RO_min = 50. - 0.25 * (slag["Al2O3"] / slag["total"] * 100.) +
             3. * 2. * (slag["S/2"] / slag["total"] * 100.) -
             (0.3 * condition.hot_metal.Getd("Si") +
              30. * condition.hot_metal.Getd("S")) /
                 slag["total"];
    auto RO_elements = array<string, 4>{"CaO", "MgO", "MnO", "FeO"};
    RO = {};
    for (auto x : RO_elements) {
      RO += 100. * slag[x] / slag["total"];
    }
    cout << "RO: " << RO << "   ";
    cout << "RO_min: " << RO_min << "\n";
    return *this;
  }
  BF_Info &get_blast_volumn() {
    auto C_total =
        burdens.get_total_input_content("C") +
        0.01 * coke.OMEGA * (coke.CF + M_C / M_CH4 * coke.CONTENT.Getd("CH4")) +
        0.01 * coal.OMEGA * coal.CF;
    auto C_dFe = M_C / M_Fe * condition.org.Getd("RD") *
                 (10 * condition.hot_metal.Getd("Fe") -
                  burdens.partical_Fe.get_element_volumn("TFe"));
    auto C_dothers = 2 * M_C / M_Si * 10 * condition.hot_metal.Getd("Si") +
                     M_C / M_Mn * 10 * condition.hot_metal.Getd("Mn") +
                     5 * M_C / (2 * M_P) * 10 * condition.hot_metal.Getd("P") +
                     M_C / M_S * slag["S/2"];
    auto C_hot = 10 * condition.hot_metal.Getd("C");
    auto C_CO_2 = 0.5 * M_C / M_CO2 * burdens.lime.OMEGA * 0.6;
    auto oreH2O = burdens.get_ores_content_volumn("H2O");
    auto C_H2O_r = 0.3 * M_C / M_H2O * oreH2O;
    auto C_blast = C_total - C_dFe - C_dothers - C_hot - C_CO_2 - C_H2O_r;
    cout << "C??????" << C_blast << "\n";
    cout << "\nC_total: " << C_total << "\nC_dFe: " << C_dFe
         << "\nC_dothers: " << C_dothers << "\nC_hot: " << C_hot
         << "\nC_CO_2: " << C_CO_2 << "\nC_H2O_r: " << C_H2O_r << "\n";
    auto blastHumi = 0.01 * condition.org.Getd("BLAST_HUMIDITY");
    auto blastOxy = 0.01 * condition.org.Getd("BLAST_O2");
    auto concentOxy = (1 - blastHumi) * blastOxy + 0.5 * blastHumi;
    auto blastVolumn = (IdealR / (M_C * 2) * C_blast -
                        0.01 * coal.OMEGA *
                            (IdealR / (M_O * 2) * coal.CONTENT.Getd("O2") +
                             IdealR / M_H2O * coal.CONTENT.Getd("H2O"))) /
                       concentOxy; // V???
    // ????????????
    auto gammaBlast = 32. / IdealR * (1 - blastHumi) * blastOxy +
                      28. / IdealR * (1 - blastHumi) * (1 - blastOxy) +
                      18. / IdealR * blastHumi;
    // ????????????
    auto weightBlast = blastVolumn * gammaBlast;
    // ????????????
    auto humiVolumnBlast = blastHumi * blastVolumn;
    // ????????????
    auto dryVolumnBlast = (1 - blastHumi) * blastVolumn;
    // ????????????
    auto dryWeightBlast = dryVolumnBlast * (32. / IdealR * blastOxy +
                                            28. / IdealR * (1 - blastOxy));
    // ????????????C
    auto directC = C_dFe + C_dothers + 2 * C_CO_2 + C_H2O_r;
    blast.insert({"C??????", C_blast});
    blast.insert({"??????????????????", concentOxy});
    blast.insert({"V???", blastVolumn});
    blast.insert({"????????????", gammaBlast});
    blast.insert({"????????????", weightBlast});
    blast.insert({"????????????", humiVolumnBlast});
    blast.insert({"????????????", dryVolumnBlast});
    blast.insert({"????????????", dryWeightBlast});
    blast.insert({"????????????C", directC});

    cout << "\n";
    for (auto x : blast) {
      cout << x.first << ": " << x.second << endl;
    }

    return *this;
  }
  BF_Info &get_topgas_volumn() {
    auto blastHumi = 0.01 * condition.org.Getd("BLAST_HUMIDITY");
    auto blastOxy = 0.01 * condition.org.Getd("BLAST_O2");
    auto oreH2O = burdens.get_ores_content_volumn("H2O");
    auto insertH2 = blastHumi * blast["V???"] +
                    0.01 * 11.2 * coke.OMEGA *
                        (coke.CONTENT.Getd("H2") +
                         2 * M_H2 / M_CH4 * coke.CONTENT.Getd("CH4")) +
                    0.01 * 11.2 * coal.OMEGA *
                        (coal.CONTENT.Getd("H2") +
                         M_H2 / M_H2O * coal.CONTENT.Getd("H2O")) +
                    0.3 * IdealR / M_H2O * oreH2O;
    double etaH2 = 0.4;
    while (1) {
      auto reductH2 = insertH2 * etaH2;
      auto V_H2 = (1 - etaH2) * insertH2;

      // ????????????????????????????????? CO2
      auto indirect_CO2_from_oxid =
          IdealR / M_Fe2O3 * burdens.get_total_input_content("Fe2O3") +
          IdealR / M_MnO2 * burdens.get_total_input_content("MnO2");

      // CO ????????????FeO?????? CO2
      auto riH2 = M_Fe / IdealR * reductH2 /
                  (10. * condition.hot_metal.Getd("Fe") -
                   burdens.partical_Fe.get_element_volumn("TFe"));
      auto indirect_CO2_from_Iron =
          IdealR / M_Fe * (1 - condition.org.Getd("RD") - riH2) *
          (10. * condition.hot_metal.Getd("Fe") -

           burdens.partical_Fe.get_element_volumn("TFe"));

      // ????????????????????? CO2
      auto oreCO2 = burdens.get_ores_content_volumn("CO2");
      auto demixore_CO2 =
          IdealR / M_CO2 * (oreCO2 - burdens.dust.get_content_volumn("CO2"));
      // ?????????????????? CO2
      auto lampCO2 =
          IdealR / M_CO2 * (1 - 0.5) * burdens.lime.get_content_volumn("CO2");
      // ?????????????????? CO2
      auto cokeCO2 = IdealR / M_CO2 * coke.get_content_volumn("CO2");
      // ???????????????????????? CO2
      auto V_CO2 = indirect_CO2_from_oxid + indirect_CO2_from_Iron +
                   demixore_CO2 + lampCO2 + cokeCO2;

      // ???????????? CO
      auto V_CO = IdealR / M_C * blast["C??????"] +
                  IdealR / M_C * blast["????????????C"] +
                  IdealR / M_CO * coke.get_content_volumn("CO") -
                  indirect_CO2_from_Iron - indirect_CO2_from_oxid;
      // ??? N2
      auto V_N2 =
          blast["V???"] * (1 - blastHumi) * (1 - blastOxy) +
          IdealR / M_N2 *
              (coke.get_content_volumn("N2") + coal.get_content_volumn("N2"));
      // ?????????
      auto V_Gas = V_CO + V_CO2 + V_N2 + V_H2;
      auto gammaGas = M_CO2 / IdealR * V_CO2 / V_Gas +
                      M_CO / IdealR * (V_CO + V_N2) / V_Gas +
                      2 / IdealR * V_H2 / V_Gas;
      auto weightGas = V_Gas * gammaGas;
      // ?????????
      auto gasH2O = IdealR / M_H2O *
                    (M_H2O / IdealR * reductH2 +
                     (1 - 0.3) * burdens.get_ores_content_volumn("H2O"));

      auto etaCO =
          (indirect_CO2_from_Iron + indirect_CO2_from_oxid) /
          (IdealR / M_C * blast["C??????"] + IdealR / M_C * blast["????????????C"]);

      if (abs(etaCO - etaH2) < 0.001) {
        topgas = {};
        topgas.insert({"CO", V_CO});
        topgas.insert({"CO2", V_CO2});
        topgas.insert({"H2", V_H2});
        topgas.insert({"N2", V_N2});
        topgas.insert({"???????????????(???)", V_Gas});
        topgas.insert({"????????????(???)", weightGas});
        topgas.insert({"???????????????", gasH2O});
        topgas.insert({"CO?????????", etaCO});
        topgas.insert({"H2?????????", etaH2});
        topgas.insert({"CO??????????????????",
                       indirect_CO2_from_Iron + indirect_CO2_from_oxid});
        topgas.insert({"????????????H2", reductH2});
        cout << "\n";
        for (auto x : topgas) {
          cout << x.first << ": " << x.second << endl;
        }
        cout << "etaCO: " << etaCO << "\netaH2: " << etaH2 << "\n";
        break; // ????????????
      } else {
        etaH2 += 0.0005 * (etaCO - etaH2);
      }
    };
    cout << "\n";
    return *this;
  }
  BF_Info &check_material_balance() {
    double oreWeight = {};
    for (auto x : burdens.ores) {
      oreWeight += x.OMEGA;
    }
    auto in = oreWeight + burdens.lime.OMEGA + burdens.partical_Fe.OMEGA +
              coal.OMEGA + coke.OMEGA + blast["????????????"];
    auto out = 1000 + slag["total"] + topgas["????????????(???)"] +
               topgas["???????????????"] + burdens.dust.OMEGA;
    auto absLoss = abs(in - out);
    // auto absLoss = in - out;
    auto relLoss = absLoss / in;
    // cout << "????????????: " << in << "  ????????????: " <<out <<endl;
    cout << "????????????  ????????????: " << absLoss << "  ????????????: " << relLoss * 100
         << "%\n\n";
    return *this;
  }
  BF_Info &check_rd() {
    auto C_gas = 0.01 * coke.CF * coke.OMEGA + 0.01 * coke.CF * coal.OMEGA +
                 burdens.get_total_input_content("C") -
                 10 * condition.hot_metal.Getd("C");
    auto beta = topgas["CO?????????"] / (1 - topgas["CO?????????"]);
    auto riCO = M_Fe /
                (10 * condition.hot_metal.Getd("Fe") -
                 burdens.partical_Fe.get_element_volumn("TFe")) *
                (C_gas / M_C * topgas["CO2"] / (topgas["CO2"] + topgas["CO"]) -
                 0.5 / M_CO2 * burdens.lime.get_content_volumn("CO2") -
                 burdens.get_ores_content_volumn("CO2") / M_CO2 -
                 coke.get_content_volumn("CO2") / M_CO2 -
                 burdens.get_total_input_content("Fe2O3") / M_Fe2O3 -
                 burdens.get_total_input_content("MnO2") / M_MnO2);
    auto riH2 = M_Fe * C_gas /
                (M_C * (10 * condition.hot_metal.Getd("Fe") -
                        burdens.partical_Fe.get_element_volumn("TFe"))) *
                beta * topgas["H2"] / (topgas["CO2"] + topgas["CO"]);
    auto ri = riH2 + riCO;
    cout << "?????? rd: " << 1 - ri << "\n\n";
    return *this;
  }
  BF_Info &get_heat_balance() {
    /**
    ????????????
    */
    auto blastHumi = 0.01 * condition.org.Getd("BLAST_HUMIDITY");
    auto blastOxy = 0.01 * condition.org.Getd("BLAST_O2");
    auto qBlast_C = 9800 * blast["C??????"];
    auto qCd = 9800 * blast["????????????C"];
    auto qiCO = 12650 * topgas["CO??????????????????"];
    auto qiH2 = 10800 * topgas["????????????H2"];
    auto blastC = blastOxy * (1 - blastHumi) * CpO2 +
                  (1 - blastOxy) * (1 - blastHumi) * CpN2 + blastHumi * CpH2O;
    // ?????????????????????
    auto qBlast =
        blast["V???"] *
        (blastC * condition.org.Getd("BLAST_TEMPERATURE") - 10800 * blastHumi);
    // ?????????
    auto qSlag_F = 1130 * (slag["CaO"] + slag["MgO"]);
    // ?????????????????????, ????????????, ???0
    auto qBurden = 0;
    /**
    ????????????
    */
    auto omega_FeO_Fe = 0.2 * burdens.ores[0].get_content_volumn("FeO") +
                        coke.get_content_volumn("FeO") +
                        coal.get_content_volumn("FeO") -
                        1000 * condition.hot_metal.Getd("FeO");
    auto omega_Fe2SiO4_FeO = M_Fe2SiO4 / M_FeO *
                             (burdens.get_total_input_content("FeO") -
                              0.2 * burdens.ores[0].get_content_volumn("FeO"));
    auto omega_Fe2O3_FeO =
        burdens.get_ores_content_volumn("Fe2O3") -
        M_Fe2O3 / M_FeO *
            (burdens.get_total_input_content("FeO") -
             0.2 * burdens.ores[0].get_content_volumn("FeO"));
    // ?????????????????????
    auto qOxide = 4801. * omega_FeO_Fe + 4806. * omega_Fe2SiO4_FeO +
                  5159. * omega_Fe2O3_FeO +
                  25733 * 10 * condition.hot_metal.Getd("P") +
                  7358 * 10 * condition.hot_metal.Getd("Mn") +
                  31059 * 10 * condition.hot_metal.Getd("Si");

    // ????????????
    auto qS = 8300 * slag["S/2"] * 2;
    // ?????????????????????
    double qCarbonic = {};
    for (auto ore : burdens.ores) {
      if (ore.CONTENT.Getd("CO2") >= 1e-6) {
        qCarbonic += // ....? ?????? CO2XXX
            (4040 + 0.5 * 3770) * ore.get_content_volumn("CaO") +
            2307 * ore.get_content_volumn("MgO") +
            1918 * ore.get_content_volumn("FeO") +
            2650 * ore.get_content_volumn("MnO");
      }
    }

    // ????????????????????????
    auto qBlow = coal.OMEGA *
                 (blowHeat + (331 + 13440) * coal.CONTENT.Getd("H2O"));
    // ??????????????????
    auto qSlag_O = slag["total"] * 1780;
    // ?????????????????????
    auto qPFe = 0; // ????????????
    // ???????????????
    auto qIron = 1000 * 1240;

    // ??????????????????
    auto CpGas = CpCO * topgas["CO"] / topgas["???????????????(???)"] +
                 CpCO2 * topgas["CO2"] / topgas["???????????????(???)"] +
                 CpN2 * topgas["N2"] / topgas["???????????????(???)"] +
                 CpH2 * topgas["H2"] / topgas["???????????????(???)"];
    //
    auto tempTop = condition.org.Getd("TOP_GAS_TEMPERATURE");
    auto reductH2O = topgas["????????????H2"] * M_H2O / IdealR;
    auto chemicalH2O = burdens.get_total_input_content("H2O") +
                       coke.get_content_volumn("H2O") +
                       coal.get_content_volumn("H2O");
    double oresH2O = {};
    for (auto x : burdens.ores) {
      oresH2O += x.MOIST;
    }
    auto materialH2O =
        oresH2O + 0.01 * burdens.lime.OMEGA * burdens.lime.MOIST +
        0.01 * burdens.partical_Fe.OMEGA * burdens.partical_Fe.MOIST +
        0.01 * coke.OMEGA * coke.MOIST + 0.01 * coal.OMEGA * coal.MOIST;
    // ?????????????????????
    auto qTopgas =
        (topgas["???????????????(???)"] * CpGas + CpH2O * reductH2O) * tempTop +
        (2450 + 1.244 * CpH2O * tempTop) * materialH2O +
        (331 + 0.7 * (2450 + 1.244 * CpH2O * tempTop) + 0.3 * 6150) *
            chemicalH2O +
        0.8 * burdens.dust.OMEGA * tempTop;
    // ???????????????????????????
    auto qIn = qBlast_C + qCd + qiCO + qiH2 + qBlast + qSlag_F + qBurden;
    auto qOut =
        qOxide + qS + qCarbonic + qBlow + qSlag_O + qPFe + qIron + qTopgas;
    auto qLoss = abs(qIn - qOut);

    auto qEffect = qIn - qTopgas;
    // ??????????????????
    auto eta_t = qEffect / qIn;
    // ??????????????????
    auto eta_C = 0.293 + 0.707 * topgas["CO?????????"];
    total_heat = {};
    total_heat.insert({"???????????????????????????", qBlast_C});
    total_heat.insert({"????????????C?????????CO", qCd});
    total_heat.insert({"??????CO?????????CO2", qiCO});
    total_heat.insert({"????????????H2?????????H2O", qiH2});
    total_heat.insert({"?????????????????????", qBlast});
    total_heat.insert({"?????????", qSlag_F});
    total_heat.insert({"?????????????????????", qBurden});
    total_heat.insert({"?????????????????????", qOxide});
    total_heat.insert({"????????????", qS});
    total_heat.insert({"?????????????????????", qCarbonic});
    total_heat.insert({"????????????????????????", qBlow});
    total_heat.insert({"??????????????????", qPFe});
    total_heat.insert({"???????????????", qSlag_O});
    total_heat.insert({"???????????????", qIron});
    total_heat.insert({"?????????????????????", qTopgas});
    total_heat.insert({"???????????????????????????", qLoss});
    total_heat.insert({"????????????", qIn});
    total_heat.insert({"??????????????????", eta_t});
    total_heat.insert({"??????????????????", eta_C});
    for (auto x : total_heat) {
      cout << x.first << ": " << x.second << endl;
    }
    cout << "???????????????: " << 100 * qLoss / qIn << "%" << endl;
    cout << "\n";
    return *this;
  }
  BF_Info &get_area_heat_balance() {
    /**
    ????????????
    */
    double oresOmega = {};
    for (auto x : burdens.ores) {
      oresOmega += x.OMEGA;
    }
    // ???????????????????????????
    auto qC_blast =
        9800 * blast["C??????"] -
        coal.OMEGA * (coalHeatp + (331 + 13440) * coal.CONTENT.Getd("H2O"));
    // ?????????????????????
    auto qBlast = total_heat["?????????????????????"];
    auto qCoke =
        (coke.OMEGA -
         burdens.dust.OMEGA * burdens.dust.get_content_volumn("C") / coke.CF) *
        1.507 * 950;
    auto q_ore_and_lime =
        ((oresOmega + burdens.lime.OMEGA) -
         burdens.dust.OMEGA *
             (1 - burdens.dust.get_content_volumn("C") / coke.CF) -
         0.5 * burdens.lime.OMEGA * burdens.get_ores_content_volumn("CO2") -
         32. / 44.8 * (topgas["CO??????????????????"] + topgas["????????????H2"])) *
        1.507 * 950;
    // ?????????????????????????????????
    auto qBurden = qCoke + q_ore_and_lime;
    /**
    ????????????
    */
    // ??????????????????
    auto qDirect = 2890 * condition.org.Getd("RD") *
                       (10 * condition.hot_metal.Getd("Fe") -
                        burdens.partical_Fe.get_content_volumn("Fe")) +
                   22960 * 10 * condition.hot_metal.Getd("Si") +
                   4877 * 10 * condition.hot_metal.Getd("Mn") +
                   26520 * 10 * condition.hot_metal.Getd("P") +
                   11310 * 10 * condition.hot_metal.Getd("V") +
                   9500 * 10 * condition.hot_metal.Getd("Ti");
    // ????????????
    auto qs = 4650 * slag["S/2"]; // S
    // ?????????????????????????????????
    auto VGas_i = topgas["???????????????(???)"] -
                  IdealR / 44. *
                      (burdens.get_total_input_content("CO2") -
                       0.5 * burdens.lime.get_content_volumn("CO2") +
                       coke.get_content_volumn("CO2"));
    // ???????????????????????????
    auto q_Gas = VGas_i * 1.411 * 1000 - qCoke;

    area_heat = {};
    area_heat.insert({"???????????????????????????", qC_blast});
    area_heat.insert({"?????????????????????", qBlast});
    area_heat.insert({"?????????????????????????????????", qBurden});
    area_heat.insert({"??????????????????", qDirect});
    area_heat.insert({"????????????", qs});
    area_heat.insert({"???????????????????????????", q_Gas});
    for (auto x : area_heat) {
      cout << x.first << ": " << x.second << endl;
    }
    cout << "\n";
    return *this;
  }
  BF_Info &print_mixed_ore_content(vector<string> element_labels,
                                   vector<string> content_labels) {
    double mixOmega = {};
    double mixMoist = {};
    for (auto ore : burdens.ores) {
      mixOmega += ore.OMEGA;
    }
    for (auto ore : burdens.ores) {
      mixMoist += ore.MOIST * ore.OMEGA / mixOmega;
    }

    cout << "???????????????: " << endl;
    cout << "??????: " << mixOmega << endl;
    cout << "??????: " << endl;
    for (size_t i = 0; i < burdens.ores.size(); i++) {
      cout << "  " << oresLabel[i] << ": "
           << 100 * burdens.ores[i].OMEGA / mixOmega << "%" << endl;
    }
    cout << "????????????: " << endl;
    for (size_t i = 0; i < element_labels.size(); i++) {
      double _element_ratio = {};
      for (auto ore : burdens.ores) {
        _element_ratio +=
            ore.OMEGA / mixOmega * ore.ELEMENT.Getd(element_labels[i]);
      }
      cout << "  " << element_labels[i] << ": " << _element_ratio << "%"
           << endl;
    }
    cout << "????????????: " << endl;
    for (size_t i = 0; i < content_labels.size(); i++) {
      double _content_ratio = {};
      for (auto ore : burdens.ores) {
        _content_ratio +=
            ore.OMEGA / mixOmega * ore.CONTENT.Getd(content_labels[i]);
      }
      cout << "  " << content_labels[i] << ": " << _content_ratio << "%"
           << endl;
    }
    cout << "??????: " << mixMoist << endl;
    cout << "\n";
    return *this;
  }
};
