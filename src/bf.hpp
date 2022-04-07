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
    // Fe平衡方程
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

    // P 平衡方程
    auto tP = 10. * condition.hot_metal.Getd("P") +
              10. * condition.hot_metal.Getd("P") *
                  condition.diffusion.SLAG.Getd("P") /
                  condition.diffusion.HOT_METAL.Getd("P") +
              burdens.dust.get_element_volumn("P") - 0. -
              burdens.partical_Fe.get_element_volumn("P") -
              2. * M_P / M_P2O5 * coke.get_content_volumn("P2O5") -
              2. * M_P / M_P2O5 * coal.get_content_volumn("P2O5");

    // Mn 平衡方程
    auto tMn = 10. * condition.hot_metal.Getd("Mn") +
               10. * condition.hot_metal.Getd("Mn") *
                   condition.diffusion.SLAG.Getd("Mn") /
                   condition.diffusion.HOT_METAL.Getd("Mn") +
               burdens.dust.get_element_volumn("Mn") - 0. -
               burdens.partical_Fe.get_element_volumn("Mn") -
               M_Mn / M_MnO * coke.get_content_volumn("MnO") -
               M_Mn / M_MnO * coal.get_content_volumn("MnO");

    // S 平衡方程
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

    // 矩阵运算线性方程, 求解配矿
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
      cout << "系数矩阵: " << endl;
      printMatrix(vect, 5);
      vector<vector<double>> vect_inv =
          vector<vector<double>>(vect.size(), vect[0]);

      // if (condition.org.Get<bool>("RATIO") == true){
      //   for(auto x:)
      // }
      get_M_Inverse(vect, vect_inv);
      auto ret = mCross(vect_inv, b);
      cout << "结果矩阵: " << endl;
      printMatrix(ret, 5);

      // 保存计算结果
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
    // 碱度平衡方程
    auto R0 = condition.org.Getd("SLAG_BASIC"); // 炉渣碱度
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
      Q = 0.0; // 不需要熔剂 石灰石
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

    // 处理 "CaO", "SiO2", "MgO", "Al2O3","(K+Na)2O"
    for (size_t i = 0; i < 5; i++) {
      slag_volumn[i] = burdens.get_total_input_content(slag_labels[i]) +
                       coke.get_content_volumn(slag_labels[i]) +
                       coal.get_content_volumn(slag_labels[i]);
    }
    // SiO2 另外做处理
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
      cout << x.first << ": " << x.second << endl;
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
    cout << "C风口" << C_blast << "\n";
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
                       concentOxy; // V风
    // 鼓风比重
    auto gammaBlast = 32. / IdealR * (1 - blastHumi) * blastOxy +
                      28. / IdealR * (1 - blastHumi) * (1 - blastOxy) +
                      18. / IdealR * blastHumi;
    // 鼓风质量
    auto weightBlast = blastVolumn * gammaBlast;
    // 风中水分
    auto humiVolumnBlast = blastHumi * blastVolumn;
    // 干风体积
    auto dryVolumnBlast = (1 - blastHumi) * blastVolumn;
    // 干风质量
    auto dryWeightBlast = dryVolumnBlast * (32. / IdealR * blastOxy +
                                            28. / IdealR * (1 - blastOxy));
    // 直接还原C
    auto directC = C_dFe + C_dothers + 2 * C_CO_2 + C_H2O_r;
    blast.insert({"C风口", C_blast});
    blast.insert({"鼓风氧气浓度", concentOxy});
    blast.insert({"V风", blastVolumn});
    blast.insert({"鼓风比重", gammaBlast});
    blast.insert({"鼓风质量", weightBlast});
    blast.insert({"风中水分", humiVolumnBlast});
    blast.insert({"干风体积", dryVolumnBlast});
    blast.insert({"干风质量", dryWeightBlast});
    blast.insert({"直接还原C", directC});

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
    auto insertH2 = blastHumi * blast["V风"] +
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

      // 高价氧化物间接还原得到 CO2
      auto indirect_CO2_from_oxid =
          IdealR / M_Fe2O3 * burdens.get_total_input_content("Fe2O3") +
          IdealR / M_MnO2 * burdens.get_total_input_content("MnO2");

      // CO 间接还原FeO得到 CO2
      auto riH2 = M_Fe / IdealR * reductH2 /
                  (10. * condition.hot_metal.Getd("Fe") -
                   burdens.partical_Fe.get_element_volumn("TFe"));
      auto indirect_CO2_from_Iron =
          IdealR / M_Fe * (1 - condition.org.Getd("RD") - riH2) *
          (10. * condition.hot_metal.Getd("Fe") -

           burdens.partical_Fe.get_element_volumn("TFe"));

      // 混合矿分解得到 CO2
      auto oreCO2 = burdens.get_ores_content_volumn("CO2");
      auto demixore_CO2 =
          IdealR / M_CO2 * (oreCO2 - burdens.dust.get_content_volumn("CO2"));
      // 熔剂分解得到 CO2
      auto lampCO2 =
          IdealR / M_CO2 * (1 - 0.5) * burdens.lime.get_content_volumn("CO2");
      // 焦炭挥发得到 CO2
      auto cokeCO2 = IdealR / M_CO2 * coke.get_content_volumn("CO2");
      // 总共从炉顶出来的 CO2
      auto V_CO2 = indirect_CO2_from_oxid + indirect_CO2_from_Iron +
                   demixore_CO2 + lampCO2 + cokeCO2;

      // 接下来算 CO
      auto V_CO = IdealR / M_C * blast["C风口"] +
                  IdealR / M_C * blast["直接还原C"] +
                  IdealR / M_CO * coke.get_content_volumn("CO") -
                  indirect_CO2_from_Iron - indirect_CO2_from_oxid;
      // 算 N2
      auto V_N2 =
          blast["V风"] * (1 - blastHumi) * (1 - blastOxy) +
          IdealR / M_N2 *
              (coke.get_content_volumn("N2") + coal.get_content_volumn("N2"));
      // 干煤气
      auto V_Gas = V_CO + V_CO2 + V_N2 + V_H2;
      auto gammaGas = M_CO2 / IdealR * V_CO2 / V_Gas +
                      M_CO / IdealR * (V_CO + V_N2) / V_Gas +
                      2 / IdealR * V_H2 / V_Gas;
      auto weightGas = V_Gas * gammaGas;
      // 含水量
      auto gasH2O = IdealR / M_H2O *
                    (M_H2O / IdealR * reductH2 +
                     (1 - 0.3) * burdens.get_ores_content_volumn("H2O"));

      auto etaCO =
          (indirect_CO2_from_Iron + indirect_CO2_from_oxid) /
          (IdealR / M_C * blast["C风口"] + IdealR / M_C * blast["直接还原C"]);

      if (abs(etaCO - etaH2) < 0.001) {
        topgas = {};
        topgas.insert({"CO", V_CO});
        topgas.insert({"CO2", V_CO2});
        topgas.insert({"H2", V_H2});
        topgas.insert({"N2", V_N2});
        topgas.insert({"炉顶煤气量(干)", V_Gas});
        topgas.insert({"煤气质量(干)", weightGas});
        topgas.insert({"煤气含水量", gasH2O});
        topgas.insert({"CO利用率", etaCO});
        topgas.insert({"H2利用率", etaH2});
        topgas.insert({"CO间接还原消耗",
                       indirect_CO2_from_Iron + indirect_CO2_from_oxid});
        topgas.insert({"还原部分H2", reductH2});
        cout << "\n";
        for (auto x : topgas) {
          cout << x.first << ": " << x.second << endl;
        }
        cout << "etaCO: " << etaCO << "\netaH2: " << etaH2 << "\n";
        break; // 迭代结束
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
              coal.OMEGA + coke.OMEGA + blast["鼓风质量"];
    auto out = 1000 + slag["total"] + topgas["煤气质量(干)"] +
               topgas["煤气含水量"] + burdens.dust.OMEGA;
    auto absLoss = abs(in - out);
    // auto absLoss = in - out;
    auto relLoss = absLoss / in;
    // cout << "收入质量: " << in << "  支出质量: " <<out <<endl;
    cout << "物料平衡  绝对误差: " << absLoss << "  相对误差: " << relLoss * 100
         << "%\n\n";
    return *this;
  }
  BF_Info &check_rd() {
    auto C_gas = 0.01 * coke.CF * coke.OMEGA + 0.01 * coke.CF * coal.OMEGA +
                 burdens.get_total_input_content("C") -
                 10 * condition.hot_metal.Getd("C");
    auto beta = topgas["CO利用率"] / (1 - topgas["CO利用率"]);
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
    cout << "校检 rd: " << 1 - ri << "\n\n";
    return *this;
  }
  BF_Info &get_heat_balance() {
    /**
    热收入项
    */
    auto blastHumi = 0.01 * condition.org.Getd("BLAST_HUMIDITY");
    auto blastOxy = 0.01 * condition.org.Getd("BLAST_O2");
    auto qBlast_C = 9800 * blast["C风口"];
    auto qCd = 9800 * blast["直接还原C"];
    auto qiCO = 12650 * topgas["CO间接还原消耗"];
    auto qiH2 = 10800 * topgas["还原部分H2"];
    auto blastC = blastOxy * (1 - blastHumi) * CpO2 +
                  (1 - blastOxy) * (1 - blastHumi) * CpN2 + blastHumi * CpH2O;
    // 鼓风带入物理热
    auto qBlast =
        blast["V风"] *
        (blastC * condition.org.Getd("BLAST_TEMPERATURE") - 10800 * blastHumi);
    // 成渣热
    auto qSlag_F = 1130 * (slag["CaO"] + slag["MgO"]);
    // 炉料带入物理热, 冷矿入炉, 为0
    auto qBurden = 0;
    /**
    热支出项
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
    // 氧化物分解耗热
    auto qOxide = 4801. * omega_FeO_Fe + 4806. * omega_Fe2SiO4_FeO +
                  5159. * omega_Fe2O3_FeO +
                  25733 * 10 * condition.hot_metal.Getd("P") +
                  7358 * 10 * condition.hot_metal.Getd("Mn") +
                  31059 * 10 * condition.hot_metal.Getd("Si");

    // 脱硫耗热
    auto qS = 8300 * slag["S/2"] * 2;
    // 碳酸盐分解耗热
    auto qCarbonic = // ....? 啥是 CO2XXX
        (4040 + 0.5 * 3770) * burdens.get_ores_content_volumn("CaO") +
        2307 * burdens.get_ores_content_volumn("MgO") +
        1918 * burdens.get_ores_content_volumn("FeO") +
        2650 * burdens.get_ores_content_volumn("MnO");
    // 喷吹燃料分解耗热
    auto qBlow = coal.OMEGA *
                 (blowHeat + (331 + 13440) * coal.get_content_volumn("H2O"));
    // 炉渣带走热量
    auto qSlag_O = slag["total"] * 1780;
    // 废铁升温熔化热
    auto qPFe = 0; // 太少不算
    // 铁水带走热
    auto qIron = 1000 * 1240;

    // 总的体积比热
    auto CpGas = CpCO * topgas["CO"] / topgas["炉顶煤气量(干)"] +
                 CpCO2 * topgas["CO2"] / topgas["炉顶煤气量(干)"] +
                 CpN2 * topgas["N2"] / topgas["炉顶煤气量(干)"] +
                 CpH2 * topgas["H2"] / topgas["炉顶煤气量(干)"];
    //
    auto tempTop = condition.org.Getd("TOP_GAS_TEMPERATURE");
    auto reductH2O = topgas["还原部分H2"] * M_H2O / IdealR;
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
    // 炉顶煤气带走热
    auto qTopgas =
        (topgas["炉顶煤气量(干)"] * CpGas + CpH2O * reductH2O) * tempTop +
        (2450 + 1.244 * CpH2O * tempTop) * materialH2O +
        (331 + 0.7 * (2450 + 1.244 * CpH2O * tempTop) + 0.3 * 6150) *
            chemicalH2O +
        0.8 * burdens.dust.OMEGA * tempTop;
    // 冷却带走和散热损失
    auto qIn = qBlast_C + qCd + qiCO + qiH2 + qBlast + qSlag_F + qBurden;
    auto qOut =
        qOxide + qS + qCarbonic + qBlow + qSlag_O + qPFe + qIron + qTopgas;
    auto qLoss = abs(qIn - qOut);

    auto qEffect = qIn - qTopgas;
    // 热能利用系数
    auto eta_t = qEffect / qIn;
    // 碳素利用系数
    auto eta_C = 0.293 + 0.707 * topgas["CO利用率"];
    total_heat = {};
    total_heat.insert({"风口前碳素燃烧放热", qBlast_C});
    total_heat.insert({"直接还原C氧化成CO", qCd});
    total_heat.insert({"间接CO氧化成CO2", qiCO});
    total_heat.insert({"间接还原H2氧化成H2O", qiH2});
    total_heat.insert({"热风带入物理热", qBlast});
    total_heat.insert({"成渣热", qSlag_F});
    total_heat.insert({"炉料带入物理热", qBurden});
    total_heat.insert({"氧化物分解耗热", qOxide});
    total_heat.insert({"脱硫耗热", qS});
    total_heat.insert({"碳酸盐分解耗热", qCarbonic});
    total_heat.insert({"喷吹燃料分解耗热", qBlow});
    total_heat.insert({"废铁升温耗热", qPFe});
    total_heat.insert({"炉渣带走热", qSlag_O});
    total_heat.insert({"铁水带走热", qIron});
    total_heat.insert({"炉顶煤气带走热", qTopgas});
    total_heat.insert({"冷却带走和散热损失", qLoss});
    total_heat.insert({"总热收入", qIn});
    total_heat.insert({"热能利用系数", eta_t});
    total_heat.insert({"碳素利用系数", eta_C});
    for (auto x : total_heat) {
      cout << x.first << ": " << x.second << endl;
    }
    cout << "\n";
    return *this;
  }
  BF_Info &get_area_heat_balance() {
    /**
    热收入项
    */
    double oresOmega = {};
    for (auto x : burdens.ores) {
      oresOmega += x.OMEGA;
    }
    // 风口前碳素燃烧放热
    auto qC_blast =
        9800 * blast["C风口"] -
        coal.OMEGA *
            (coalHeatp + (331 + 13440) * coal.CONTENT.Getd("H2O"));
    // 热风带入物理热
    auto qBlast = total_heat["热风带入物理热"];
    auto qCoke =
        (coke.OMEGA -
         burdens.dust.OMEGA * burdens.dust.get_content_volumn("C") / coke.CF) *
        1.507 * 950;
    auto q_ore_and_lime =
        ((oresOmega + burdens.lime.OMEGA) -
         burdens.dust.OMEGA *
             (1 - burdens.dust.get_content_volumn("C") / coke.CF) -
         0.5 * burdens.lime.OMEGA * burdens.get_ores_content_volumn("CO2") -
         32. / 44.8 * (topgas["CO间接还原消耗"] + topgas["还原部分H2"])) *
        1.507 * 950;
    // 中温区炉料带入的物理热
    auto qBurden = qCoke + q_ore_and_lime;
    /**
    热支出项
    */
    // 直接还原耗热
    auto qDirect = 2890 * condition.org.Getd("RD") *
                       (10 * condition.hot_metal.Getd("Fe") -
                        burdens.partical_Fe.get_content_volumn("Fe")) +
                   22960 * 10 * condition.hot_metal.Getd("Si") +
                   4877 * 10 * condition.hot_metal.Getd("Mn") +
                   26520 * 10 * condition.hot_metal.Getd("P") +
                   11310 * 10 * condition.hot_metal.Getd("V") +
                   9500 * 10 * condition.hot_metal.Getd("Ti");
    // 脱硫耗热
    auto qs = 4650 * slag["S/2"]; // S
    // 进入间接还原区的煤气量
    auto VGas_i = topgas["炉顶煤气量(干)"] -
                  IdealR / 44. *
                      (burdens.get_total_input_content("CO2") -
                       0.5 * burdens.lime.get_content_volumn("CO2") +
                       coke.get_content_volumn("CO2"));
    // 煤气从高温区带走热
    auto q_Gas = VGas_i * 1.411 * 1000 - qCoke;

    // 焦比
    auto cokeRate =
        (10 * condition.hot_metal.Getd("C") + blast["直接还原C"] +
         blast["C风口"]) /
        burdens.dust.get_content_volumn("C");

    area_heat = {};
    area_heat.insert({"风口前碳素燃烧放热", qC_blast});
    area_heat.insert({"热风带入物理热", qBlast});
    area_heat.insert({"中温区炉料带入的物理热", qBurden});
    area_heat.insert({"直接还原耗热", qDirect});
    area_heat.insert({"脱硫耗热", qs});
    area_heat.insert({"煤气从高温区带走热", q_Gas});
    area_heat.insert({"焦比", cokeRate});
        for (auto x : area_heat) {
      cout << x.first << ": " << x.second << endl;
    }
    cout << "\n";
    return *this;
  }
};
