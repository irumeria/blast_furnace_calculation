
#include "src/bf.hpp"
#include <iostream>
#include <string>

using std::string;

int main() {

  // 第一张配料单所使用的输入参数
  // string filename("input1.json"); //读取输入文件 -- 第一张
  // vector<string> ores_name = {"sinter", "aust", "pellet"}; // 使用的矿名
  // vector<string> equ_elements = {"TFe", "Mn", "P"}; // 配比计算使用的元素

  //  第二张配料单所使用的输入参数
  string filename("input2.json"); //读取输入文件 -- 第二张
  vector<string> ores_name = {"sinter", "aust"}; // 使用的矿名
  vector<string> equ_elements = {"TFe", "P"};   // 配比计算使用的元素

  string file_contents = readFileIntoString(filename);
  tiny::TinyJson root;
  root.ReadJson(file_contents);

  BF_Info bf = {};

  bf.processData(root, ores_name)     // 处理数据
      .get_burden_ratio(equ_elements) // 计算炉料配比以及混合矿成分
      .get_basic_volumn()             // 计算熔剂量
      .get_slag_volumn()              // 计算渣量
      .deSulfur()                     // 核算脱硫能力
      .get_blast_volumn()             // 计算鼓风参数
      .get_topgas_volumn()            // 计算炉顶煤气量
      .check_material_balance()       // 检查物料平衡
      .check_rd()                     // 检查直接还原率
      .get_heat_balance()             // 计算第一总热平衡
      .get_area_heat_balance();       // 计算区域热平衡
}