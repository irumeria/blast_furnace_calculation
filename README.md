# blast_furnace_calculation

### 简介
北科大 "钢铁冶金学" 课程的高炉大计算程序

使用c++编写, 计算过程为设计高炉

输入参数参见 "钢铁冶金学.pdf"

其中, "原燃料及炉尘成分" 版块对应了输入的json文件的 "BURDEN" 部分

"燃料" 对应了输入的json文件的 "FUEL" 部分

"冶炼条件" 对应了输入的json文件的 "CONDITION" 部分

### 运行

+ 需要支持C++ 11 以上的编译器

GNU
```bash
g++ ./main.cpp -o ./main
./main
```

其他如VS / devC++
```
在本文件夹内编译运行main.cpp即可
```

### 示例

内置了前两张配料单的内容作为示例, 物料平衡计算误差均小于0.3%

```c
// in main.cpp

  // 第一张配料单所使用的输入参数
  // string filename("input1.json"); //读取输入文件 -- 第一张
  // vector<string> ores_name = {"sinter", "aust", "pellet"}; // 使用的矿名
  // vector<string> equ_elements = {"TFe", "Mn", "P"}; // 配比计算使用的元素

  //  第二张配料单所使用的输入参数
  string filename("input2.json"); //读取输入文件 -- 第二张
  vector<string> ores_name = {"sinter", "aust"}; // 使用的矿名
  vector<string> equ_elements = {"TFe", "P"};   // 配比计算使用的元素
```

通过变更注释内容的方法切换即可