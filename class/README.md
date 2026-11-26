# class

分类算法（占位目录）。

和 `geant4/`、`gasdigital/`、`track fit/`、`vertex fit/`、`picture draw/` 平级，
以后把分类器（事件分类、粒子鉴别之类）放进来。

按仓库约定：**代码只放在子目录里**，这个文件的旁边只放说明。

建议的布局（还没有代码，先写在这里）：

```
class/
  app/            命令行入口（读数据 -> 调分类器 -> 写结果）
  class/          Python 包：分类器实现（或者 include/ + src/ 的 C++ 版本）
  README.md       本文件
```
