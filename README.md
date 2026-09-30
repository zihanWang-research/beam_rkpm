# Beam / RKPM 研究工作区

本目录用于研究基于 RKPM 的曲线梁、缆索与空间杆。唯一研究入口是
[`projects/`](projects/README.md)。

## 现在从哪里开始

当前阶段是：

```text
02_nonlinear_static_2d
```

目标是以 Raknes 连续理论为基础，在已经验证的二维线性代码上实现二维非线性静力。

## 目录职责

- `projects/`：01–07 研究主线、阶段目标和主代码入口。
- `code/matlab/`：MATLAB 几何、基准算例、求解器和学习代码。
- `code/python/`：Python 原型、测试和早期分步骤学习版本。
- `literature/docs/`：理论笔记和推导材料。
- `literature/references/papers/`：参考论文和导师资料。
- `weekly/`：每周汇报快照；不在这里继续开发主代码。
- `archive/`：已经被替代的历史代码与结果。

根目录中的 `.venv/`、`.codex/`、`.agents/` 和 `.git/` 属于运行环境或工具配置，不放研究文件。

## 文件管理规则

1. 每个研究阶段只在 `projects/` 中登记一个主代码入口。
2. `weekly/` 中保存当周稳定快照，避免与主代码双向修改。
3. 新结果放在对应阶段的 `results/`，图片放在 `figures/`；阶段真正启动时再建立。
4. 不再为一次参数扫描复制完整主程序，优先写独立的扫描脚本调用主求解器。
5. 原始论文、已验证结果和本人学习代码不删除。

## MATLAB 目录

```text
code/matlab/
├── benchmarks/    # 四分之一圆环和 Raknes 4.1 线性基准
├── geometry/      # PCA、RKPM 参数化和空间曲线标架
├── learning/      # 本人逐步改写的学习版本
└── solvers/       # 可复用的参考求解器
```

## 当前研究位置

当前学习笔记位于：

`projects/02_nonlinear_static_2d/01_configuration_and_membrane_strain.md`
