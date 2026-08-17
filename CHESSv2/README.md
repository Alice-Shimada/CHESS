# CHESSv2

CHESSv2 是统一的谱传播程序包。它把原 canonical CHESS 与
`CHESS_non_canonical` 的实现放在同一个加载入口下，并由唯一的公开函数
`SpectralPropagate` 自动选择算法。默认值仍是纯 Mathematica；规则路径上的
canonical 与只含 `B0` 的方程可以显式切换到可选 C++ 高精度后端。

改动保持模块化：历史 canonical/non-canonical 数值内核继续原样复用，C++
传播、FORM/FLINT 多点表达式求值以及公开分派是互不依赖的模块。没有编译
native 后端时，默认 Mathematica 行为和返回结构不变。

## 加载

```wl
Get["/path/to/CHESS/CHESSv2/CHESSv2.wl"];
```

加载后直接使用 `SpectralPropagate`，不需要写 `CHESS`` 或 `CHESSv2``
前缀。当前版本号保存在 `$CHESSv2Version`。

## 统一方程约定

程序求解

```text
dY(t,eps)/dt = (B0(t) + eps B1(t) + eps^2 B2(t) + ...) Y(t,eps).
```

显式系数列表的第 `p+1` 项表示 `Bp`：

```wl
SpectralPropagate[{B0, B1, B2, ...}, boundary, {t0, t1}, options]
```

自动分派规则是：

| 输入的非零系数 | 采用的算法 |
|---|---|
| 只有 `B0` | 伪 delta 方法，自动或固定截断阶数 |
| 只有 `B1` | 原 canonical 顺序 epsilon 算法 |
| 混合 `B0,B1,...` 或只有更高次项 | 原 non-canonical 算法 |

这里的“零”只作结构判断：字面上的 `0` 或精确全零常数矩阵会被认作零；
程序不会在若干数值点采样一个函数后，凭容差猜测它是否恒为零。

为了兼容旧 canonical 笔记本，裸矩阵求值器仍解释为 `B1`：

```wl
SpectralPropagate[B1, boundary, {t0, t1}, options]
```

一个二维数值矩阵也按常数 `B1` 处理；若要给出多个常数系数，应增加列表
层级，例如 `{{{b0}}, {{b1}}}` 是一维系统的 `{B0,B1}`。

## 矩阵求值器

显式系数列表中的每个 `Bp` 可以是常数数值矩阵，也可以采用以下任一调用
约定：

```wl
Bp[t_, "Precision" -> p_] := ...
Bp[t_, p_] := ...
Bp[t_] := ...
```

返回值必须是数值方阵，而且尺寸必须与 `boundary` 的行数一致。适配层只负责
统一这三种已有约定，不进行符号化简，也不会通过低精度试探改变路由。

若大矩阵来自 FORM 生成的巨大有理表达式，可用
`CHESSNodeEvaluator[scalar,batch]` 提供一次计算所有节点的 batch evaluator。
`Core/FLINTBatchEvaluation.wl` 提供通用的
`CHESSFLINTRunSLP[executable,slp,pointRows,precision,threads]`：它复用实验性
non-canonical 后端的 `NFLOAT01` 协议，把 FORM 的 straight-line program 交给
FLINT `nfloat` 多线程求值。具体变量、SLP 和稀疏矩阵装配仍属于各物理问题，
不会硬编码进 CHESS。完整接口示例见 `FLINT/README.md`。
默认额外 20 位只是可调的经验 guard digits，不是严格误差界；近极点、严重
消去或巨大中间量必须提高 `"GuardDigits"`/`"WorkingBits"`，并用更高精度
重算检查稳定性。任何 FLINT infinity 或 NaN 记录都会失败关闭。

若 batch evaluator 一次返回完整的 `{B0,B1,...}`，通过显式正整数
`"EpsilonDegree"` 让统一入口选择 non-canonical 路线，同时仍只调用一次
多点评值：

```wl
SpectralPropagate[
  CHESSNodeEvaluator[scalarFallback, flintBatch],
  boundary,
  {0, 1},
  "EpsilonDegree" -> 2
]
```

## 可选 C++ 数值传播

编译：

```bash
make -C CHESSv2/Native/Propagation
```

默认 Makefile 使用 Mathematica 14.3 的 WSTP、Boost 1.83、MPFR/MPC/GMP 和
OpenMP；安装位置不同时可覆盖 `BOOST_INCLUDE` 与 `WSTP_DIR`。传播线性代数
使用 MPFR/MPC，而不是 FLINT；FLINT 专用于上一节那类巨大 FORM 表达式的
多点评值。两者如此拆开，是为了让每种库只负责其有优势的热区。

用户仍只调用统一入口：

```wl
rules = {
  "Nodes" -> 32,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "ParallelEvaluation" -> False
};

(* 准备一次，重复传播时复用矩阵采样与 LU 分解。 *)
handle = CHESSNativePrepare[b0, Length[boundary], {0, 1}, rules];

result = SpectralPropagate[
  {b0}, boundary, {0, 1},
  Sequence @@ rules,
  "NumericalBackend" -> "Native",
  "NativeHandle" -> handle,
  "ResultData" -> "Endpoint"
];
```

句柄绑定到准备时的 evaluator 表达式、区间、维数、节点数和精度，并递归
记录 evaluator 依赖的 Wolfram 符号定义；普通函数、辅助函数或选项定义变化
都会使旧句柄自动失效。只有不体现在 Wolfram 定义中的外部可变状态发生变化
时，用户才需要显式重新运行 `CHESSNativePrepare`。

native 选项如下：

| 选项 | 默认值 | 作用 |
|---|---|---|
| `"NumericalBackend"` | `"Mathematica"` | `"Native"` 显式启用；默认不改变历史行为 |
| `"NativeHandle"` | `Automatic` | 复用 `CHESSNativePrepare` 的缓存；`Automatic` 每次自行准备 |
| `"ResultData"` | `"Full"` | `"Endpoint"` 不缓存/回传节点值，适合重复边界传播 |

当前 C++ 模块支持规则 canonical 与规则 `B0` 伪 delta 路由。多段 `B0` 会逐段
准备；奇异端点和真正混合的 `{B0,B1,...}` 仍由已验证的 Mathematica 内核
处理。请求 native 执行尚未支持的混合路线会明确失败，不会悄悄换算法。

## 三类调用示例

以下例子都使用显式系数列表，因此方程含义最清楚。

```wl
ClearAll[b0, b1];
b0[t_, p_] := N[{{1/3}}, p];
b1[t_, p_] := N[{{2/5}}, p];

(* 只有 B0：自动进入伪 delta 路由。 *)
r0 = SpectralPropagate[
  {b0}, N[{{1}}, 70], {0, 1},
  "Nodes" -> 24,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90,
  "DeltaTolerance" -> 10^-40
];

(* 只有 B1：进入 canonical 路由。boundary 的列是 epsilon 系数。 *)
r1 = SpectralPropagate[
  {0, b1}, N[{{1, 0, 0, 0}}, 70], {0, 1},
  "Nodes" -> 24,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90
];

(* B0 与 B1 同时存在：进入原 non-canonical 路由。 *)
r01 = SpectralPropagate[
  {b0, b1}, N[{{1, 0, 0, 0}}, 70], {0, 1},
  "Nodes" -> 24,
  "Precision" -> 70,
  "WorkingPrecisionA" -> 90
];
```

## 只含 B0 时的伪 delta 方法

对

```text
Y'(t) = B0(t) Y(t)
```

引入辅助参数 `delta`：

```text
Y'(t,delta) = delta B0(t) Y(t,delta),
Y(t,delta) = C0(t) + delta C1(t) + delta^2 C2(t) + ... .
```

于是

```text
C0'(t) = 0,
Ck'(t) = B0(t) C(k-1)(t),  k >= 1,
```

恰好成为 canonical CHESS 能直接处理的递推。传播完成后取 `delta=1`，即
求和 `C0+...+CK`。若 `boundary` 有多列，每一列被独立传播，输出仍保持
原来的矩阵形状，不构造巨大的 Kronecker 耦合系统。

自动阶数使用系数与当前部分和定义尺度无关指标

```text
a_k = max_i |C_k(i)| / (1 + |S_k(i)|),   S_k = C0 + ... + Ck.
```

取最近三个系数比的上包络 `q_k`，当 `0 <= q_k < 1` 时估计余项

```text
E_k = 2 a_k q_k / (1-q_k).
```

第一次低于目标容差后，还必须让紧接着的一层也通过，才接受截断阶数。
求解阶数按 `8,16,32,...` 增长，直到收敛或达到 `"MaxDeltaOrder"`。
达到上限仍未通过时返回 `$Failed`，不会把未收敛的部分和当成答案。

需要强调：这个比例余项是经过数值实验检验的经验停止准则，不是对所有未来
系数的严格数学上界。重要计算仍应改变 `"DeltaTolerance"`、提高
`"MaxDeltaOrder"`，并与更高固定阶结果交叉检查。

伪 delta 专用选项：

| 选项 | 默认值 | 作用 |
|---|---:|---|
| `"DeltaOrder"` | `Automatic` | 自动选择；非负整数表示固定阶数 |
| `"DeltaTolerance"` | `Automatic` | 自动值约为 `10^(-0.7 Precision)`，至少 `10^-10` |
| `"MaxDeltaOrder"` | `128` | 自动搜索允许的最高阶 |
| `"Segments"` | `1` | 将规则路径分段传播，每段独立估阶 |

## 返回值

沿用原 CHESS 的六部分结果：

```wl
{nodes, nodeValues, finalValue, endpointBlocks, endpointInfo, methodMetadata}
```

`finalValue` 是最常用的终点值，即结果的第 3 部分。第 6 部分记录实际路由。
伪 delta 自动路由还记录 `SelectedOrder`、`ComputedOrder`、`TailEstimate`、
`Tolerance` 和保护层数，便于复核停止条件。

## 当前边界与失败策略

- `B0` 伪 delta 路由目前只支持规则端点。若提供
  `"RegularizedEndpoints"` 或 `"EndpointData"`，会明确返回 `$Failed`，
  不会静默忽略。
- 只有 `B1` 时继续使用原 canonical 端点正则化实现。
- 一般 non-canonical 情形继续使用原包已有的端点数据与耦合算法。
- 一般 non-canonical 路由若给出 `"EndpointData"`，必须同时指定实际的
  `"RegularizedEndpoints"`；孤立端点数据会失败关闭。
- 一般 non-canonical 路由不能在 `"Segments" > 1` 时重复施加奇异端点
  正则化；这种组合会失败关闭，规则子路径应由调用者明确拆分。
- canonical 分段只支持规则路径；不能把奇异端点选项机械地重复施加到每段。
- 传入非法分段数、非法 delta 阶数或自动阶数达到上限时均失败关闭。
- 算符尺寸错误以及 canonical 内核发出的逐层求解失败都会由统一接口提升为
  整体 `$Failed`，不会返回只填充了一部分的成功形结果。

## 文件结构

```text
CHESSv2.wl                 唯一加载器
Core/Canonical.wl          原 canonical 数值内核（逐字复用）
Core/NonCanonical.wl       原 non-canonical 数值内核（加载路径及失败传播修补）
Core/MatrixAdapters.wl     三种矩阵调用约定和常数矩阵适配
Core/FLINTBatchEvaluation.wl  FORM/FLINT nfloat 多点评值协议
Core/NativeBackend.wl      C++ 传播接口、缓存句柄和结果装配
Core/FakeDelta.wl          B0 伪 delta 传播及自动阶数估计
Core/Dispatch.wl           统一 SpectralPropagate 与结构化分派
Native/Propagation/        MPFR/MPC/OpenMP 传播源代码与 Makefile
FLINT/                     FLINT evaluator 接口说明
Tests/run_all.wls          fresh-kernel 回归测试
Tests/test_native_backend.wls  编译后端专项测试
Tests/test_flint_evaluation.wls  FLINT 二进制协议及可选真实 SLP 测试
Benchmarks/pbb64_native_vs_mathematica.wls  316 维、64 位传播基准
Benchmarks/f3_zoia_flint_benchmark.wls  131 维巨大符号 F3/FLINT 基准
```

新增模块中保留了较密集的设计注释，特别标出数学约定、索引关系、失败边界
以及哪些步骤只是经验判据。这样后续若替换某个模块，不必重新猜测其他内核的
隐含约定。

## 测试

```bash
wolframscript -file CHESSv2/Tests/run_all.wls
```

测试覆盖裸 `B1`、显式 `B1`、稀疏常数 `B1`、`B0` 自动/固定阶、非对易
`B0(t)`、多边界列、分段传播、三种 evaluator 调用约定、混合
non-canonical 路由、后段失败传播、算符尺寸检查、`RuleDelayed` 分段选项、
端点失败关闭行为以及原始源码哈希未变。
成功时最后打印 `CHESSV2_ALL_TESTS_PASSED`。

编译 native 后端后再运行：

```bash
wolframscript -file CHESSv2/Tests/test_native_backend.wls
wolframscript -file CHESSv2/Tests/test_flint_evaluation.wls
```

316 维 PBB 基准严格把数据加载和 native 准备排除在传播计时之外，同时单独
报告准备成本：

```bash
OMP_NUM_THREADS=8 \
  wolframscript -file CHESSv2/Benchmarks/pbb64_native_vs_mathematica.wls
```

Simone Zoia 的 131 维 F3、96 节点、100 位巨大符号方程复现见
`Benchmarks/F3_ZOIA_FLINT_BENCHMARK.md`。在本机冻结源码的 fresh run 中，
纯 Mathematica 节点评值墙钟为 435.87 秒；三次 FLINT 完整节点后端为
36.15--41.11 秒（中位数 39.49 秒），统一入口端到端为 45.86--50.77 秒
（中位数 49.21 秒）。
