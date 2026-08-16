# CHESSv2

CHESSv2 是一个纯 Mathematica 的统一谱传播程序包。它把原 canonical
CHESS 与 `CHESS_non_canonical` 的实现放在同一个加载入口下，并由唯一的
公开函数 `SpectralPropagate` 自动选择算法。

这一版刻意保持范围很小：不包含 C++、FORM、FLINT、Fermat 或 Fermatica
后端，也没有后端注册系统。canonical 与 non-canonical 的数值内核被尽量
原样复用；新增代码只负责输入适配、自动分派以及只含 `B0` 时的伪 delta
传播。

## 加载

```wl
Get["/home/liuyuanche/MMA_Package/CHESSv2/CHESSv2.wl"];
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
Core/FakeDelta.wl          B0 伪 delta 传播及自动阶数估计
Core/Dispatch.wl           统一 SpectralPropagate 与结构化分派
Tests/run_all.wls          fresh-kernel 回归测试
```

新增模块中保留了较密集的设计注释，特别标出数学约定、索引关系、失败边界
以及哪些步骤只是经验判据。这样后续若替换某个模块，不必重新猜测其他内核的
隐含约定。

## 测试

```bash
wolframscript -file /home/liuyuanche/MMA_Package/CHESSv2/Tests/run_all.wls
```

测试覆盖裸 `B1`、显式 `B1`、稀疏常数 `B1`、`B0` 自动/固定阶、非对易
`B0(t)`、多边界列、分段传播、三种 evaluator 调用约定、混合
non-canonical 路由、后段失败传播、算符尺寸检查、`RuleDelayed` 分段选项、
端点失败关闭行为以及原始源码哈希未变。
成功时最后打印 `CHESSV2_ALL_TESTS_PASSED`。
