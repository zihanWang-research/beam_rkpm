# Beam / RKPM project handoff

This file is a portable checkpoint for continuing the project on another computer.

## Workspace

- Current workspace: `D:/python/beam`
- Main research roadmap: `projects/README.md`
- MATLAB code map: `code/matlab/README.md`
- Python code map: `code/python/README.md`
- Weekly records: `weekly/`
- Literature and notes: `literature/`

## Research status

### Completed or validated

- PCA–RKPM local frame convergence for circular, conical, and variable helices.
- Quarter-ring learning model and deformation/analytical comparison.
- Raknes 4.1.1, 4.1.2, and 4.1.3 linear nodal-integration benchmarks.
- Two- and three-point Gauss-integration experiments.
- Nodal-integration membrane stabilization and SCNI experiments.
- Three-dimensional rank diagnostics for one-direction versus two-direction bending.
- Workspace reorganization into `projects/`, `code/`, `literature/`, `weekly/`, and `archive/`.

### Current numerical checkpoint

The W38 nodal-integration study contains:

- baseline nodal integration at `alpha = 2.6`;
- membrane-stabilized nodal integration at `alpha = 2.6`;
- a fine alpha scan with a narrow candidate near `alpha = 2.61`;
- convergence CSV files and publication-style figures under
  `weekly/2026/2026-W38_09-14_to_09-20/`.

The observed fine-grid slopes are approximately:

| case | fine-grid slope |
|---|---:|
| nodal baseline, alpha = 2.6 | 2.74 |
| membrane stabilized, alpha = 2.6 | 2.75 |
| nodal baseline, alpha = 2.61 | 3.08 |
| membrane stabilized, alpha = 2.61 | 3.15 |

The alpha = 2.61 result is sensitive and should not be interpreted as a stabilization-only improvement; the baseline also improves at that support radius.

## Current theory-learning checkpoint

The nonlinear study has restarted from the membrane kinematics. Already reviewed:

\[
\mathbf x=\mathbf X+\mathbf u,
\qquad
\mathbf A=\mathbf X_{,\eta},
\qquad
\mathbf a=\mathbf x_{,\eta},
\]

\[
G=\mathbf A^T\mathbf A,
\qquad
g=\mathbf a^T\mathbf a,
\qquad
\mathcal E=\frac{g-G}{2G}.
\]

The next learning step is to derive, in order:

1. \(\delta\mathcal E\);
2. the nonlinear membrane matrix \(\mathbf B_m\);
3. the internal force vector;
4. material and geometric tangent stiffness;
5. Newton–Raphson iteration.

The working rule is: explain the mathematics first, then write code step by step, then check the user's implementation.

## Suggested first message on the new computer

> 我已经打开 beam 工作区，请读取 PROJECT_HANDOFF_2026-09-18.md。我们继续从二维非线性 membrane 的 \(\delta\mathcal E\) 推导开始；我自己写代码，你负责解释公式和检查。

## Environment note

Do not copy `.venv` as the primary environment. Install Python, MATLAB, and required packages on the new computer, then open the copied workspace. The absolute path may change, so replace hard-coded `D:/python/beam` paths in scripts when necessary.
