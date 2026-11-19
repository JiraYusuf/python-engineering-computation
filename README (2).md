# Newton–Raphson Load Flow Analysis Tool (Python)

## Overview
A reusable Python tool for solving steady-state AC power-flow problems in
electrical networks using the Newton–Raphson method. It constructs the bus
admittance matrix (Y-bus), solves bus voltage magnitudes and phase angles,
calculates bus power injections, branch power flows and network losses.

The included example is a three-bus network with one slack, one PV and one PQ
bus. The solver supports arbitrary bus counts subject to the assumptions below.

## Engineering problem
Given network topology, line impedances, generator setpoints and loads,
determine bus voltages, active/reactive injections, branch flows and losses.
Load-flow analysis supports power-system planning, operation and voltage studies.

## Mathematical model

For a branch between buses i and k, series impedance is
`z_ik = r_ik + j*x_ik` and series admittance is `y_ik = 1/z_ik`.
The Y-bus diagonal terms include series admittance and half the total line
charging susceptance at each end; off-diagonal terms are `-y_ik`.

With `Ybus = G + jB` and `theta_ik = theta_i - theta_k`:

`P_i = V_i * sum(V_k * (G_ik*cos(theta_ik) + B_ik*sin(theta_ik)))`

`Q_i = V_i * sum(V_k * (G_ik*sin(theta_ik) - B_ik*cos(theta_ik)))`

Newton–Raphson solves `J(x) * delta_x = mismatch`. The state vector contains
angles at all non-slack buses and voltage magnitudes at PQ buses. The Jacobian
uses the standard H-N-M-L block structure. Convergence is declared when the
maximum absolute P/Q mismatch is below the specified tolerance.

## Bus types and sign convention
- Slack: voltage magnitude and angle specified; P and Q are outputs.
- PV: active power and voltage magnitude specified; Q and angle are outputs.
- PQ: active and reactive power specified; voltage magnitude and angle are outputs.
- Positive P/Q means net injection into the network; loads are entered as negative.
- Bus powers and network quantities are in per-unit. `base_mva` converts powers
  to MW and MVAr.

## Repository structure
```text
python-engineering-computation/
├── power_flow.py
├── run_power_flow.py
├── requirements.txt
├── README.md
├── data/case_3bus.json
├── results/
│   ├── bus_results.csv
│   ├── branch_flows.csv
│   └── summary.json
└── tests/test_power_flow.py
```

## Requirements and execution
Python 3.10+ and NumPy are required.

```bash
python -m pip install -r requirements.txt
python run_power_flow.py
python -m unittest discover -s tests -v
```

## Sample network
The demonstration uses a 100 MVA base, a slack bus at 1.04 pu, a PV bus with
0.50 pu active-power injection and 1.02 pu voltage setpoint, and a PQ bus
with 0.60 pu active load and 0.30 pu reactive load. Branch data are illustrative,
not data from an actual utility network.

## Output
The runner writes bus voltage and injection results, branch power flows, and
a JSON summary into `results/`. The program reports iteration count, final
mismatch and total active/reactive losses.

## Verification
Automated tests check Y-bus symmetry, shunt row sums, convergence, slack/PV
voltage constraints, active-power balance, branch-loss consistency, invalid
slack-bus input rejection and finite calculated power values.

For stronger engineering validation, compare against an independent
power-system analysis package or a published benchmark case. Record the
benchmark source, base values, tolerance, solver settings and numerical
differences. Do not claim physical-grid validation without measured or
authoritative reference data.

## Scope and limitations
This version supports one slack bus, PV/PQ buses, constant specified powers,
series branch impedance and line charging. It does not yet implement generator
reactive limits with PV-to-PQ switching, transformer tap ratios, phase shifters,
distributed slack, voltage-dependent loads, thermal ratings, contingency
analysis, unbalanced networks or voltage-stability continuation.

## Future extensions
1. Generator reactive-power limits and PV/PQ switching.
2. Transformer off-nominal taps and phase shifts.
3. Branch loading percentages and thermal-limit checks.
4. Voltage-profile and loss plots.
5. Benchmark comparisons and larger-network performance tests.

## Reproducibility
The case data, sign convention, equations, code, numerical tolerance and
automated tests are included. This is an educational benchmark-style example,
not a validated model of a real power network.
