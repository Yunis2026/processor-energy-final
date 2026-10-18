# Coupled Electrothermal Modeling of Processor Operating Strategies

Environmental Modeling course project comparing high-frequency execution
followed by idle with low-frequency execution for an equal workload.

## Research question

Under what idle-power and thermal conditions does high-frequency execution
followed by idle reduce processor energy use and associated use-phase carbon
emissions? How sensitive is that choice to the active-current inputs?

## Requirements and quick start

- R with its standard packages and PNG graphics support.
- No additional R packages are required.
- Keep all five R scripts and paper_data.csv in the repository root.

From the repository folder, run:

```bash
Rscript run_all.R
```

The runner copies the required inputs into a clean temporary folder, executes
each analysis in a separate R process, and checks newly generated outputs.
After all checks pass, regenerated files are copied into outputs/.

Existing source files and paper_data.csv are not modified.

Successful completion prints:

```text
All four analyses and output checks passed.
```

Full stage logs are saved in outputs/logs/.
The execution summary is outputs/run_manifest.csv.

If a calculation fails, the terminal reports the temporary folder containing
diagnostics. Existing results are not published as a successful new run.

## Project files

| File | Purpose |
| --- | --- |
| paper_data.csv | Two published STM32L0 operating points |
| pde_model.R | Baseline simulations, figures and numerical checks |
| sensitivity_analysis.R | Four-factor one-at-a-time screening |
| strategy_map.R | Joint idle-power/thermal-resistance scenarios and break-even points |
| input_robustness.R | Assumed active-current perturbations and strategy crossings |
| run_all.R | Clean-workspace reproduction and cross-script output checks |

Each analysis can also run separately with Rscript followed by its filename.
For the standalone baseline script, run from the repository folder.

## Data source

The two selected operating points are transcribed from Table 1 of:

Zidar, J., Matić, T., Aleksi, I., and Hocenski, Z. (2024).
Dynamic Voltage and Frequency Scaling as a Method for Reducing Energy
Consumption in Ultra-Low-Power Embedded Systems.
Electronics, 13(5), 826.
https://doi.org/10.3390/electronics13050826

| Mode | Frequency (MHz) | Voltage (V) | Current (mA) | Task duration (ms) |
| --- | ---: | ---: | ---: | ---: |
| High | 16 | 3.3 | 3.73185 | 2.83 |
| Low | 2.097 | 3.3 | 0.51297 | 22.15 |

The workload is SHA-256 on STM32L0.

These selected points use constant voltage and therefore constitute a DFS
comparison. Both strategies complete 162528 tasks in a 3600-second observation
window with the supplied inputs.

## Model and assumptions

- A one-dimensional transient heat equation represents conduction, heat
  generation and effective distributed cooling, with insulated domain ends.
- The active reference source occupies the middle third; idle heat is uniform.
- An assumed exponential temperature-dependent increment couples power to
  temperature. It is exploratory, not calibrated leakage data.
- Reference electrical power is voltage multiplied by measured current.
- Idle power, geometry, thermal properties and feedback coefficients are assumed.
- Baseline values include idle power 0.10 mW, thermal resistance 60 K/W,
  alpha 0.20 and beta 0.05 per K.
- The fixed carbon intensity is an assumed 0.40 kg CO2e/kWh, not an estimate
  of a specified electricity grid.
- Emissions equal energy in kWh multiplied by the carbon intensity.
- The boundary includes processor use-phase electricity only.
- Manufacturing, replacement, conversion losses, peripheral power and end of
  life are excluded.

## Analysis design

1. Baseline:
   Compare equal workloads with and without temperature feedback.

2. OAT screening:
   Vary idle power 0-0.30 mW, thermal resistance 30-120 K/W,
   alpha 0-0.40 and beta 0.025-0.10 per K, with seven levels each.

3. Joint scenarios:
   Evaluate 13 idle-power values and seven thermal resistances.
   Calculate coupled break-even idle power using bracketed numerical
   root finding.

4. Input perturbations:
   Change each current alone by -5, -3, -1, 0, 1, 3 and 5 percent.
   Evaluate four joint sign combinations at magnitudes 1, 3 and 5 percent.
   Execution times and task count remain fixed.

Parameter ranges and current changes are assumed exploration scenarios,
not paper-reported error bars or confidence intervals.

No probability of strategy superiority is estimated.
The factor ranking depends on the selected ranges.

## Main results with the supplied inputs

- High-frequency execution followed by idle saves about 1.90% without
  feedback and 1.29% with feedback.
- At 60 K/W, the break-even idle power is about 0.13684 mW without feedback
  and 0.12514 mW with feedback.
- Coupled break-even idle power decreases as thermal resistance increases
  over the explored range.
- A high-mode current increase of about 1.3723%, or a low-mode current
  decrease of about 1.2927%, individually removes the baseline energy advantage.
- Some opposing joint 1% current perturbations reverse the strategy ranking.

These are conditional model results, not a universal processor operating rule.

A common fixed positive carbon intensity preserves the energy ranking.

## Numerical checks and limitations

The baseline uses conservative finite volumes, backward Euler and Picard
iteration.

Checks include:

- Energy balance.
- Agreement with the no-feedback energy formula.
- Uniform-temperature analytical comparisons.
- Grid and time-step refinement.
- Consistent baselines across analysis scripts.
- Expected output-table sizes.
- Finite, nonnegative energies.
- Energy agreement at calculated break-even points.

The main 1-second step does not resolve initial subsecond temperature dynamics.
Short-time baseline figures use a finer step.

The checks do not establish full convergence for every scenario.

Thermal parameters and feedback coefficients have not been calibrated against
STM32L0 temperature measurements.

Numerical consistency and reproducibility do not establish physical
prediction accuracy.

## Output folders

- outputs/pde/: baseline, thermal figures and numerical checks.
- outputs/sensitivity/: four-factor screening tables and figures.
- outputs/strategy_map/: joint scenarios and coupled break-even results.
- outputs/input_robustness/: current perturbation tables and figures.
- outputs/logs/: logs from the unified reproduction run.

## Course submission

The GitHub repository supplements the submission.

The course also requires runnable code uploaded to Canvas and commented code
as a PDF appendix to the final report.

A repository URL alone does not replace these requirements.