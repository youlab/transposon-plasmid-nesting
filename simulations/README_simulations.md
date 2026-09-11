# Simulation notebooks

These notebooks contain the simulations used for Figure 2, Figure 5a, and Supplementary Figures 2–4 and 11–14. Keep all notebooks in the same folder. Figure-specific notebooks load the shared model with `%run ./simulation_core.ipynb`.

- **`simulation_core.ipynb`** — Defines the coarse-grained TCN model, growth functions, transition structure, parameter sets, serial-passage protocol, initialization functions, convergence criteria, response-rate calculations, and shared simulation utilities.
- **`Figure2.ipynb`** — Generates the modeled fixed-environment and fluctuating-environment results for Figure 2 using the primary model with maximum PCN capacity `n = 100`.
- **`Figure5.ipynb`** — Generates the modeled maximum-PCN-capacity sweep for Figure 5a for transposition-on and transposition-off configurations.
- **`SF2.ipynb`** — Repeats the Figure 2 simulations with the chromosome-free boundary condition, allowing total TCN to reach zero.
- **`SF3.ipynb`** — Tests alternative transposition-associated forward-transition increments across several maximum PCN capacities.
- **`SF4.ipynb`** — Tests alternative PCN-dynamics-associated transition increments across several maximum PCN capacities.
- **`SF11.ipynb`** — Tests the maximum-PCN-capacity response using alternative normal initial TCN-state distributions.
- **`SF12.ipynb`** — Tests the maximum-PCN-capacity response using a uniform initial distribution and five independently generated random initial distributions.
- **`SF13.ipynb`** — Tests alternative TCN-growth parameterizations and the four combinations of transposition and PCN-dynamics redistribution settings.
- **`SF14.ipynb`** — Tests a specified twofold coupling between transposition status and maximum PCN capacity.
