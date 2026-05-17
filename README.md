# Oil-Pipeline-Pump-Operations-Optimization

To be updated further (the actual work has been completed).

This work is implemented in GAMS (.gms files), a programming language specifically designed for large-scale optimization problems. One advantage of GAMS is that the decision variables across the entire optimization horizon (i.e., all variables over all time steps) can be formulated and visualized more systematically compared with Python, where handling such large-scale indexed optimization structures can become considerably more cumbersome.

However, GAMS itself is not well suited for visualization and plotting. Therefore, all figures and result visualizations are generated through Python scripts that read the .gdx output files produced by GAMS.

The general framework and methodology are generally illustrated in the presentation slides (202509 Presentation.pptx). Three cases are presented there. The fixed 6 ppm and 12 ppm cases (where ppm represents the concentration of DRA) are primarily included for comparison purposes, demonstrating the advantage of an optimized schedule for pump operations and DRA injections compared with traditional fixed-injection industrial practices.

After reformulation and linearization, the model becomes a MILP (Mixed-Integer Linear Programming) problem solved using the Gurobi solver. The formulation involves two major types of linearization:

Piecewise linearization for nonlinear quadratic-type relationships (e.g., squared terms);
SOS2-based linearization for colinear bilinear relationships involving variable–variable interactions.
