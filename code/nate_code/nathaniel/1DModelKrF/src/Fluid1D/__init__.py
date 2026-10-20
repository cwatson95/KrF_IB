"""Declarative composition of the 1D charged-fluid, chemistry, and photon model.

Compose an experiment from component specs, numeric collectors, and run-owned plots::

    from Fluid1D import Simulation, Grid1D, GasMixture, Plasma, Cavity, Numerics
    from Fluid1D import BoltzmannSpec, CircuitSpec
    from Fluid1D import DataRef, LinePlot, OpticalGain

    Simulation(
        name="krf", grid=Grid1D(cells=128, length=2e-2, area_m2=1.0),
        gas=GasMixture({"Kr": 0.99, "F2": 0.01}), plasma=Plasma(density=1e18),
        cavity=Cavity(), boltzmann=BoltzmannSpec(backend="stub"), circuit=CircuitSpec(),
        reactions="data/reactions/krf.toml",
        diagnostics=[OpticalGain()],
        plots=[LinePlot(
            filename="gain.png",
            x=DataRef("optical_gain/times_ps", "time (ps)"),
            y=(DataRef("optical_gain/gain_per_m", "gain"),),
            ylabel="optical gain (1/m)",
        )],
    ).run(Numerics(steps=2000, dt=1e-12, fluid_every=10))

The production F2-concentration sweep (the top-level use case) is
:func:`Fluid1D.sweep.sweep_f2_parallel`; :func:`Fluid1D.sweep.sweep_f2` is the serial
validation path.
"""

from __future__ import annotations

from diagnostics import DataRef, HeatmapPanel, HeatmapPlot, LinePlot, Plot, data_limits

from .boltzmann import (
    AnalyticStubSolver,
    BoltzmannSolver,
    CoefficientTable,
    build_lut,
)
from .charged_fluid import ChargedFluid
from .components import (
    BoltzmannSpec,
    Cavity,
    CircuitSpec,
    GasMixture,
    Grid1D,
    Numerics,
    Plasma,
    Preionization,
    ProgressReport,
    Units,
)
from .convergence import (
    ConvergenceMetrics,
    ConvergenceResult,
    ConvergenceTolerance,
    compare_refinement,
    pulse_fwhm,
)
from .diagnostics import (
    Absorption,
    ChargedDensityXT,
    CircuitWaveforms,
    Diagnostic,
    EEDF,
    EnergyClosure,
    Efficiency,
    ExcimerRates,
    FieldLine,
    LaserOutput,
    OpticalGain,
    PlasmaScalarTraces,
    SpeciesDensities,
)
from .field import DischargeCircuit, Poisson1DSolver, SemiImplicitField, StubCircuit
from .energy_closure import ReconciledEnergy, reconcile_continuous_energy
from .mesh import Mesh1D
from .plasma_model import DriftDiffusionPlasmaModel, PlasmaModel, StepReport, StubPlasmaModel
from .simulation import Simulation
from .sweep import SweepFailure, SweepPoint, SweepResult, sweep_f2, sweep_f2_parallel

__all__ = (
    # specs
    "Grid1D",
    "GasMixture",
    "Plasma",
    "Preionization",
    "Cavity",
    "Numerics",
    "BoltzmannSpec",
    "CircuitSpec",
    "ProgressReport",
    "Units",
    "ConvergenceMetrics",
    "ConvergenceResult",
    "ConvergenceTolerance",
    "compare_refinement",
    "pulse_fwhm",
    # engine
    "Simulation",
    "Mesh1D",
    "ChargedFluid",
    # plasma model
    "DriftDiffusionPlasmaModel",
    "PlasmaModel",
    "StubPlasmaModel",
    "StepReport",
    "ReconciledEnergy",
    "reconcile_continuous_energy",
    # boltzmann
    "BoltzmannSolver",
    "AnalyticStubSolver",
    "CoefficientTable",
    "build_lut",
    # field
    "DischargeCircuit",
    "StubCircuit",
    "Poisson1DSolver",
    "SemiImplicitField",
    # diagnostics
    "Diagnostic",
    "EEDF",
    "EnergyClosure",
    "SpeciesDensities",
    "OpticalGain",
    "Absorption",
    "LaserOutput",
    "Efficiency",
    "ChargedDensityXT",
    "FieldLine",
    "PlasmaScalarTraces",
    "CircuitWaveforms",
    "ExcimerRates",
    # plots (filenames and datasets are declared by experiment files)
    "Plot",
    "DataRef",
    "LinePlot",
    "HeatmapPanel",
    "HeatmapPlot",
    "data_limits",
    # sweep
    "sweep_f2",
    "sweep_f2_parallel",
    "SweepResult",
    "SweepPoint",
    "SweepFailure",
)
