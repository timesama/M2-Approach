from dataclasses import dataclass, field
import numpy as np

@dataclass
class SpectrumState:
    frequency: np.ndarray | None = None
    re_spectra: np.ndarray | None = None
    im_spectra: np.ndarray | None = None

@dataclass
class AppState:
    current_tab: str = 'SE'
    spectrum: SpectrumState = field(default_factory=SpectrumState)
    se_files: list[str] = field(default_factory=list)
    dq_files: list[str] = field(default_factory=list)
    t1t2_files: list[str] = field(default_factory=list)
    gs_files: list[str] = field(default_factory=list)
    dqmq_files: list[str] = field(default_factory=list)
    glycerol_files: list[str] = field(default_factory=list)
    baseline_files: list[str] = field(default_factory=list)
    dq_t2: dict = field(default_factory=dict)
    tau_dictionary: dict = field(default_factory=dict)
    gs_dictionary: dict = field(default_factory=dict)
