# hodgkin-huxley

**Interactive Hodgkin–Huxley Neuron Model** — browser-based educational simulator.

Part of the [A. Mirza academic tools portfolio](https://dthornz.github.io/website-cv-tools/).

🌐 **Live:** [dthornz.github.io/hodgkin-huxley](https://dthornz.github.io/hodgkin-huxley/)

---

Fully interactive browser simulation of the Nobel Prize-winning 1952 Hodgkin–Huxley model. Implements 4th-order Runge–Kutta integration (dt = 0.01 ms) on the absolute voltage scale (V_rest = −65 mV, peak AP ≈ +40 mV).

**Simulator features:**
- 6 real-time Canvas plots: V(t), gating variables m/h/n, ionic currents I_Na/I_K/I_L, steady-state curves x∞(V), time constants τx(V), V–n phase portrait
- Stimulus modes: step, pulse train, ramp, sinusoidal
- All parameters adjustable: ḡNa, ḡK, ḡL, ENa, EK, EL, Cm, temperature (Q₁₀ scaling)
- Channel blockers: Na⁺ (simulate TTX), K⁺ (simulate TEA)
- Speed control (⅛× to 14×), auto-pause, AP spike counter, peak voltage, firing rate
- Phase Demo preset — auto-runs a limit-cycle orbit in the V–n phase portrait

**Educational content:**
- Complete derivation of all 4 ODEs with equivalent circuit diagram (SVG)
- Rate functions α/β for m, h, n gates with empirical constants
- Parameter table with physiological interpretations
- 7 cited references (Hodgkin & Huxley 1952, Goldman UC Davis, Cambridge DAMTP, FSU Bertram, etc.)

**Tech:** Vanilla HTML/CSS/JS · Canvas API · RK4 ODE solver · No build step · GitHub Pages

---

© 2026 Asad Mirza, Ph.D. · Research Assistant Professor · FIU Biomedical Engineering
