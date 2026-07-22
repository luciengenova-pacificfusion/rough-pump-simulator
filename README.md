# Pump Vibration Isolation Simulation (Octave)

Two-degree-of-freedom (2-DOF) vibration isolation model for a pump assembly mounted on passive isolators above a receiving structure (lab optical table, DPM airbox, or DPM). The script compares candidate isolators by simulating time-domain acceleration responses and Welch RMS velocity spectra against VC (Vibration Criterion) limits.

## Requirements

- GNU Octave
- `control` package (`pkg load control`)
- `signal` package (`pkg load signal`)

Install packages if needed:

```octave
pkg install -forge control
pkg install -forge signal
```

## Model Overview

The system is modeled as two lumped masses connected in series by spring–damper pairs:

```
        F(t)  (pump forcing)
         │
   ┌─────▼─────┐
   │    m1     │   Pump + base + heat exchanger
   └─────┬─────┘
       k1, c1     Isolator under test (Sorbothane / air shock / wire rope)
   ┌─────▼─────┐
   │    m2     │   Receiving structure (optical table / airbox / DPM)
   └─────┬─────┘
       k2, c2     Support stiffness/damping (optical table legs, or 0 for worst case)
   ═════════════   Ground
```

- **m1** — pump (75 kg) + heat exchanger (7 kg) + base plates (2x at 22 kg each, edit to 1\*22 to simulate single baseplate)
- **m2** — selected by `setup` (optical table top, airbox, or DPM)
- **k1, c1** — isolator stiffness and damping (per-isolator lists, multiplied by mount count)
- **k2, c2** — optical table leg stiffness/damping (4 legs); can be zeroed for a worst-case rigid-mount comparison

Forcing is a sum of three sinusoids representing pump excitation:

| Component | Amplitude | Frequency |
|-----------|-----------|-----------|
| F1 | 60 N | 112 Hz |
| F2 | 20 N | 62 Hz |
| F3 | 60 N | 224 Hz |

### Transfer functions

Transfer functions are built symbolically with `tf('s')` from the coupled equations of motion:

```
m1·ẍ1 + c1·(ẋ1 − ẋ2) + k1·(x1 − x2) = F(t)
m2·ẍ2 + c1·(ẋ2 − ẋ1) + k1·(x2 − x1) + c2·ẋ2 + k2·x2 = 0
```

The script implements:

**Receiving structure displacement per pump displacement (transmissibility):**

```
X2(s)          c1·s + k1
───── = ─────────────────────────────
X1(s)   m2·s² + (c1+c2)·s + (k1+k2)
```

**Receiving structure displacement per unit pump force:**

```
X2(s)                         c1·s + k1
───── = ────────────────────────────────────────────────────────────────
F(s)    (m1·s² + c1·s + k1)·(m2·s² + (c1+c2)·s + (k1+k2)) − (c1·s + k1)²
```

**Derived transfer functions:**

| Variable | Definition | Description |
|----------|------------|-------------|
| `x2A` | X2/F above | Receiving-structure displacement per unit pump force |
| `x1A` | `x2A / x2x1` | Pump displacement per unit pump force |
| `v1A`, `v2A` | `s·x1A`, `s·x2A` | Velocity transfer functions |
| `a1A`, `a2A` | `s²·x1A`, `s²·x2A` | Acceleration transfer functions |

Time responses are computed with `lsim` over a 10 s window at 10 kHz sampling.

## Configuration

Edit these variables at the top of the script:

### `isoSelect` — isolator choice

| Value | Isolator |
|-------|----------|
| 1 | Sorbothane (6 mounts) |
| 2 | Air shock (4 mounts) |
| 3 | Wire rope (4 mounts) |
| 4 | All three (overlaid on plots) |

### `setup` — receiving structure (m2)

| Value | Setup | Mass |
|-------|-------|------|
| 1 | Lab optical table test setup | ~399 kg (877.8 lb / 2.2) |
| 2 | DPM airbox only | 4000 kg |
| 3 (or other) | Full DPM | 16000 kg |

### Isolator parameter sets

`k1_list` (stiffness) and `z1_list` (damping ratio) each have three variants in the script — researched values, static-deflection-derived values, and match-tuned values. Uncomment the desired row. Match-tuned values are uncommented by default as these most closely match real world test data. Damping coefficients `c1_list` are derived from the damping ratios via `c = 2·ζ·√(k·m1)` (valid single-mass approximation since m1 ≪ m2).

### Worst-case model

Uncomment `k2 = 0; c2 = 0;` to remove the optical-table-leg isolation and treat m2 as free floating in space (a worst case assumption for large structure modeling).

## Outputs

### Figure 2 — Acceleration time response
- Top subplot: pump (m1) acceleration vs. time
- Bottom subplot: receiving structure (m2) acceleration vs. time
- One trace per selected isolator

### Figure 5 — Welch RMS velocity spectra
- Top subplot: pump velocity spectrum
- Bottom subplot: receiving-structure velocity spectrum
- Computed via `pwelch` (Hanning window, NFFT = 2¹⁵, 50 % overlap), scaled to RMS velocity per bin (`√(PSD·Δf)`) in µm/s
- Dashed reference lines mark VC limits:
  - **VC-C** = 12.5 µm/s RMS
  - **VC-D** = 6.25 µm/s RMS
  - **VC-E** = 3.12 µm/s RMS
- X-axis: 0–1000 Hz, linear

### Optional (commented-out) figures
- Figure 1 — displacement time response
- Figure 3 — velocity time response with VC-C peak markers
- Figure 4 — frequency response (Bode magnitude) with operating-frequency marker
- Console printout of gain/phase at the operating frequency

Uncomment the corresponding blocks at the bottom of the script to enable them.

## Usage

1. Open the script in Octave.
2. Set `isoSelect` and `setup` as desired.
3. Choose the isolator parameter set (comment/uncomment `k1_list` / `z1_list` rows).
4. Run the script. Figures 2 and 5 will be generated.

## Interpretation

An isolator passes the vibration criterion if the receiving-structure (m2) RMS velocity spectrum stays below the relevant VC line (e.g., VC-C at 12.5 µm/s) across the band of interest, particularly at the pump forcing frequencies (62, 112, 224 Hz) and their harmonics.

## Notes & Assumptions

- All units SI unless noted (µm/s used only for plotting VC comparisons).
- Damping coefficients calculated from damping ratios assume m1 dominates the effective isolator mass (m_eff = m1·m2/(m1+m2) ≈ m1 when m1 ≪ m2).
- Forcing amplitudes and frequencies are estimates tuned to match accelerometer test data
- The Welch spectrum resolution is Δf = fs/NFFT = 10000/32768 ≈ 0.305 Hz.
