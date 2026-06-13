# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 12:39:51 2026

@author: anvifer2

Solución analítica del reactor placa 1G heterogéneo (reflector-combustible-reflector).
Calcula TODOS los modos (pares e impares) resolviendo dos ecuaciones de criticidad:
  - Paridad par  (cos en el núcleo): Dc·Bc·tan(Bc·a) = Dr·κ_r·coth(κ_r·b)
  - Paridad impar (sin en el núcleo): -Dc·Bc·cot(Bc·a) = Dr·κ_r·coth(κ_r·b)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import warnings

warnings.filterwarnings('ignore')
plt.close('all')

# ==========================================
# 1. Define Material Properties & Dimensions
# ==========================================
N_plot = 200          # Points for plotting

# FUEL properties
Dc = 0.776            # Diffusion coefficient in core (cm)
Sigma_ac = 0.0244     # Absorption cross-section in core (cm^-1)
nu_Sigma_f = 0.0260   # Production cross-section in core (cm^-1)

# Reflector properties
Dr = 1.446            # Diffusion coefficient in reflector (cm)
Sigma_ar = 0.0077     # Absorption cross-section in reflector (cm^-1)

# Geometry (symmetric about x=0):
#   REFLECTOR | CORE | REFLECTOR
#   <- b ->   <-a->a  <- b ->
a = 150.0             # Core half-width (cm)
b = 25.0              # Reflector thickness (cm)
L = 2 * (a + b)       # Total reactor length (cm)

# ==========================================
# 2. Reflector Parameters
# ==========================================
kappa_r = np.sqrt(Sigma_ar / Dr)
coth_kb = 1.0 / np.tanh(kappa_r * b)
rhs = Dr * kappa_r * coth_kb          # Common right-hand side

# ==========================================
# 3. Criticality Equations
# ==========================================
# Even parity (cos modes → modes 1, 3, 5, …):
#   Dc·Bc·tan(Bc·a) = rhs
def f_even(Bc):
    return Dc * Bc * np.tan(Bc * a) - rhs

# Odd parity (sin modes → modes 2, 4, 6, …):
#   -Dc·Bc·cot(Bc·a) = rhs
def f_odd(Bc):
    return -Dc * Bc / np.tan(Bc * a) - rhs

# ==========================================
# 4. Root Finding (scan for sign changes)
# ==========================================
num_modes_wanted = 3
Bc_max = num_modes_wanted * np.pi / a   # enough range for desired modes
Bc_scan = np.linspace(1e-8, Bc_max, 50000)

def find_roots(func, Bc_grid):
    """Find roots of func by detecting sign changes, avoiding poles."""
    roots = []
    for i in range(len(Bc_grid) - 1):
        v1, v2 = func(Bc_grid[i]), func(Bc_grid[i + 1])
        if np.isfinite(v1) and np.isfinite(v2) and v1 * v2 < 0:
            # At a genuine root the function crosses smoothly (- → +  or + → -).
            # At a pole the function jumps (e.g. +∞ → -∞), so |v1|,|v2| are huge.
            # We also verify monotonicity: at a root |v| should be small near zero.
            if abs(v1) < 50 and abs(v2) < 50:
                try:
                    root = brentq(func, Bc_grid[i], Bc_grid[i + 1],
                                  xtol=1e-14, rtol=1e-14)
                    # Double-check: the function value should be near zero
                    if abs(func(root)) < 1e-8:
                        roots.append(root)
                except ValueError:
                    pass
    return roots

even_roots = find_roots(f_even, Bc_scan)
odd_roots = find_roots(f_odd, Bc_scan)

# ==========================================
# 5. Combine and Sort by k_eff (descending)
# ==========================================
# k_eff = nu_Sigma_f / (Dc·Bc² + Sigma_ac)
all_modes = []
for Bc in even_roots:
    k = nu_Sigma_f / (Dc * Bc**2 + Sigma_ac)
    all_modes.append(('even', Bc, k))
for Bc in odd_roots:
    k = nu_Sigma_f / (Dc * Bc**2 + Sigma_ac)
    all_modes.append(('odd', Bc, k))

all_modes.sort(key=lambda m: m[2], reverse=True)

print(f"{'Mode':>4s} | {'Parity':>6s} | {'Bc':>10s} | {'k_eff':>14s}")
print("-" * 45)
for idx, (parity, Bc, k) in enumerate(all_modes):
    print(f"  {idx+1:2d}  | {parity:>6s} | {Bc:10.6f} | {k:14.10f}")

# ==========================================
# 6. Build Flux Profiles
# ==========================================
x_sym = np.linspace(-(a + b), a + b, N_plot)   # symmetric domain

plt.figure(figsize=(12, 6))
num_to_plot = min(len(all_modes), 6)

for idx in range(num_to_plot):
    parity, Bc, k = all_modes[idx]

    A_core = 1.0
    if parity == 'even':
        C_R = A_core * np.cos(Bc * a) / np.sinh(kappa_r * b)
        C_L = C_R                        # symmetric
        core_phi = A_core * np.cos(Bc * x_sym)
    else:
        C_R = A_core * np.sin(Bc * a) / np.sinh(kappa_r * b)
        C_L = -C_R                       # antisymmetric
        core_phi = A_core * np.sin(Bc * x_sym)

    left_refl = C_L * np.sinh(kappa_r * (x_sym + a + b))
    right_refl = C_R * np.sinh(kappa_r * (a + b - x_sym))

    phi = np.where(x_sym <= -a, left_refl,
          np.where(x_sym >= a,  right_refl, core_phi))

    # Normalize to max |φ| = 1
    phi /= np.max(np.abs(phi))

    # Shift domain to [0, L]
    x_plot = x_sym + a + b
    plt.plot(x_plot, phi, linewidth=2,
             label=f'Mode {idx+1} ({parity}), $k_{{eff}}$={k:.7f}')

# Interface lines
plt.axvline(b, color='gray', linestyle='--', alpha=0.7, label='Core-Reflector Interface')
plt.axvline(2 * a + b, color='gray', linestyle='--', alpha=0.7)

plt.title("Analytical Neutron Flux — Reflected Slab Reactor (All Modes)", fontsize=14)
plt.xlabel("Position x (cm)", fontsize=12)
plt.ylabel("Normalized Flux $\\phi(x)$", fontsize=12)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=9)
plt.tight_layout()
plt.show()
