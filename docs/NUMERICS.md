# Numerical methods, references and known limits

Run `pnpm test` (Node 24+) for the numerical verification suite. The browser itself does not require Node. Tests use independently known results where possible, not only comparisons between implementations. Passing tests establish these bounded properties; they do not establish general engineering accuracy.

## Existing Hess–Smith solver

Source: `lib/foil.ts`. Constant source strength on each straight panel, one constant vortex strength per element, midpoint no-penetration equations and a trailing-edge Kutta equation. The linear solve uses partial pivoting. Cp = (1 − Vt²)/sqrt(1 − M²) applies the approximate Prandtl–Glauert pressure correction only (Mach ≤ 0.5). Velocity fields remain incompressible. Lift and pressure-drag residual are pressure integrals; moment is nose-up about quarter chord. Input contours must be finite, closed, counterclockwise and non-intersecting.

Verified: symmetric zero lift/moment, odd lift/moment symmetry, even drag symmetry, cambered lift/moment signs, thin symmetric lift close to 2π radians, mirror invariance, source-flux refinement, pressure scaling and exact single-element agreement with the generalized solver.

Independent circular-cylinder check at 200 panels: Cp RMS error approximately 0.00060 relative to 1 − 4 sin² θ, with zero lift and drag to roundoff.

**Known accuracy limit:** a mapped Joukowski airfoil (circle center −0.08 + 0.08i, cusp mapping z = ζ + 1/ζ, normalized chord) at α = 4° gives CL approximately 0.85457 / 0.89840 / 0.91694 at 100 / 200 / 300 panels, versus analytical CL approximately 0.96943. At 300 panels the relative error remains 5.4%. This test checks that the error decreases and remains below 6%; it is a documented limit, not evidence of precision. Sharp cusps and poorly resolved trailing edges need a higher-order panel method or better Kutta discretization. Pressure drag is still approximately −0.00625 in this case and must never be interpreted as profile drag. No experimental data comparison is claimed.

## Generalized inviscid model

`lib/advanced.ts` adds a coupled source/vortex system with a separate circulation unknown and Kutta condition for each of up to three elements. Loads reference the main element's unit chord and quarter-chord; total panels are limited to 500. Intersecting element edges are rejected. This is experimental: no independent multi-element load benchmark has yet been run. Nested elements are rejected; nearly touching elements are outside the validated scope.

Ground effect uses equal source and opposite vortex image strengths across y = 0. Freestream must be parallel to the wall; incidence is produced by rotating the geometry clockwise. Geometry must stay strictly above the wall. Tests check zero normal wall velocity and recovery of free-air loads at large clearance. This is an inviscid slip wall, not a ground boundary layer or moving ground model.

Streamlines integrate the incompressible velocity field with normalized midpoint steps. They stop at solids/domain bounds and have fixed work limits. They do not model separated wakes or turbulence.

## Experimental boundary layers

Implemented independently; **XFOIL source code has not been copied or ported**. XFOIL is a reference for architecture, applicable physics and limits, not a numerical-equivalence claim.

- Laminar momentum thickness: Thwaites, θ² = 0.45/Re · Ue⁻⁶ ∫ Ue⁵ ds.
- Laminar shape/friction use Thwaites correlations. λ ≤ −0.09 flags laminar separation.
- Free transition uses Michel's empirical Reθ criterion; optional chordwise forced trips. This is **not XFOIL's envelope e^N transition method**, and no Ncrit control is offered.
- Turbulent momentum equation and Head entrainment/shape relation, with Ludwieg–Tillmann skin friction. H ≥ 2.4 flags turbulent separation. Integration is one-way using the inviscid external velocity; no viscous displacement feedback is solved.
- Attached flow drag uses Squire–Young trailing-edge extrapolation: each side contributes 2 θ Ue^((H+5)/2).
- Separation terminates that surface's integration and withholds total viscous drag. No fabricated stall correction, bubble reattachment or post-stall lift is supplied. Lift/moment remain inviscid even when a viscous estimate is requested.

Tests compare laminar flat-plate momentum thickness/Cf with Blasius and turbulent flat-plate drag with the empirical 0.074 Re^−1/5 correlation (20% tolerance), check symmetry and physical ranges, and require withheld drag after separation. These broad checks do not validate pressure-gradient airfoil polars. In particular the current model often flags trailing-edge separation, so many runs intentionally have no drag result. Do not use these estimates to rank airfoils for engineering decisions without independent validation.

Reference consulted: [XFOIL source mirror](https://github.com/RobotLocomotion/xfoil), revision `d11a1544b53623c01bb3b0cceb5862311be0e1f8`, `xfoil_doc.txt` viscous formulation/transition/drag sections and `src/xbl.f`, `src/xblsys.f`. [Mark Drela's XFOIL page](https://web.mit.edu/drela/Public/web/xfoil/) and Drela & Giles (1987), “Viscous-Inviscid Analysis of Transonic and Low Reynolds Number Airfoils.” The MIT host was inaccessible from this environment; documentation was read from the source mirror. Full coupled XFOIL or MSES fidelity is not implemented.

## Geometry and design

Non-reflex NACA five-digit families use standard tabulated camber parameters; reflex families are rejected. CST uses a sqrt(x)(1 − x) class function and Bernstein coefficients, with a closed trailing edge. A blended flap is a geometry approximation, not an exact hinged flap or slotted element. Smoothing changes y only and preserves the leading/trailing edges.

Lift-target design adjusts one parabolic camber mode using finite-difference Newton updates, with a 0.001 CL residual and an 8% added-camber bound. It does not solve full inverse design for a specified pressure distribution. JSON batches are validated bounded data jobs; no JavaScript is evaluated.

## Follow-up work needed for full capability parity

A coupled viscous/inviscid solver or audited XFOIL WebAssembly port (including GPL obligations), envelope e^N transition, separated-flow/stall validation, higher-order trailing-edge treatment, full pressure-target inverse design, and independent multi-element/ground-effect datasets remain necessary. A near-stall XFOIL result itself is not ground truth.
