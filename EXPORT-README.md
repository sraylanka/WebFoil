# Foil Studio source export

React + Material UI airfoil analysis web app, using the Vinext/Vite Sites starter.
This archive contains the complete tracked source of the deployed first version.
Dependencies, generated builds, Git history, credentials, and local runtime state
are excluded. The original Sites hosting manifest is retained for reference.

## Run locally

Install Node.js 22.13 or newer (Node.js 24 recommended) and pnpm 11.25.0.
Extract the archive and open a terminal in `foil-studio`:

```sh
pnpm install --frozen-lockfile
pnpm dev
```

Open http://localhost:5173. On Windows, use PowerShell or WSL.
The exported source defaults to its portable execution profile because managed
runtime state is excluded. No API key is required for airfoil calculations.

```sh
pnpm build
pnpm start
```

The starter uses Cloudflare's local runtime. `pnpm start` prints its local URL.
The source is independently editable; publishing outside Sites requires configuring
hosting for this Vinext/Cloudflare Worker application.

## Main files

- app/page.tsx: Material UI workspace, geometry import/export, charts, sweeps,
  comparisons, browser-local project storage.
- app/globals.css: responsive application styling.
- lib/foil.ts: NACA generation, DAT parsing, Hess–Smith potential-flow solver.
- app/layout.tsx: page metadata and root layout.
- package.json and pnpm-lock.yaml: dependencies and reproducible installation.

## Numerical scope

Single-element, 2D, inviscid potential flow with a Kutta condition and optional
Prandtl–Glauert pressure correction (Mach 0–0.5). Pressure drag is a numerical
residual, not total drag. Boundary layers, viscous drag, transition, separation,
stall, multi-element analysis, ground effect, and full JavaFoil parity are not
implemented. Projects saved in the hosted app's browser storage are not included;
export those separately using the Projects tab.

Production build, TypeScript, and numerical sanity checks passed when this version
was created. Browser visual QA and experimental validation were not performed.
