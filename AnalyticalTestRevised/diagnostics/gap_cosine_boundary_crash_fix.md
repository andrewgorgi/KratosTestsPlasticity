# Gap-SBM geometry crash, 2026-10-13

Case: cosine_boundary, plastic, degree 1, 10x10 spans, Kratos/bin/Release.

The native backtrace of the original input stops in
`SnakeGapSbmProcess::CreateConditions`. The JSON requested unregistered
`CutSbmSolidCondition` / `CutSbmLoadSolidCondition` operators. The current
Kratos application registers `GapSbm…` operators instead. A FullDebug run
reports the missing component explicitly; Release dereferences the missing
condition prototype and crashes. The active JSON now uses GapSbm condition,
element and interface names.

With those names corrected, initialization exposed a second native crash in
`GapSbmLoadSolidCondition::InitializeSbmMemberVariables`: the 2D geometry has
no PROJECTION_NODE. The condition now uses the same null-safe fallback to the
boundary geometry center as GapSbmSolidCondition. This change was applied in
Kratos/applications/IgaApplication/custom_conditions/gap_sbm_load_solid_condition.cpp,
built with target KratosIgaCore, and installed in Kratos/bin/Release/libs.
The portable diff is ../kratos_gap_load_projection_fix.patch. The previous
installed library is backed up at /tmp/gap_crash_audit/libKratosIgaCore.before.so.
FullDebug was not rebuilt.

The current GapSbmLoadSolidCondition reads FORCE from the condition rather
than hardcoding the manufactured traction. FORCE expressions were therefore
added on bottom/top using the existing analytical stress, with outward
normals (0,-1) and (0,1). The displacement and body-force expressions,
material, solver tolerances, timestep and quadrature were preserved.

Verification:

- Degree 1, 10x10 initialization passes with the corrected Release library.
- Prescribed FORCE was checked against analytic_fields.evaluator at all 60
  bottom/top condition centers at t=0.00125 (tolerance 1e-10).
- Full 24-thread solve reaches Newton, then returns false after 100 iterations
  at t=0.00125. See gap_cosine_boundary_p1_crash_fix.log. This is a separate
  numerical convergence failure; no successful plastic solve is claimed.
- study.Analysis now raises on a failed SolveSolutionStep before finalization
  and sampling, preventing invalid convergence/error output.
- --refinements is respected; the existing effective default 10..59 remains.

Reproduce from gap_sbm (with the Release PYTHONPATH and LD_LIBRARY_PATH):

```bash
python3 convergence.py --physics plastic --field cosine_boundary --degree 1 --threads 24 --refinements 10
```

This command currently ends with `Newton non convergente al tempo 0.00125`,
not the geometry initialization segmentation fault. Successful end-of-step
finalization and a convergence study have not been validated.
