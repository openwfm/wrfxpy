# Project Instructions

## Basic Memory

- Use the `basic-memory` MCP server for cross-thread continuity in this project.
- At the start of a new task or thread, read Basic Memory for relevant project context before making assumptions.
- After meaningful progress, decisions, or handoff-worthy work, write a concise update to Basic Memory so future threads can pick up the context.
- Keep memory entries short and durable: branch, goal, decisions, constraints, and next steps.

## Active Hard Constraints Checklist (Non-Negotiable)

Before any state-changing action, the agent MUST complete this checklist in the response first.

State-changing actions include:
- Editing files (`apply_patch`, writes, deletes, moves)
- Running commands that mutate repo or outputs (`git add/commit`, generators, training runs that write artifacts)
- Any command that changes environment or persisted state

Required pre-action output (max 3 bullets):
1. Failure mode 1 (how action could fail / cause unintended outcome)
2. Failure mode 2
3. Failure mode 3

Then include:
- Why this action is still justified for current project goals
- One-line rollback/containment plan if failure occurs

No exceptions unless user explicitly says: "skip failure-mode check".

BEFORE making any changes the agent MUST stop, provide an outline of any changes, and ask to proceed after providing it. Do not skip the stop/approval gate unless I explicitly say: "do not stop before changes".

Commit policy:
- Do not commit unless I explicitly say "commit" or "amend"
- When I say "amend", move baseline with: `git commit -a --amend`
- When I say "commit", suggest a brief informative message, ask which files to change, stage thenm,upon approval do `git commit -m "<message>"`. "updated" only is not informative, say what was updated.

Do NOT change anything in submodule ml_tign.

## Ground rules of behavior (REQUIRED)

	•	Missing pre-change failure-mode check = policy violation; do not execute the change.
	•	No guessing. If you cannot verify something from the code or logs, say so and stop.
	•	Always verify the current git commit and branch before reasoning from repo history or workspace state. Do not assume what commit is checked out.
	•	One change at a time. Keep diffs minimal and scoped.
	•	Equivalent changes first. Refactors must preserve behavior unless explicitly requested.
	•	Always test before finalizing. Show the exact commands used and the key outputs.
	•	Clean pythonic python only. Do not mix python and shell with heredocs, options management, etc.


## Development Rules - REQUIRED

Changes and rewrites
	•	Make small, correct, reproducible changes.
	•	Prefer fixes that improve robustness and debuggability over adding features.
	•	Do not “rewrite” other code unless explicitly requested.
	•	No unnecesary "improvements" unless explicity requested and approved
	•	Do NOT make one-off code. Reuse or extract a common core and call it with proper values.
	•	Simplify shared core code but keep changes minimal necessary  before adding wrappers or options.
        *       Test if changes to shared core code did not break existing applications of it.
	•	Do NOT duplicate code without explicit permission.

Documentation workflow
	•	Keep sequential comprehensive notes first, even if redundant or messy.
                Record in particular the full paths of any code and artifacts
	•	Keep parallel git commit messages as a second chronology of intentional code changes.
	•	When the thinking becomes clear, promote stable understanding into continuously maintained documentation.
	•	When usage stabilizes, promote it into README recipes with no implicit knowledge and no unstated assumptions.

Build the code in small steps:
	•	Equivalent refactor to prepare base for extension
	•	Regression pass on base
	•	Minimal extension reusing core logic
	•	Regression proves base still works
	•	Independent validation of new behavior
	•	Extend regression to cover new behavior

Editing workflow (required)
	0. Pre-change gate: list up to three failure modes and why the change is still justified.
	1.	Identify the exact file(s) and line ranges to change.
	2.	Make the smallest possible patch.
	3.	Run the most local/fast test that exercises the change.
	4.	Only then consider follow-up improvements.

Changes/Commit discipline
	•	Commit only intentional changes with clear messages.
	•	Each commit should represent a coherent unit of work.
	•	Avoid mixing formatting-only changes with functional changes.
	•	Maintain up to date top README.md and separate README.md files in subdirectories
	•	Record progress with date in neural-network-notes.txt at level of detail as existing, add short commit hash.

What NOT to do
	•	Do not add complex CLI menus or new dependencies unless requested.
	•	Do not “clean up” unrelated code while fixing a bug. Minimal auditable changes only.
	•	Do not change hyperparameters (especially defaults) without a stated goal and a baseline comparison.
	•	do not duplicate logic that already has a stable owner
	•	do not add wrappers unless an existing seam is unusable
	•	do not create parallel data flows
	•	do not introduce new config surfaces without approval
	•	do not widen interfaces to solve one-off cases
	•	do not replace a simple explicit path with a generic framework
	•	preserve existing ownership boundaries between modules

# Agent operating instructions: scientific visualization projects

## Primary principle
Treat the scientific data and scientifically valid transforms as the source of truth.
Treat the visualization strictly as a human inspection and communication layer over that truth.

A visually improved result is not sufficient evidence of correctness.

## Layer hierarchy
Work in these layers and keep them distinct:

1. Scientific data layer
   - raw inputs
   - metadata, provenance, units, uncertainty
   - domain-derived quantities

2. Scientific transformation layer
   - filtering
   - interpolation
   - coordinate transforms
   - projections
   - feature extraction
   - domain-valid computations

3. Visualization mapping layer
   - mapping scientifically valid quantities to position, color, opacity, size, motion, labels, camera framing

4. UI / presentation layer
   - controls
   - layout
   - interaction
   - annotations
   - aesthetic polish

## Hard rules
- Do not change scientific meaning to improve appearance.
- Do not compensate in rendering for errors originating in data or scientific transforms.
- Fix problems at the earliest valid layer.
- Do not duplicate transform logic across scientific and rendering layers.
- Do not introduce a second path for the same scientific behavior without explicit approval.
- Do not add ad hoc exceptions to make one case “look right.”
- Do not hide outliers, discontinuities, missingness, or uncertainty unless explicitly approved.
- Do not add silent normalization, clipping, smoothing, filtering, or rescaling.
- Do not change default scaling, normalization, projection, interpolation, filtering, or coordinate conventions without explicit justification.
- Do not encode scientific policy in UI code.
- Do not introduce speculative abstractions, wrappers, or generic frameworks unless clearly necessary.

## Structural conservation rules
Optimize for minimal structural disturbance.

Preferred order of action:
1. Extend the current owner of the behavior.
2. Refactor the current owner if needed.
3. Extract a shared helper if duplication is real.
4. Only then consider a new abstraction.

Never create a parallel pipeline, shadow transform path, or convenience wrapper just to make a local task easier.

## Before making changes, state explicitly
For each proposed change, identify:

- the scientific quantity or invariant involved
- units and comparability requirements
- which layer is wrong
- which layer you will modify
- why that is the earliest valid layer
- whether the change affects scientific meaning or only presentation
- duplication risk
- any architectural boundary crossed
- the simplest rejected alternative and why it was rejected
