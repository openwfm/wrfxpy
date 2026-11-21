Read the account-level AGENTS.md first: $CODEX_HOME/AGENTS.md if CODEX_HOME
is set, otherwise ~/.codex/AGENTS.md. Then read any AGENTS.md files in parent
directories, from the outermost directory toward this repository, before
applying the instructions below.

# Documentation synchronization

A code change is not complete until its documentation impact has been reviewed.

Before editing:
- Identify existing documentation and notes affected by the proposed change.
- Include their proposed updates in the same approval gate as the code.
- Classify each affected document as current guidance, historical record,
  or explicitly preserved text.

When implementing:
- Update current guidance in place.
- Preserve historical records; append a dated update when their conclusions
  or instructions need qualification.
- Do not alter explicitly preserved text.
- Document changed inputs, outputs, defaults, units, interfaces, timing,
  configuration, and limitations.
- Prefer updating an existing authoritative document over creating another.

Before reporting completion:
- Check documentation against the final implementation.
- Distinguish tested behavior from intended or unverified behavior.
- Report the documentation files updated, or explain why none were needed.
- Include implementation and corresponding documentation in the same
  approved commit, subject to existing decision-provenance rules.
