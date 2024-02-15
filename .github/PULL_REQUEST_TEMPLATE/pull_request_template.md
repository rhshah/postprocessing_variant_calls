## PR checklist

Closes #XXX <!-- If this PR fixes an issue, please link it here! -->

- [ ] This comment contains a description of changes (with reason).
- [ ] You've followed the CONTRIBUTING.md.
- [ ] You have fully documented `--help` for any new CLI commands / sub-commands.
- [ ] If you've fixed a bug or added code that should be tested, add tests!
- [ ] Remove all TODO statements.
- Ensure that following tests pass locally and in git-actions:
  - [ ] `pytest tests`
  - [ ] `black --check .`
