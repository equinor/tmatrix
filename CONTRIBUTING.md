# Contribution Guidelines

Instructions and guidelines for how you can contribute to this repository.

## Issues

For reporting bugs, requesting new features and more do so by [creating a new issue][new-issue]. Make sure that:

- The issue clearly states expected vs current behavior
- Includes steps to reproduce bugs
- Look through the [existing issues][existing-issues] to see if your issue or suggestion has already been reported, if it has, add your own context to that issue instead of creating a new one.

### Pull-Requests

- should only solve a single problem/issue
- can be in the form of a single or multiple commits
- must always be up to date with the `main` branch when:
  - the PR is created
  - before the PR is merged
- must explain **why** and **how** of the proposed changes in the PR description.
- ensure the commit history is clean

#### Commits

- Use [conventional commit messages][conventional-commits], as they are used to generate the changelog and perform releases using [release please][release-please].
- A single commit should only contain a *single logical change*
  - Commit should be as small as possible
  - Commit should contain a complete change

#### Review

A PR must be approved by a *Code Owner*. Once a PR is marked as ready for review _and_ passes all required checks a code owner will look through the PR as soon as they are able to.

If the PR is approved you can either `squash and merge` or `rebase and merge` it with main depending on commits/changes in the PR. If the reviewer had comments, or the request is not approved, address and resolve these comments before re-requesting a review.

<!-- External Links -->
[new-issue]: https://github.com/equinor/tmatrix/issues/new/choose
[existing-issues]: https://github.com/equinor/tmatrix/issues
[release-please]: https://github.com/googleapis/release-please?tab=readme-ov-file#release-please
[conventional-commits]: https://www.conventionalcommits.org/en/v1.0.0/
