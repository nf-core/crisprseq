# nf-core/crisprseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v2.1.0dev]

### Added

- Template update v2.9 ([#52](https://github.com/nf-core/crisprseq/pull/52))
- Use `Channel.fromSamplesheet()` from `nf-validation` to validate input sample sheets and create an input channel ([#58](https://github.com/nf-core/crisprseq/pull/58))

### Fixed

- Summary processes don't modify the input file anymore, allowing resuming these processes ([#66](https://github.com/nf-core/crisprseq/pull/66))
- Do not stash unexistent files, use empty lists instead. Fixes AWS tests ([#67](https://github.com/nf-core/crisprseq/pull/67))
- Rename process `merging_summary` to `preprocessing_summary` to improve clarity ([#69](https://github.com/nf-core/crisprseq/pull/69))

### Deprecated

## [v2.0.0 - Paprika Lovelace](https://github.com/nf-core/crisprseq/releases/tag/2.0.0) - [05.07.2023]

### Added

- Crisprseq screening analysis : mageck mle, mageck rra, mageck count and crisprcleanr-normalize ([#22]https://github.com/nf-core/crisprseq/pull/22))
- Add new parameter `--analysis` to select analysis type (screening/targeted) ([#27](https://github.com/nf-core/crisprseq/pull/27))
- Tests to run screening analysis ([#926]https://github.com/nf-core/test-datasets/pull/926)
- Metro map for targeted analysis ([#35](https://github.com/nf-core/crisprseq/pull/35))
- Add new parameters `--reference` and `--protospacer` ([#45](https://github.com/nf-core/crisprseq/pull/45))
- Add UMI clustering to crisprseq-targeted ([#24](https://github.com/nf-core/crisprseq/pull/24))
- Template update v2.8 ([#21](https://github.com/nf-core/crisprseq/pull/21))

### Fixed

- Fix warning "module used more than once" ([#25](https://github.com/nf-core/crisprseq/pull/25))
- Remove profile `aws_tower`, update from nf-core/tools v2.9 to make full tests pass ([#50](https://github.com/nf-core/crisprseq/pull/50))
- Remove `quay.io` from all containers ([#50](https://github.com/nf-core/crisprseq/pull/50))
- Fix resources on test screening full ([#50](https://github.com/nf-core/crisprseq/pull/56))
- Fix id in crisprcleanr-normalize ([#50](https://github.com/nf-core/crisprseq/pull/56))

## [v1.0 - Salted Hypatia](https://github.com/nf-core/crisprseq/releases/tag/1.0) - [02.02.2023]

Initial release of nf-core/crisprseq, created with the [nf-core](https://nf-co.re/) template.
