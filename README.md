# htslib-install-script

A convenience script to install any version of `samtools` on an Ubuntu (`22.04`) or Debian system (`bookworm`).

## Table of Contents

- [Description](#description)
- [Installation - quick start](#installation---quick-start)
- [Installation - quick start with libdeflate](#installation---quick-start-with-libdeflate)
- [Installation - controlling the install location](#installation---controlling-the-install-location)
- [Requirements](#requirements)
- [Testing](#testing)
- [Development](#development)

## Description

The script `install_samtools.sh` is a convenience script to install
`samtools` ([GitHub: samtools/samtools](https://github.com/samtools/samtools)),
a suite of programs for interacting with high-throughput sequencing data.

The script encapsulates the steps to download, configure, compile and install
`samtools` to a specified location, for versions `1.10` to `1.20`. It does not install htslib (or tabix).

The script is tested via a private GitLab CICD against Ubuntu 22.04 and Debian
bookworm with popular Docker images.

## Installation - quick start

With a single command the script can be downloaded and installed. For details on how to install to a custom location, 
see [Installation - controlling the install location](#installation---controlling-the-install-location) section.

It is recommended to install both `libdeflate` for faster decompression of BAM/CRAM files and install `htslib` so that you have access to the `tabix` and `bgzip` tools. See the [quick start with recommended extras](#installation---quick-start-with-recommended-extras) section.

For a Dockerised example of how to install `samtools`, see the `docker/Dockerfile.ubuntu22.via_github`.

Various version of the script can be downloaded from the [releases page](https://github.com/team113sanger/htslib-install-script/releases).

```bash
SAMTOOLS_VERSION="1.16"
SAMTOOLS_SCRIPT_URL="https://github.com/team113sanger/samtools-install-script/releases/download/1.0.1/install_samtools.sh"

curl -sSL $SAMTOOLS_SCRIPT_URL | bash -s -- $SAMTOOLS_VERSION

# or with wget if curl is not available
wget -qO- $SAMTOOLS_SCRIPT_URL | bash -s -- $SAMTOOLS_VERSION
```

## Installation - quick start with recommended extras

Typically you want samtools with tabix and bgzip (which are provided by
installing htslib) and libdeflate for faster decompression of BAM/CRAM files.

For a Dockerised example of how to install `samtools` with both `libdeflate` and `htslib`, see the `docker/Dockerfile.ubuntu22.via_github.with_htslib_and_libdeflate`.


```bash
LIBDEFLATE_VERSION="v1.9"
HTSLIB_VERSION="1.16"
SAMTOOLS_VERSION="1.16"
LIBDEFLATE_SCRIPT_URL="https://github.com/team113sanger/libdeflate-install-script/releases/download/1.0.1/install_libdeflate.sh"
HTSLIB_SCRIPT_URL="https://github.com/team113sanger/htslib-install-script/releases/download/1.0.1/install_htslib.sh"
SAMTOOLS_SCRIPT_URL="https://github.com/team113sanger/samtools-install-script/releases/download/1.0.1/install_samtools.sh"

curl -sSL $LIBDEFLATE_SCRIPT_URL | bash -s -- $LIBDEFLATE_VERSION
curl -sSL $HTSLIB_SCRIPT_URL | bash -s -- $HTSLIB_VERSION
curl -sSL $SAMTOOLS_SCRIPT_URL | bash -s -- $SAMTOOLS_VERSION --use-installed-htslib
```

**Note 1**: The `libdeflate` script is run first as `htslib` will use `libdeflate` if it is available.
**Note 2**: The `--use-installed-htslib` flag is used to tell the `samtools` script to re-use the installed `htslib`.

For more information on how to install `libdeflate`, see the [libdeflate-install-script](https://github.com/team113sanger/libdeflate-install-script) and on how to install `htslib`, see the [htslib-install-script](https://github.com/team113sanger/htslib-install-script).

## Installation - controlling the install location

**The easiset way is to look at the Dockerfiles in the repository** as this is tested and under CI.

But in general, you can run the following commands to install libdeflate which will install to `/usr/local`:

```bash
bash install_samtools.sh 1.16
```

Or you can specify a different install location e.g. `/path/to/install`:
```bash
DEST_DIR=/path/to/install

bash install_samtools.sh 1.16 --install-dir $DEST_DIR

export PATH=$DEST_DIR/bin:$PATH
```

## Requirements

The common requirements for the script and the installation of `samtools` are documented
clearly in the Dockerfiles in the repository. **It is strongly recommended to install `libdeflate` and `htslib`**
though it is optional. See the [quick start with recommended extras](#installation---quick-start-with-recommended-extras) section.

`samtools` officially documents its required and optional system requirements in
the [INSTALL](https://github.com/samtools/samtools/blob/develop/INSTALL) file of
its repository.

If a required system library is missing compilation will fail with a
self-explanatory error message e.g. `configure: error: libbzip2 development files not found ...`.


## Testing

The testing of script is done using Docker images to capture the minimal installation requirements.

| Samtools Version | Environment | Default install `/usr/local` | Custom install `/opt/install` |
| --------------- | ----------- | ---------------------------- | ----------------------------- |
| 1.10            | Ubuntu 22.04                               | ✅ | ✅ |
| 1.10            | R-Base 4.2.3 (*Debian* bookworm)           | ✅ | ✅ |
| 1.10            | Python 3.11.9 (*Debian* bookworm)          | ✅ | ✅ |
| 1.11            | Ubuntu 22.04                               | ✅ | ✅ |
| 1.11            | R-Base 4.2.3 (*Debian* bookworm)           | ✅ | ✅ |
| 1.11            | Python 3.11.9 (*Debian* bookworm)          | ✅ | ✅ |
| 1.12            | Ubuntu 22.04                               | ✅ | ✅ |
| 1.12            | R-Base 4.2.3 (*Debian* bookworm)           | ✅ | ✅ |
| 1.12            | Python 3.11.9 (*Debian* bookworm)          | ✅ | ✅ |
| 1.13            | Ubuntu 22.04                               | ✅ | ✅ |
| 1.13            | R-Base 4.2.3 (*Debian* bookworm)           | ✅ | ✅ |
| 1.14            | Ubuntu 22.04                               | ✅ | ✅ |
| 1.14            | R-Base 4.2.3 (*Debian* bookworm)           | ✅ | ✅ |
| 1.14            | Python 3.11.9 (*Debian* bookworm)          | ✅ | ✅ |
| 1.15            | Ubuntu 22.04                               | ✅ | ✅ |
| 1.15            | R-Base 4.2.3 (*Debian bookworm*)           | ✅ | ✅ |
| 1.15            | Python 3.11.9 (*Debian* bookworm)          | ✅ | ✅ |
| 1.16            | Ubuntu 22.04                               | ✅ | ✅ | 
| 1.16            | R-Base 4.2.3 (*Debian* bookworm)           | ✅ | ✅ |
| 1.16            | Python 3.11.9 (*Debian* bookworm)          | ✅ | ✅ |
| 1.17            | Ubuntu 22.04                               | ✅ | ✅ |
| 1.17            | R-Base 4.2.3 (*Debian* bookworm)           | ✅ | ✅ |
| 1.17            | Python 3.11.9 (*Debian* bookworm)          | ✅ | ✅ |
| 1.18            | Ubuntu 22.04                               | ✅ | ✅ |
| 1.18            | R-Base 4.2.3 (*Debian* bookworm)           | ✅ | ✅ |
| 1.18            | Python 3.11.9 (*Debian* bookworm)          | ✅ | ✅ |
| 1.19            | Ubuntu 22.04                               | ✅ | ✅ |
| 1.19            | R-Base 4.2.3 (*Debian* bookworm)           | ✅ | ✅ |
| 1.19            | Python 3.11.9 (*Debian* bookworm)          | ✅ | ✅ |
| 1.20            | Ubuntu 22.04                               | ✅ | ✅ |
| 1.20            | R-Base 4.2.3 (*Debian* bookworm)           | ✅ | ✅ |
| 1.20            | Python 3.11.9 (*Debian* bookworm)          | ✅ | ✅ |


## Development

To build
```bash
docker build -f docker/Dockerfile.ubuntu22.usr_local -t example:local .

# or to build with a specific version
VERSION=1.16
docker build -f docker/Dockerfile.ubuntu22.usr_local --build-arg HTSLIB_VERSION=$VERSION -t example:local .

```

To run
```bash
docker run -it --rm example:local bash

# or if wanting to bind mount the repo
docker run -it --rm -v $(pwd):/opt/repo example:local bash
```

To test
```bash
VERSION=1.16
docker run --rm example:local bash run_tests_samtools.sh $VERSION
docker run --rm example:local bash run_tests_htslib.sh $VERSION
```

## Cutting a release

To cut a release, update the version in the script and the README.md. This
repository uses [semantic versioning](https://semver.org/spec/v2.0.0.html).

Tests are automatically run in the GitLab CI pipeline.

Tags will automatically create releases on GitHub.