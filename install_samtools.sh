#!/bin/bash
# Description: Install samtools from source
# Usage: install_samtools.sh --help

set -euo pipefail

### GLOBAL VARIABLES ###
SETUP_DIR=""
INSTALL_DIR=""
PROGRAM_VERSION=""
IS_DEFAULT_INSTALL_DIR=1 # 1 == true, 0 == false
IS_HTSLIB_BUILD_DIR_SET=0 # 1 == true, 0 == false
HTSLIB_BUILD_DIR="na"

### CONSTANTS ###
PROGRAM_NAME="samtools"
TARBALL_SUFFIX=".tar.bz2"
DEFAULT_INSTALL_DIR="/usr/local"
SCRIPT_VERSION="1.0.0"
REQUIRED_PROGRAMS=(curl make autoconf autoheader gcc tar sed)
URL_TEMPLATE="https://github.com/samtools/samtools/releases/download/{}/samtools-{}.tar.bz2"

### FUNCTIONS ###

function print_info() {
  # Print a message in green to stderr
  local message="${1}"
  echo -e "\033[0;32m INFO: ${message}\033[0m" >&2
}

function print_warning() {
  # Print a message in yellow to stderr
  local message="${1}"
  echo -e "\033[0;33mWARNING: ${message}\033[0m" >&2
}

function print_error() {
  # Print a message in red to stderr
  local message="${1}"
  echo -e "\033[0;31mERROR: ${message}\033[0m" >&2
}

function print_usage() {
  echo "Usage: $0  <SAMTOOLS-VERSION> [--install-dir <INSTALL-DIR>] [--setup-dir <SETUP-DIR>] [ --installed-htslib-dir <HTSLIB-DIR>] [-h|--help] [--version]"
  echo ""
  echo "Installs ${PROGRAM_NAME:?} only, compiling it from source. This script does not install htslib."
  echo ""
  echo "Arguments:"
  echo "  SAMTOOLS-VERSION: The version of ${PROGRAM_NAME:?} to install e.g. 1.14"
  echo ""
  echo "Options:"
  echo "  --install-dir:    [Default: ${DEFAULT_INSTALL_DIR:?}] The directory to install "
  echo "                    ${PROGRAM_NAME:?} into, must be an absolute path and writable e.g. /opt/${PROGRAM_NAME:?}"
  echo "  --setup-dir:      [Default: (set by mktemp)] The directory to download and unpack "
  echo "                    the ${PROGRAM_NAME:?} tarball (must be writable). Deleted after installation."
  echo "  --installed-htslib-dir:   If specified, the path to an existing htslib installation is used to "
  echo "                    compile samtools. Otherwise, samtools will compile with the bundled "
  echo "                    htslib source (this WILL NOT install htslib)."
  echo "  -h, --help:       Print this message"
  echo "  --version:        Print the version of this script"
  echo ""
  echo "Required programs:"
  for program in "${REQUIRED_PROGRAMS[@]}"; do
    echo "  ${program}"
  done
}

function parse_args() {
  # Parse the command-line arguments, setting the INSTALL_DIR and PROGRAM_VERSION
  # global variables. If the arguments are invalid, print an error message and
  # exit.
  local user_setup_dir=""
  local install_dir=""
  local user_install_dir=""
  local program_version=""
  local is_default_install_dir=""
  local htslib_dir=""

  # Parse the command-line arguments
  while [ "$#" -gt 0 ]; do
    case "$1" in
      -h|--help)
        print_usage
        exit 0
        ;;
      --version)
        echo "${SCRIPT_VERSION}"
        exit 0
        ;;
      --setup-dir)
        if [ -z "${2}" ]; then
          print_error "Missing argument for --setup-dir"
          print_usage
          exit 1
        fi
        user_setup_dir="${2}"
        shift
        ;;
      --install-dir)
        if [ -z "${2}" ]; then
          print_error "Missing argument for --install-dir"
          print_usage
          exit 1
        fi
        user_install_dir="${2}"
        shift
        ;;
      --installed-htslib-dir)
        if [ -z "${2}" ]; then
          print_error "Missing argument for --installed-htslib-dir"
          print_usage
          exit 1
        fi
        htslib_dir="${2}"
        shift
        ;;
      *)
        if [ -z "${program_version}" ]; then
          program_version="${1}"
        else
          print_error "Invalid argument: ${1}"
          print_usage
          exit 1
        fi
        ;;
    esac
    shift
  done

  # Check that the required arguments are set
  if [ -z "${program_version}" ]; then
    print_error "Missing required arguments"
    print_usage
    exit 1
  fi

  # Normalize paths by removing trailing slashes
  local normalized_user_install_dir="${user_install_dir%/}"
  local normalized_default_dir="${DEFAULT_INSTALL_DIR%/}"

  # `install_dir`
  # Parse the user-specified install directory, or use the default if not set
  if [ -z "${user_install_dir}" ]; then
    install_dir="${DEFAULT_INSTALL_DIR}"
    is_default_install_dir=1
  elif [ "${normalized_user_install_dir}" = "${normalized_default_dir}" ]; then
    install_dir="${DEFAULT_INSTALL_DIR}"
    is_default_install_dir=1
  else
    install_dir="${user_install_dir}"
    is_default_install_dir=0
  fi
  
  # `install_dir`
  # Check that the install directory is valid: it is a directory and can be
  # written to and is an absolute path.
  if [ ! -d "${install_dir}" ]; then
    print_error "Invalid install directory: ${install_dir}"
    exit 1
  fi
  if [ ! -w "${install_dir}" ]; then
    print_error "Cannot write to install directory: ${install_dir}"
    exit 1
  fi
  if [[ ! "${install_dir}" = /* ]]; then
    print_error "Install directory must be an absolute path: ${install_dir}"
    exit 1
  fi

  # `setup_dir`
  # Parse the user-specified setup directory, or use a temporary directory if not set
  if [ -z "${user_setup_dir}" ]; then
    # append the program version to the setup directory and a random string
    setup_dir=$(mktemp -d -t "${PROGRAM_NAME:?}-${program_version:?}-XXXXX")
  else
    setup_dir="${user_setup_dir}"
  fi

  # `setup_dir`
  # Check that the setup directory is valid: it is a directory and can be
  # written to.
  if [ ! -d "${setup_dir}" ]; then
    print_error "Invalid setup directory: ${setup_dir}"
    exit 1
  fi
  if [ ! -w "${setup_dir}" ]; then
    print_error "Cannot write to setup directory: ${setup_dir}"
    exit 1
  fi

  # `program_version`
  # Check that the program version is valid format e.g. 1.14 otherwise throw a warning
  if [[ ! "${program_version}" =~ ^[0-9]+\.[0-9]+$ ]]; then
    print_warning "Unexpected ${PROGRAM_NAME:?} version: ${program_version}, expected format is X.Y"
  fi

  # `htslib_dir`
  # Check that the htslib directory is valid: it is a directory
  if [[ -n "${htslib_dir}" ]] && [[ ! -d "${htslib_dir}" ]]; then
    print_error "Invalid htslib directory: ${htslib_dir}"
    exit 1
  fi

  # Set the global variables
  INSTALL_DIR="${install_dir}"
  print_info "Setting INSTALL_DIR=${INSTALL_DIR}"
  SETUP_DIR="${setup_dir}"
  print_info "Setting SETUP_DIR=${SETUP_DIR}"
  PROGRAM_VERSION="${program_version}"
  print_info "Setting PROGRAM_VERSION=${PROGRAM_VERSION}"
  IS_DEFAULT_INSTALL_DIR="${is_default_install_dir}"
  print_info "Setting IS_DEFAULT_INSTALL_DIR=${IS_DEFAULT_INSTALL_DIR} (1 == true, 0 == false)"
  if [[ -n "${htslib_dir}" ]]; then
    IS_HTSLIB_BUILD_DIR_SET=1
    HTSLIB_BUILD_DIR="${htslib_dir}"
    print_info "Setting HTSLIB_BUILD_DIR=${HTSLIB_BUILD_DIR}"
  fi

}

function assert_programs_exists() {
  # Check that all required programs are installed by iterating over an array of
  # program names and accumulating the missing ones in a list. If any are
  # missing, print an error message and exit.
  local required_programs=("$@")
  local missing_programs=()
  for program in "${required_programs[@]}"; do
    if ! command -v "${program}" &> /dev/null; then
      missing_programs+=("${program}")
    fi
  done

  if [ "${#missing_programs[@]}" -gt 0 ]; then
    print_error "Cannot complete installation because the following programs are missing: ${missing_programs[*]}"
    print_error "Please install them and try again."
    exit 1
  fi
}

function format_url() {
  # Format the URL for downloading the program download asset by replacing the
  # placeholders with the version number.
  local version="${1}"
  local template_url="${URL_TEMPLATE}"
  local url=$(echo "$template_url" | sed "s|{}|$version|g")
  echo "$url"
}

function test_url() {
  local url="${1}"
  print_info "Checking URL: ${url}"
  
  # Use curl to make a HEAD request, following redirects, and capture both the final URL and status code
  local result=$(curl -sI -o /dev/null -w "%{http_code}" -L "$url")
  local status_code=$(echo "$result")

  case $status_code in
    200)
      print_info "URL is valid."
      return 0
      ;;
    301|302|307|308)
      print_warning "URL is a redirect."
      return 0
      ;;
    *)
      print_error "Unexpected status code ${status_code} for URL: ${url}"
      return 1
      ;;
  esac
}

function download_url() {
  # Download the URL to a file using curl, retrying up to 5 times
  local url="${1}"
  local file="${2}"
  print_info "Downloading from ${url} to ${file}"
  curl -o "$file" -sSL --retry 5 --location "$url"
}

function unpack_tarball() {
  # Unpack the tarball into the specified directory
  local tarball="${1}"
  local unpack_dir="${2}"

  print_info "Unpacking ${tarball:?} to ${unpack_dir:?}"
  mkdir -p "${unpack_dir:?}"
  tar --strip-components 1 -C "${unpack_dir:?}" -jxf "${tarball:?}"
}

function clean_path_strings() {
  # Remove any trailing colons from the path strings
  local path="${1}"
  # `%%:` is fancy syntax to remove a trailing colon from the variable, if present.
  path="${path%%:}"
  echo "${path}"
}

function get_cpu_count() {
    # Get the number of CPUs available on the system, with an optional threshold
    local threshold=${1:-0}  # Default threshold is 0 (no limit)
    local cpu_count=1  # Default to 1 CPU if unable to determine

    # Try to get CPU count from lscpu
    if command -v lscpu > /dev/null 2>&1; then
        cpu_count=$(lscpu -p | egrep -v '^#' | sort -u -t, -k 2,4 | wc -l)
        print_info "CPU count determined using lscpu: ${cpu_count}"
    # If lscpu fails, try nproc
    elif command -v nproc > /dev/null 2>&1; then
        cpu_count=$(nproc)
        print_info "CPU count determined using nproc: ${cpu_count}"
    # If both fail, try counting from /proc/cpuinfo
    elif [ -f /proc/cpuinfo ]; then
        cpu_count=$(grep -c ^processor /proc/cpuinfo)
        print_info "CPU count determined from /proc/cpuinfo: ${cpu_count}"
    else
        print_warning "Unable to determine CPU count. Defaulting to 1 CPU."
    fi

    # Check if cpu_count is a number
    if ! [[ "$cpu_count" =~ ^[0-9]+$ ]]; then
        print_error "Invalid CPU count obtained: ${cpu_count}. Defaulting to 1 CPU."
        cpu_count=1
    fi

    # If threshold is set and cpu_count exceeds it, return threshold
    if [ "$threshold" -gt 0 ] && [ "$cpu_count" -gt "$threshold" ]; then
        print_warning "CPU count (${cpu_count}) exceeds threshold (${threshold}). Using threshold value."
        echo "$threshold"
    else
        echo "$cpu_count"
    fi
}

function install() {
  # Args
  local install_dir="${1}"
  local work_dir="${2}"
  local is_default_install_dir="${3}"
  local htslib_build_dir="${4}" # This value could be "na" if not set -- rely on is_htslib_build_dir_set to determine if it's set
  local is_htslib_build_dir_set="${5}"  # 1 == true, 0 == false
  
  # Local variables
  local cmd=""
  local cmd_non_default_arg__prefix="--prefix=${install_dir:?}"
  local cmd_non_default_arg__cppflags="CPPFLAGS=\"-I${install_dir:?}/include\""
  local cmd_non_default_arg__ldflags="LDFLAGS=\"-L${install_dir:?}/lib -Wl,-R${install_dir:?}/lib\""
  local cmd_optional_arg__htslib="--with-htslib=${htslib_build_dir:?}"
  local cpu_count=$(get_cpu_count 6)



  print_info "Changing current directory to ${work_dir:?}"
  cd "${work_dir:?}"
  
  print_info "Running autoheader..."
  autoheader

  print_info "Running autoconf... (generating configure script)"
  autoconf -Wno-syntax

  print_info "Check that configure script exists"
  test -f configure || {
    print_error "configure script not found in ${work_dir:?}"
    exit 1
  }

  print_info "Running configure script"
  cmd="./configure"
  if [ "${is_default_install_dir:?}" -eq 0 ]; then
    cmd="${cmd:?} ${cmd_non_default_arg__prefix:?}"
    cmd="${cmd:?} ${cmd_non_default_arg__cppflags:?}"
    cmd="${cmd:?} ${cmd_non_default_arg__ldflags:?}"
  fi
  if [ "${is_htslib_build_dir_set:?}" -eq 1 ]; then
    cmd="${cmd:?} ${cmd_optional_arg__htslib:?}"
  fi

  print_info "Running configure with command: '${cmd:?}'"
  eval "${cmd:?}"

  print_info "Running make"
  make -j"${cpu_count:?}"

  print_info "Running make install"
  make install
}

function clean_up() {
  local setup_dir="${1}"
  local tarball="${2}"
  print_info "Cleaning up temporary files"
  rm -rf "${setup_dir:?}"
  rm -f "${tarball:?}"
}

function main() {
  # Positional arguments
  local install_dir="${1}"
  local setup_dir="${2}"
  local program_version="${3}"
  local is_default_install_dir="${4}"
  local htslib_build_dir="${5}"
  local is_htslib_build_dir_set="${6}"

  # Local constants
  local tarball="${install_dir:?}/${PROGRAM_NAME:?}-${program_version:?}.${TARBALL_SUFFIX:?}"
  local unpack_dir="${setup_dir:?}/${PROGRAM_NAME:?}"
  local url=$(format_url "${program_version:?}")

  # Download and unpack the tarball
  test_url "${url:?}"
  download_url "${url:?}" "${tarball:?}"
  unpack_tarball "${tarball:?}" "${unpack_dir:?}"

  # Update environment paths for shared libraries and headers (we disable the
  # set -u check for this, as the *_PATH variables may unset)
  print_info "Updating environment variables: LD_LIBRARY_PATH, PATH, MANPATH"
  set +u
  LD_LIBRARY_PATH="${install_dir:?}/lib:${LD_LIBRARY_PATH}"
  PATH="${install_dir:?}/bin:$PATH"
  MANPATH="${install_dir:?}/man:${install_dir:?}/share/man:$MANPATH"
  set -u
  export LD_LIBRARY_PATH=$(clean_path_strings "${LD_LIBRARY_PATH:?}")
  export PATH=$(clean_path_strings "${PATH:?}")
  export MANPATH=$(clean_path_strings "${MANPATH:?}")

  # Run the installation
  install "${install_dir:?}" "${unpack_dir:?}" "${is_default_install_dir:?}" "${htslib_build_dir:?}" "${is_htslib_build_dir_set:?}"

  # Clean up the installation directory
  clean_up "${setup_dir:?}" "${tarball:?}"

  print_info "Installation complete."
}

### RUNTIME ###

parse_args "$@"
assert_programs_exists "${REQUIRED_PROGRAMS[@]}"
# Save the current working directory and restore it when the script exits
cwd=$(pwd)
trap 'exit_code=$?; cd "$cwd"; exit $exit_code' EXIT SIGHUP SIGINT SIGTERM SIGQUIT
main "${INSTALL_DIR:?}" "${SETUP_DIR:?}" "${PROGRAM_VERSION:?}" "${IS_DEFAULT_INSTALL_DIR:?}" "${HTSLIB_BUILD_DIR:?}" "${IS_HTSLIB_BUILD_DIR_SET:?}"

# End of script