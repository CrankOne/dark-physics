# A Geant4 Monte-Carlo Modules for Simulating Dark Sector Physics

This repository is dedicated to simulating elementary particle physics from
so called "dark sector": hypothetical dark matter mediators, axion-like
particles, light thermal dark matter and so on.

# Models Included in Repo

The purpose of this repo is, of course not to cover every aspect of the dark
sector theory, but rather to help people to gain some numerical predictions for
relevant phenomenology concerning few well-elaborated models that appeared up
to the time.

## Dark Photon Particle and Corresponding Production Processes

"Dark photon" is massive vector-like mediator particle supposingly
participating in electromagnetic processes to which it is weakly coupled. We
expect it to appear mostly in high-energetic processes. The numerical
estimation for the production cross-section process relies on the
Weitzacker-Williams approximation.

# Project Structure

The primary binary artifact is the shared library called `libdphmc.so`.
Following shared library layout convention that is standard for UNIX systems,
the versioning symlinks will be created after build and installation stages.
All the relevant code is located in top level `include/` and `src/` directories
as a common practice for single-library projects.

Optionally, testing toolkit may be built. It provides some basic testing
fixtures for testing and debugging the procedures (unit tests are provided as
well). See corresponding readme file(s) within `tests/` directory for
instructive descriptions.

## Notes on File Structure and Code Guidelines

* Following to original Geant4 conventions, the strict correspondance between
header/source files has to be respected: one class -- one header.
* However, as it is still an extension library, custom classes has to be
denoted either by naming prefix (in Geant4 it is `G4`), or be in the namespace.
Namespace was choosenas it is natural for modern C++ and more laconic then name
suffixes.
* The code within the package is written on C++ and pure C (to decrease the
build time for some numerical procedures that does not require much expressive
power). For C++ the common namespace is `dphmc`. For C sources all the public
declarations must be prefixed with `dphmc_`.
* As the main purpose of this library is to extend the Gean4 process lib, all the
Geant4-relevant entities have to be first-level file entries in `include/`
and `src/`
* We also tend to minimize dependencies list, besides standard Geant4
installtion to expand the use case and portability.

