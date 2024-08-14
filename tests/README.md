@page testing Testing Guide

[TOC]

[Integration testing](https://en.wikipedia.org/wiki/Integration_testing) is an essential part of software development. It is performed when integrating code changes into the main development branch to verify that the code works as expected. Below is a quick guide on how to run and add integration tests for `svFSI`.


## Prerequisites
There are two things you need to do before you can run a test case: Build `svFSI` and install `Git LFS` to download the test cases. To run certain test cases, you also need `svZeroDSolver` (see below).

### Build svFSI
Follow the build instructions outlined [here](https://simvascular.github.io/svFSIplus/index.html#autotoc_md52). Importantly, to automatically run test cases with `pytest` (see below), you need to build `svFSI` in the folder
```
./build
``` 
in the repository root.


### Install Git LFS
You need to install `Git LFS` ([*Large File Storage*](https://git-lfs.com/)) to run any test, which we use to track large files. Tracking large files with `Git` can significantly add to the repository size. These large files include meshes and boundary conditions for all test cases. They are stored on `GitHub`, but the files themselves just contain a hash or Object ID that is tracked with `Git`. All file extensions currently tracked with `Git LFS` are listed under [in this file](../.gitattributes).

When using `Git LFS` for the first time, you need to follow these simple steps:
1. Install on your platform by following [this guide](https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage).
2. Initialize in your `svFSIplus` repository with
    ```
    git lfs install
    ```
3. Download all large files with
    ```
    git lfs pull
    ```
After performing these steps once, you never need to worry about Git LFS again. All large files are handled automatically during all Git operations, like `push`, `pull`, or `commit`.


### Build svZeroDSolver
Some test cases require svZeroDSolver built in the svFSIplus directory.
Importantly, to automatically run test cases with `pytest` (see below), you need to build `svZeroDSolver` in the folder
```
./svZeroDSolver/build
``` 
in the repository root.

To do so, you can run the following in the svFSIplus repository root:
```
git clone https://github.com/SimVascular/svZeroDSolver.git
cd svZeroDSolver
mkdir build
cd build
cmake ..
make -j2
``` 


## Running tests with pytest
You can run an individual test by navigating to the `./tests/cases/<physics>/<test>` folder you want to run and execute `svFSIplus` with the `svFSI.xml` input file as an argument. A more elegant way, e.g., to run a whole group of tests, is using [`pytest`](https://docs.pytest.org/). By default, it will run all tests defined in the `test_*.py` files in the [./tests](https://github.com/SimVascular/svFSIplus/tree/main/tests) folder. Tests and input files in [./tests/cases](https://github.com/SimVascular/svFSIplus/tree/main/tests/cases) are grouped by physics type, e.g., [struct](https://github.com/SimVascular/svFSIplus/tree/main/tests/cases/struct), [fluid](https://github.com/SimVascular/svFSIplus/tree/main/tests/cases/fluid), or [fsi](https://github.com/SimVascular/svFSIplus/tree/main/tests/cases/fsi) (using the naming convention from `EquationType`). Here are a couple of useful `Pytest` commands:

- Run only tests matching a pattern (can be physics or test case name):
    ```
    pytest -k ustruct
    pytest -k block_compression
    ```
- List individual test cases that were run:
    ```
    pytest -v
    ```
- Print the simulation output of all tests:
    ```
    pytest -rP
    ```

For more options, simply call `pytest -h`.

## Code coverage
We expect that new code is fully covered with at least one integration test. We also strive to increase our coverage of existing code. You can have a look at our current code coverage [with Codecov](https://codecov.io/github/SimVascular/svFSIplus). It analyzes every pull request and checks the change of coverage (ideally increasing) and if any non-covered lines have been modified. We avoid modifying untested lines of codeas there is no guarante that the code will still do the same thing as before.

## Create a new test
Here are some steps you can follow to create a new test for the code you implemented. This will satisfy the coverage requirement (see above) and help other people who want to run your code. A test case is a great way to show what your code can do! Ideally, you do this early in your development. Then you can keep running your test case as you are refractoring and optimizing your code.

1. Create an **example** (mesh and input file) that showcases what your code can do. If it doesn't cover everything, consider adding other tests.
2. **Verify** the results: Since an analytical solution would be ideal but is rarely available, find other ways to ensure that your code is doing the right thing (reference solutions from other codes, manufactured solutions, convergence analysis, ...).
3. Crank up the **tolerances**: Be as strict in your linear and nonlinear solver as you can be with the test still converging. This will avoid problems when running the test on different machines or with multiple processors.
4. **Reduce** the computational effort: Try to make the test run in a few seconds by coarsening spatial and temporal discretization (but keeping the core function of the test).
5. Set the **maximum number** of nonlinear iterations to the one where it currently achieves the set tolerance. This way, the test will fail if the linearization gets broken in the future (but the test still slowly converges to the correct solution).
6. **Test the test**: Does it fail if you change parts of your input file that the test should be sensitive to (e.g., material parameters if you implemented a new solid material)?
7. Put your files in the appropriate **folder** under `./tests/cases` and append the test to the Python `test_*.py` file. If you created a new physics type, create new ones for both.
8. Check that the test is **executed** correctly in `GitHub Actions` when opening your pull request. You should see in your coverage report that your new code is covered.

If you want to parameterize values in your test case, you can use `@pytest` [fixtures](https://docs.pytest.org/en/6.2.x/fixture.html). We currently use them to automatically loop different numbers of processors, meshes, or input files.