# svFSI Testing Guide

[Integration testing](https://en.wikipedia.org/wiki/Integration_testing) is an essential part of software development. Essentially, assume that every untested line in `svFSIplus` is not working correctly. Below is a quick guide on how to run and add integration tests for `svFSI`.

## Prerequisites
There are two things you need to do before you can run a test case: Build `svFSI` and install `Git LFS` to download the test cases.

### Build svFSI
Follow the build instructions outlined [here](https://simvascular.github.io/svFSIplus/index.html#autotoc_md52). Importantly, to automatically run test cases with `pytest` (see below), you need to build `svFSI` in the folder
```
./build
``` 
in the repository root.

### Install Git LFS
You need to install [`Git LFS`](https://git-lfs.com/) (*Large File Storage*) to run any test, which we use to track large files. Tracking large files with `Git` can significantly add to the repository size. These large files include meshes and boundary conditions for all test cases. They are stored on `GitHub`, but the files themselves just contain a hash or Object ID that is tracked with `Git`. All file extensions currently tracked with Git LFS are listed under [in this file](../.gitattributes).

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

## Running tests with pytest
You can run an individual test by navigating to the `cases/<physics>/<test>` folder you want to run and execute `svFSIplis` witt the `svFSI.xml` input file as an argument. A more elegant way, e.g., to run a whole group of tests, is using `pytest`. By default, it will run all tests defined in the `test_.py` files in the `./tests` folder. Tests and input files in `./tests/cases` are grouped by physics type, e.g., `struct`, `fluid`, or `fsi` (using the naming convention from `EquationType`). Here are a couple of useful `Pytest` commands:

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
We expect that new code is fully covered with at least one integration test. We also strive to increase our coverage of existing code. You can have a look at our current code coverage [with Codecov](https://codecov.io/github/SimVascular/svFSIplus). It analyzes every pull request and checks the change of coverage (ideally increasing) and if any non-covered lines have been modified. We avoid modifying untested lines of code as without a test, there is no guarante that the code will still do the same thing as before.

## Creating new tests


### fixtures