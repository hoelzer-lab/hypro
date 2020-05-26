**Date**    27/05/2020
**Author**    Maximilian Arlt

# How to write a conda package
This protocol supports you in writing your first conda package from scratch. 

## General Remarks

1. The first important note: Conda packaging is very-well documented. I tried to wrap some of the most important and handy links in this document. Those shall further support the explanations that are given in the respective section.
2. As you will quickly experience yourself, one conda package is not like the other. There is a lot you can and have to customize and decide, depending on the code source, programming language and composition of the program you want to package. 
For details, see [building the meta.yaml file](#Building-the-meta.yaml-file). Here, this becomes more clear. Don't expect everything you need to be written in this protocol, that's impossible!
3. The protocol was written from the experience in packaging [HyPro](https://github.com/hoelzer-lab/hypro), a python-based tool performing homology-based searches for protein sequences to extend hypothetical protein annotations of Prokka. All files necessary for building the conda package can be revisited [here](https://github.com/hoelzer-lab/hypro) as I found many examples in the conda documentation rather advanced.

## Main requirements
To build a conda package, you first have to write a conda recipe. Here, conda finds all the information necessary to construct the package. The recipe is composed of:

* [**meta.yaml**](#building-the-metayaml-file) the script containing all the metadata 
* [**build.sh** or **bld.bat**](#buildsh) - a build script installing all files for the package. Shell script for macOS and Linux, executed using `bash` command. .bat for Windows, executed using `cmd`. In case of a Python tool, this executable calls the `setup.py`.
* **run_test.[py,pl,sh,bat]** - a Python test script testing proper functioning of the program in the package. Runs automatically when defined.

* Optional patches and resources, not introduced in this protocol.
* If you write a recipe for a Python tool, the build-executable simply calls your `setup.py`. Things to know about Python packaging are summarized in [this chaper](#Defining-the-`setup.py`-for-conda)

It is recommended to put those files together in a folder "conda-recipe" or something like that. 

**Note:**
A basic package can already be constructed defining only the meta.yaml. However, to meet requirements of certain repositories like the Anaconda Cloud or Bioconda, you have to meet further criterions. In the following sections, this will become more clear.
For more remarks on constructing recipes, visit the [conda recipe documentation](https://docs.conda.io/projects/conda-build/en/latest/concepts/recipe.html).

## Walking through the build process
### **Build the package**

Once your recipe is written it will be used to form the conda package. For that, you need the `conda-build` command. 
Additionally, you may install [`conda-verify`](https://github.com/conda/conda-verify). This is a handy tool that checks for the right path configuration and points on conflicts or other malfunctions in your package. It is a nice extension to the basic `conda-build` command and starts automatically after `conda-build` finished the package creation.

1. Set up the conda environment you want to build the package in.
Make sure to set your channel configuration properly. It is recommended to set conda-forge, bioconda and the defaults channel for most purposes:

```
conda create -n my_build_env --python==<version>
conda activate my_build_env
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
`conda install` will search these channels whenever you ask for installing a tool (first conda-forge, then bioconda, then defaults). 

1. Install `conda-build`.
2. Install `conda-verify`.
3. Run `conda-build` to create a package from your directory.

```
conda install conda-build
conda install conda-verify
conda-build <path-to-recipe-dir>
```

Optional, you may set the `--prefix-length` parameter (requires an INT). Otherwise `conda-build` gives a loooong placeholder to certain directories created during building process.

Generally, `conda-build` creates two environments for different purposes: build and test. 

`conda-build` performs the following steps ([source](https://docs.conda.io/projects/conda-build/en/latest/concepts/recipehtml#conda-build-process)):

1. Reads the metadata (meta.yaml). On writing the meta.yaml, see the details [here](#building-the-metayaml-file)
2. Downloads the source into a cache.
3. Extracts the source into the source directory.
4. Applies any patches.
5. Re-evaluates the metadata, if source is necessary to fill any metadata values.
6. Creates a build environment and then installs the build dependencies there. Those are dependencies necessary to build the package (e.g. `git` or `pip` for source loading, python for installing a python package)
7. Runs the build script. The current working directory is the source directory with environment variables set. The build script installs into the build environment. In case of a python tool, this simply calls the `setup.py` in the source directory.
8. Performs some necessary post-processing steps, such as shebang and rpath.
9. Creates a conda package containing all the files in the build environment that are new from step 5, along with the necessary conda package metadata.
10. Tests the new conda package if the recipe includes tests (required for upload to bioconda later on):
    1. Deletes the build environment and source directory to ensure that the new conda package does not inadvertantly depend on artifacts not included in the package.
    2. Creates a test environment with the package and its dependencies.
    3. Runs the test scripts.

If the build process finishes correctly, a compressed tarball file had been created (looks something like ``/home/<path-to-conda>/conda-bld/<os>/<mypkg>-<version>-<interpreter>.tar.bz2``, with 'os' being e.g. "linux-64"). This is your conda package. `conda-build` prints the exact path. Keep the path, you might need it for local testing and upload.

Clean up your conda-bld directory from failed packages - the build environments will still be present - delete them!

**Note:**
1. Most likely, your first `conda-build` command will fail. Maybe the second also. Keep on trying and read the error messages carefully, they might help. When your package was built, under `/home/path-to-conda/conda-bld/` there will be as many broken environments as tries you had. They appear as directories, named like your conda package name plus a long digit code. Make sure to delete all those directories, **before** running your local install test. In my case, those impaired the installation process and led to some conflicts.
2. `conda-build` defines many environment variables. A list of those is provided [here](https://docs.conda.io/projects/conda-build/en/latest/user-guide/environment-variables.html#environment-variables-set-during-the-build-process). One of those set wrongly is a common source for failing the build process. You should know or define where they are pointing to. More on that in the example in sections [meta.yaml](#building-the-metayaml-file) and [`setup.py`](#Defining-the-`setup.py`-for-conda).

### **Test package locally**
1. Create a new conda environment:
```
conda create -n test_env
conda activate test_env
```
In this environment you can install and test your conda package locally, before loading it to the cloud and opening it to everybody.
2. Install the package locally
The conda documentation recommends the command:
```
conda install --use-local <my-pkg-name>
```
However, this did not work for me since conda was not able to find my package.
I suggest to specify the local channel:

```
conda install -c file:///home/<path-to-conda>/envs/<env-name>/conda-bld/<os>/ <my-pkg-name>
```

The installer will ask you to install the dependencies for your package like during any other conda tool installation process. 
1. Check, the requirements. If all are there, continue. If not, you may check your meta.yaml. Most likely your `conda-build` process would break anyway.
2. After the installation finished, you are free to test your conda environment. Yippi!

### **Upload package**

#### Upload to Anaconda Cloud

This is the easiest way of distributing your package. An Upload to anaconda makes your package generally available under a channel you define by yourself. For convenience, this can be your anaconda username.

1. Create an account at [anaconda](https://anaconda.org/)
2. Login: In the dashboard you get an overview on packages, projects, notes and so on. For more information click on [upload a package](https://docs.anaconda.com/anaconda-cloud/user-guide/tasks/work-with-packages/#uploading-packages)

Follow the instructions provided.

#### Upload to Bioconda Recipe Repository

Bioconda provides a documentation for contributing to Bioconda. You can find it [here](https://bioconda.github.io/contributor/index.html). For contributing to bioconda, you need a github account. If you do not have one already, now is the time.

0. Be as sure as you can that your package and tool functions correctly. You can never withdraw from whats in the cloud! 
1. Read the guidlines of bioconda on [this page](https://bioconda.github.io/contributor/guidelines.html). If something is missing, you might include it in your recipe.
  Most important here is including a sha256 checksum from a stable URL. As stable bioconda considers a packed tarball file from a github release. On this, you can calculate the checksum using:
```
wget -O- <URL> | shasum -a 256
```
There are also other sources considered as stable, for more information, consider the guidlines mentioned before. 
3. Bioconda recipes are stored in a github repository called "bioconda-recipes". For contribution, you are obliged to fork the bioconda-recipes: Login to your github account and click [here](https://github.com/bioconda/bioconda-recipes/fork) to create a remote fork which is a repository in your account.
4. Now, that you have a copy of bioconda-recipes, create a local clone on your computer running:
   ```
      git clone https://github.com/<USERNAME>/bioconda-recipes.git
   ```
5. Set the bioconda-recipes repo as upstream branch to be tracked by your local branch. Doing so, your local branch keeps track with the actual bioconda repo. More about upstream branches: [here](https://devconnected.com/how-to-set-upstream-branch-on-git/).

```
cd bioconda-recipes
git remote add upstream https://github.com/bioconda/bioconda-recipes.git
```

6. After pulling the latest changes in the bioconda repo, create a new branch with an appropriate name.
```
# Make sure our master is up to date with Bioconda
git checkout master
git pull upstream master
git push origin master

# Create and checkout a new branch for our work, replace <packagename> accordingly
git checkout -b <packagename>_recipe
```
7. Now, you may add your recipe to the branch. In my case, this was simply the ``meta.yaml``. For that, you need to create a directory in the ``recipes`` folder of your local branch:
```
cd ./recipes
mkdir -p <package_name>
```
Also have a look into the other packages to get familiar with the structure. Some are pretty simple, others complex. I recommend to look into segemehl's recipe directory, find it in the recipes folder.

8. Now, add, commit, and push your changes:
```
git checkout <packagename>_recipe
git add <package_name>
git commit -m 'Meaningful message'
git push origin <packagename>_recipe
```
You may always check the differences and your active branch using ``git status`` or ``git branch``.

9. You are almost there! Now its time to request merging.
  - Go to your branch on the github webpage. There should be something written like "branch is X commits ahead of bioconda-recipes:master". On the right-hand side, click on ``Pull request``.
  - When the request is successful, bioconda will automatically build your package from your recipe and makes it available over the bioconda channel.
  - To reach a successful pull request, you need to pass the bioconda test bot and at least one verified manual review by another bioconda user. 
  **Bioconda test bot**:
  It will perform three tests: linting, test-linux and test-macos. Most likely, those tests will fail the first time to try. If you are lucky, the manual review by another user may help you passing those tests.
  However, if you want to test your package locally, consider the [local testings](https://bioconda.github.io/contributor/building-locally.html#using-the-bootstrap-method) of bioconda.
  Those steps might become rather tricky and to the point as I write this protocol, HyPro did not pass the test-linux step. I suggest to get help from experienced bioconda reviewers, maybe the person that checked your recipe if the testing does not give sufficient info about the problem.

## Packaging in Practice: An Example

In this section I want to make the recipe I created transparent and comment on each file I used. This shall give you more info on how things are ment in the online documentary of conda. Note that I only used a few features and settings possible for setting up a package, especially in the meta.yaml.

### **Building the meta.yaml file**

This is the heart piece of your recipe. Please find the general structure and a first example [here](https://docs.conda.io/projects/conda-build/en/latest/resources/define-metadata.html). Further this website provides **all options** on the meta.yaml, sorted per section. Most sections have many more options available than the ones mentioned. Please visit the website for these information.
Please always mind the exact format - code indentations, whitespaces etc. might all have a meaning so keep the format straight!
For comparison, feel free to visit the [HyPro meta.yaml](https://github.com/hoelzer-lab/hypro/blob/master/conda_recipe/meta.yaml).

#### **The Header**
```
{% set name = "hypro" %}
{% set version = "0.1" %}
```

Those are the first two lines - here, you specify the package name and package version. Mind the format!

#### **The Package Section**
```
package:
  name: hypro
  version: {{ version }}
```
Specify the package information here. Set the version number in double-quotes, so it cannot be interpretet as float, whatever number you choose.

#### **The Source Section**
```
  url: https://github.com/hoelzer-lab/hypro/archive/v{{ version }}.tar.gz
  sha256: 55510113190bc3c1154e7a67f2046970dd34b4388e88e81936fa73737bc3879e
```

Specify the code source of your package. `conda-build` will load all tool binaries from here and try to wrap them in the package. You may choose your source github repository; in this case the code is loaded from a githup repository. For providing a stable source, prepare a github release first and enter the url to the packaed tarball. Other sources are: local, hg, svn. Loading multiple source files is also possible - more details in the doc mentioned above.
The folder option specifies the destination path. `conda-build` will create this within the source-path and save the source information here. Default is the `<packagename>-<version>`. 

#### **The Build Section**
```
build:
    number: 0  
    noarch: python     
    script: {{ PYTHON }} hypro-{{ version }}/`setup.py` install --single-version-externally-managed --record=record.txt
```

Define all build information in here. The build number should be incremented for every build command. If using 'script', it replaces your `build.sh` and `bld.bat`. You should remove those files in this case. Noarch should be used for pure python packages.
Many more options like defining RPATHs or entry points of python scripts can be found on the mentioned website.

#### **The Requirements Section**

```
  build:
    - git
    - pip==19.3 
    - python==3.7.6
    - setuptools
```

Specify all requirements that are necessary for building the package. This means tools required loading the source code (git) and installing the tool binaries (python, setuptools). You may specify the exact version of the programs as min/max required (>/<=) or an specific version (==).

```
  host:         
    - python==3.7.6
    - pip==19.3 
    - prokka>=1.14.6
```

```
  run:
    - python==3.7.6
    - mygene==3.1.0
    - pandas==0.25.2
    - mmseqs2==10.6d92c
    - prokka>=1.14.6
```

Define the tools for executing your program. Conda will load them in advance, as happened with the build requirements. E.g. mmseqs2 and prokka. In case of python, this includes all non-standard libraries (standard are e.g. sys, os, argparse). In my case, the python script utilizes pandas and the mygene python API. The libraries can then be found by setuptools when executing your `setup.py` script.
Note, that the host section was skipped for the HyPro package - it was not required.

#### **The About Section**

```
about:
  home: https://github.com/hoelzer-lab/hypro.git
  summary: Extend hypothetical prokka protein annotations using additional homology searches against larger databases
  license: GPL
```

Here you add specifying information about the package. Where is the package maintained and under which license it was published, in this case.

Note: If you want to make your package available on different systems, you should make use of selectors. `conda-build` notices on which system it is run, so the selectors help to decide which information to use and which not. This is a convinient option I recommend to you.

### `Build.sh`

The `build.sh` can be a simple bash script containing only the one line you saw in [the build section](#the-build-section). However, this did not work for me because `conda-build` did not pass the environment variables to the bash correctly. Maybe you figure out how this works (pull requests are always welcome). 


### `Bld.bat`
This is the executable build script for Windows. In the most basic form it should contain exactly the following command:

```
"%PYTHON%" `setup.py` install
if errorlevel 1 exit 1
```
See the sample recipes and the conda doc for more information.

**Note:**
The build scripts can also contain a little more code and options - see the scripts provided with the [conda example recipes](https://docs.conda.io/projects/conda-build/en/latest/user-guide/recipes/sample-recipes.html). 

### Defining the `setup.py` for conda

Writing this script is the basis of formulating a python package. The following links contain more information about this:
[Writing a python package](https://packaging.python.org/tutorials/packaging-projects/)

[Writing the `setup.py`](https://docs.python.org/3/distutils/setupscript.html])
Here, you see the code of HyPro's `setup.py`:

```bash
from setuptools import setup, find_packages

requirements=['pandas==0.25.2','mygene==3.1.0'] 

long_description="HyPro uses mmseqs2 to search for homology of sequence segments that have been annotated by Prokka as hypothetical proteins. These sequences are searched in a nucleotide database to obtain a more specific Prokka annotation. In this way, HyPro updates the Prokka output, as the program is designed to be easily integrated into custom analysis pipelines."
setup(
    name="HyPro",
    version='0.1',
    author='Maximilian Arlt and Martin Hoelzer',
    author_email='hoelzer.martin@gmail.com',
    license='GPLv3',
    description='Extension for homology searches of hypothetical proteins to enhance Prokka annotations',
    install_requires=requirements,
    scripts=['hypro-0.1/scripts/hypro.py', 'hypro-0.1/scripts/mmseqs2.sh'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/hoelzer-lab/hypro.git',       
    packages=find_packages(), 
    classifiers=['Programming Language :: Python :: 3.7.6',  
    'Operating System :: Unix',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ]
)
```

The script firstly imports *setup* and *find_packages* from the setuptools module. Specifying the arguments for the setup function is key for successful packaging. The whole scope of options can be revisited in [Writing the `setup.py`](https://docs.python.org/3/distutils/setupscript.html). Here, we concentrate on specified ones:

`name` - *The name of your package*

`version` - *Version Number. You may stick to the [python version number conventions](https://www.python.org/dev/peps/pep-0396/).*

`author` - *Author Names*

`author_email` - *Autor email adresse*

`license` - *License to publish your tool*

`description` - *Short description of the tool - give the user an idea of what it's good for!*

`install_requires` - *List of required modules to load. For long lists, it is convenient to define a list variable above the setup command (see example above).*

`scripts` - *Define your executables here. In this case, hypro is called using the hypro.py script command line interface which requires the mmseqs2.sh as executable. Define relative paths, conda-build wioll extend them to ${SOURCE_DIR}/rpath_to_script. When finished successfully, the scripts will be available in the /bin of te conda environment you install your package to.*

`long_description` - *A more detailed description of your tool's function and use case. You may define a string variable to include it in setup() (see example above).*

`long_description_content_typ` - *Content type of the description.*
url - *Homepage of your program.*

`packages` - *Support in package finding in the working directory.*

`classifiers` - *Provide some classifers for this release. Review all supported options [here](https://pypi.org/classifiers/).*


# Useful links:
[Building a conda package from the scratch](https://docs.conda.io/projects/conda-build/en/latest/user-guide/tutorials/build-pkgs.html)

[Define a recipe](https://docs.conda.io/projects/conda-build/en/latest/concepts/recipe.html#conda-build-process)

[Options for composing the meta.yaml](https://docs.conda.io/projects/conda-build/en/latest/resources/define-metadata.html#build-script)

[Packaging python projects](https://packaging.python.org/tutorials/packaging-projects/)

[The `setup.py`](https://docs.python.org/3/distutils/setupscript.html)

[python classifier](https://pypi.org/classifiers/)

[license family for meta.yaml](https://github.com/conda/conda-build/blob/master/conda_build/license_family.py)

[Contributing to Anaconda](https://docs.anaconda.com/anaconda-cloud/user-guide/tasks/work-with-packages/#uploading-packages)

[Contributing to Bioconda](https://bioconda.github.io/contributor/index.html)
