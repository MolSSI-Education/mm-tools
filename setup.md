---
title: Setup
---

You can either install software on your computer, or you can use a website called [ChemCompute](https://chemcompute.org/) which will allow you to access computing resources and a pre-configured environment for this workshop without installing any software on your own computer. 
ChemCompute is only accessible to users with an academic email address. **Warning** If you are an academic user outside of the United States, there is a chance your email won't be properly recognized. If this happens, you may have to install software on your personal computer.

## Using a cloud computing environment

You can access a cloud computing environment with computational chemistry software installed by using [ChemCompute](https://chemcompute.org/).

1. [Register](https://chemcompute.org/register/) for a ChemCompute account. The easiest thing to do is to choose "Sign in with your University Login" on the right. You can login using your university email. You will need to confirm account registration from an email chemcompute will send you.
2. After confirming your account, you will be able to access a Jupyter notebook. The first time you start this workshop, you will need to get access to the files used in this workshop. Using the menu at the top, select "Jupyter" then "Clone a github repo to your notebook"
3. There will be a form on the page. Enter https://github.com/MolSSI-Education/mm-tools-workshop into the URL field and click "Pull Repo"
4. A Jupyter notebook environment with the files you need will start for you.

You only have to complete these steps once. Every other time you access ChemCompute, you will click "Jupyter" then "Use Psi4 or JupyterHub". Your files for this workshop will be in a folder called `mm-tools-workshop`.

## Using your personal computer

### Installing Python through Anaconda
[Python](https://python.org/) is a popular language for scientific computing, and great for general-purpose programming as well. Installing all of its scientific packages individually can be a bit difficult, however, so we recommend the all-in-one installer Anaconda.

1. Navigate to the [download page](https://www.anaconda.com/products/distribution) for Anaconda.
2. Double click the installer icon and follow the set-up instructions, keeping most of the default options. If you are Windows, make sure to choose to choose the option **Make Anaconda the default Python** during installation.

## Installing a Text Editor

You will need a text editor for this workshop. If you do not have a preferred text editor for writing code, we recommend [Visual Studio Code](https://code.visualstudio.com/). Download Visual Studio Code at the link and install on your computer.

## Installing Python Packages
We recommend using a conda environment for software installations. 

1.  Download or clone the [`mm-tools-workshop` repository](https://github.com/MolSSI-Education/mm-tools-workshop).

2. After you have Anaconda installed and have downloaded the file, you can install the software for this workshop using `conda`. Open a Terminal if you are on Mac or Linux. If you are on Windows, you should open a program called Anaconda Prompt. 

Navigate to where you have downloaded the environm,ment file. Type the following command into the terminal/anaconda prompt and press Enter. When prompted, answer 'yes' and wait for the install to finish:

~~~
conda env create -f environment.yml
~~~
{: .bash}

## Start a Jupyter notebook
Open the Anaconda Navigator, and choose your `mm-tools` environment. After switching to your environment, find the Jupyter notebook button on your Home tab. Click Launch.

It may take a few seconds to load the page, especially if it is the first time you have ever used the jupyter notebook, so don't panic if nothing loads for a few seconds.  Then a new window should open in your default internet browser. Use the file navigation window to navigate to your `mm-tools-workshop` folder.  In the upper right hand corner, click New, then choose Python 3 from the dropdown list.  You're ready to go!


{% include links.md %}
