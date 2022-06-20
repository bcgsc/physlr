import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="physlr",
    version="1.0.0",
    author="Amirhossein Afshinfard",
    author_email="aafshinfard@bcgsc.ca",
    description="Next-generation Physical Maps using Linked Reads",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bcgsc/physlr",
    license="GPLv3",
    python_requires=">=3",
    install_requires=["networkx", "numpy", "python-louvain", "scipy", "sklearn", "tqdm" ],
    packages=["physlr"],
    scripts=["bin/physlr"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
)
