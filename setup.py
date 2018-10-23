import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="physlr",
    version="0.1.0",
    author="Shaun Jackman",
    author_email="sjackman@gmail.com",
    description="Physical Mapping of Linked Reads",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bcgsc/physlr",
    license="GPLv3",
    python_requires=">=3",
    install_requires=["networkx", "pygraphviz"],
    packages=["physlr"],
    scripts=["bin/physlr"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
)
