import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="galspec", # Replace with your own username
    version="0.2.4",
    author="Tom Bakx (edits: Stefanie Brackenhoff)",
    description="Creates sample galaxy spectra in the GHz/THz range",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/deshima-dev/GalSpec",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    package_data={"galspec": ["K17_Table7", "coeffBonato", "coeff_spinoglio", "COcoeff", "coeff_Rangwala"]},
)