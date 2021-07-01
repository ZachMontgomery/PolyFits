import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name = 'polyFits',
    version = '1.0.0',
    author = 'Zach Montgomery',
    author_email = 'zachary.s.montgomery@gmail.com',
    description = 'Multivariable Polynomial Curve Fitting',
    long_description = long_description,
    long_description_content_type = 'test/markdown',
    url = 'https://github.com/ZachMontgomery/MultivariablePolynomialFits',
    packages = ['polyFits'],
    classifiers = [
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        ],
    python_requires = '>=3.6',
    install_requires = ['numpy', 'ZachsModules', 'multiprocessing']
    )
