from distutils.core import setup

setup(
    name='OmicsIntegrator',
    version='0.2.0',
    packages=['OmicsIntegrator'],
    package_dir={'OmicsIntegrator2': 'src'},
    scripts=['src/forest.py'],
    url='https://github.com/fraenkel-lab/OmicsIntegrator2',
    classifiers=[
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'],
    license='GNU General Public License',
    author='zfrenchee',
    author_email='alex@lenail.org',
    description='',
    install_requires=[
        "numpy",
        "pandas",
        "garnet",
        "networkx",
        "pcst_fast"
    ],
)
