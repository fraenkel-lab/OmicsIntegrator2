from setuptools import setup

setup(
    name='OmicsIntegrator',
    packages=['OmicsIntegrator'],
    package_dir={'OmicsIntegrator': 'src'},
    version='0.2.8',
    url='https://github.com/fraenkel-lab/OmicsIntegrator2',
    classifiers=[
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6'],
    license='MIT',
    author='zfrenchee',
    author_email='alex@lenail.org',
    description='',
    install_requires=[
        "numpy",
        "pandas",
        "garnet",
        "networkx",
        "pcst_fast",
        "python-louvain",
        "py2cytoscape==0.6.0",
        "goenrich",
        "sklearn"
    ],
    entry_points={
        'console_scripts': [
            'OmicsIntegrator = src.__main__:main',
            'MergePrizeFiles = src.prizes:main'
        ]
    },
)




