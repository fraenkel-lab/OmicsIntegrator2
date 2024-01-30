from setuptools import setup

setup(
    name='OmicsIntegrator',
    packages=['OmicsIntegrator'],
    package_dir={'OmicsIntegrator': 'src'},
    package_data={'OmicsIntegrator': ['annotation/final_annotation.pickle']},
    version='2.4.1',
    url='https://github.com/fraenkel-lab/OmicsIntegrator2',
    python_requires='>=3.7',
    classifiers=[
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12'],
    license='MIT',
    author='alexlenail',
    author_email='alex@lenail.org',
    description='',
    install_requires=[
        "numpy>=1.26.0",
        "pandas>=2.2.0",
        "networkx>=3.2",
        "pcst_fast==1.0.10",
        "python-louvain",
        "goenrich==1.7.0",
        "scikit-learn==1.4",
        "axial>=0.2.2"
    ],
    entry_points={
        'console_scripts': [
            'OmicsIntegrator = src.__main__:main',
        ]
    },
)

