from setuptools import setup

setup(
    name='OmicsIntegrator',
    packages=['OmicsIntegrator'],
    package_dir={'OmicsIntegrator': 'src'},
    package_data={'OmicsIntegrator': ['annotation/final_annotation.pickle']},
    version='2.4.0',
    url='https://github.com/fraenkel-lab/OmicsIntegrator2',
    python_requires='>=3.8',
    classifiers=[
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9'],
    license='MIT',
    author='zfrenchee',
    author_email='alex@lenail.org',
    description='',
    install_requires=[
        "numpy",
        "pandas",
        "networkx>=3.1",
        "pcst_fast>=1.0.7",
        "goenrich",
        "scikit-learn",
        "axial",
        "scipy"
    ],
    entry_points={
        'console_scripts': [
            'OmicsIntegrator = OmicsIntegrator.__main__:main',
        ]
    },
)

